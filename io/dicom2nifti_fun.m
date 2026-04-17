function outFile = dicom2nifti_fun(dicomDir, outDir, cfg)
% dicom2nifti_fun - 将 EPI BOLD DICOM 文件转换为 4D NIfTI
% 支持 Siemens MOSAIC 格式（多层拼接于单幅2D图像中）
% 完全 standalone，使用 MATLAB 内置 dicominfo/dicomread
%
% 输入:
%   dicomDir - DICOM 文件所在目录（字符串）
%   outDir   - 输出 NIfTI 目录（字符串）
%   cfg      - 配置结构体，需包含:
%               cfg.nSlices   - EPI 层数（如 36）
%               cfg.sliceSize - 每层体素大小（如 64，即 64×64）
%               cfg.TR        - 重复时间（秒）
%
% 输出:
%   outFile  - 输出 4D NIfTI 文件路径
%
% 原理:
%   MOSAIC 图像将 nSlices 层的 sliceSize×sliceSize 图像
%   拼接成一幅大图（行×列各有 ceil(sqrt(nSlices)) 个子块）。
%   本函数自动检测并解码 MOSAIC。

fprintf('[dicom2nifti_fun] 开始处理 EPI DICOM 目录: %s\n', dicomDir);

% -------- 收集所有 DICOM 文件 --------
dcmFiles = get_dcm_files(dicomDir);
if isempty(dcmFiles)
    error('[dicom2nifti_fun] 未找到 DICOM 文件: %s', dicomDir);
end
fprintf('[dicom2nifti_fun] 找到 %d 个 DICOM 文件\n', numel(dcmFiles));

% -------- 读取第一个文件判断是否为 MOSAIC --------
info1 = dicominfo(dcmFiles{1});
isMosaic = detect_mosaic(info1, cfg.nSlices, cfg.sliceSize);
fprintf('[dicom2nifti_fun] MOSAIC 格式: %s\n', mat2str(isMosaic));

% -------- 按 InstanceNumber 排序（时间维）--------
nFiles = numel(dcmFiles);
instNums = zeros(1, nFiles);
for i = 1:nFiles
    try
        tmpInfo = dicominfo(dcmFiles{i});
        if isfield(tmpInfo, 'InstanceNumber')
            instNums(i) = tmpInfo.InstanceNumber;
        else
            instNums(i) = i;
        end
    catch
        instNums(i) = i;
    end
end
[~, sortIdx] = sort(instNums);
dcmFiles = dcmFiles(sortIdx);

nT = nFiles;  % 时间点数（每个 DICOM 对应一个 TR/Volume）

% -------- 预分配输出数组 --------
nSlices   = cfg.nSlices;
sliceSize = cfg.sliceSize;
vol4D     = zeros(sliceSize, sliceSize, nSlices, nT, 'single');

% -------- 逐文件读取 --------
for t = 1:nT
    img2d = double(dicomread(dcmFiles{t}));
    if isMosaic
        vol3d = demosaic_epi(img2d, nSlices, sliceSize);
    else
        % 非 MOSAIC：每个文件直接是 sliceSize×sliceSize 的2D图
        % 需要按层序号累积（此路径适用于每层一个DICOM文件的情况）
        vol3d = reshape(img2d, sliceSize, sliceSize, 1);
    end
    vol4D(:,:,:,t) = single(vol3d);
end

% -------- 构建仿射矩阵 --------
% 从 DICOM 头提取像素间距和层厚
pixSpacing = [1 1];
sliceThick = 3;
if isfield(info1, 'PixelSpacing') && numel(info1.PixelSpacing) >= 2
    pixSpacing = double(info1.PixelSpacing(:)');
end
if isfield(info1, 'SliceThickness')
    sliceThick = double(info1.SliceThickness);
end

dx = pixSpacing(2);  % 列方向 = x
dy = pixSpacing(1);  % 行方向 = y
dz = sliceThick;
dt = cfg.TR;

% 图像方向：默认 LAS（左-前-上），与 DICOM 方向余弦推导
affine = build_affine_from_dicom(info1, dx, dy, dz, sliceSize, nSlices);

% -------- 构建 NIfTI 头 --------
hdr = nifti_default_hdr([sliceSize sliceSize nSlices nT], [dx dy dz dt]);
hdr.affine = affine;
hdr.srow_x = affine(1,:);
hdr.srow_y = affine(2,:);
hdr.srow_z = affine(3,:);
hdr.pixdim = single([1, dx, dy, dz, dt, 0, 0, 0]);
hdr.xyzt_units = 10;  % mm + sec
hdr.descrip = sprintf('EPI BOLD nSlices=%d sliceSize=%d TR=%.2f', nSlices, sliceSize, dt);

% -------- 写出 NIfTI --------
ensure_dir(outDir);
outFile = fullfile(outDir, 'bold_4d.nii');
nifti_write(outFile, vol4D, hdr);
fprintf('[dicom2nifti_fun] 已写出: %s  维度=[%d %d %d %d]\n', ...
    outFile, sliceSize, sliceSize, nSlices, nT);
end

% ======================================================================
% 内部辅助函数
% ======================================================================

function files = get_dcm_files(dir_path)
% 递归收集目录中的所有 DICOM 文件（.dcm 或无扩展名）
d = dir(dir_path);
files = {};
for i = 1:numel(d)
    if d(i).isdir, continue; end
    fname = d(i).name;
    [~,~,ext] = fileparts(fname);
    if strcmpi(ext, '.dcm') || isempty(ext)
        files{end+1} = fullfile(dir_path, fname); %#ok<AGROW>
    end
end
% 若无，尝试递归子目录
if isempty(files)
    for i = 1:numel(d)
        if d(i).isdir && ~strcmp(d(i).name,'.') && ~strcmp(d(i).name,'..')
            sub = get_dcm_files(fullfile(dir_path, d(i).name));
            files = [files, sub]; %#ok<AGROW>
        end
    end
end
end

function flag = detect_mosaic(info, nSlices, sliceSize)
% 判断 DICOM 是否为 MOSAIC 格式
% 判据: 图像行数 > sliceSize 且约为 ceil(sqrt(nSlices))*sliceSize
flag = false;
if ~isfield(info, 'Rows'), return; end
rows = info.Rows;
tilesPerRow = ceil(sqrt(nSlices));
expectedSize = tilesPerRow * sliceSize;
if rows >= expectedSize * 0.9  % 90% 容差
    flag = true;
end
% 也可检查 CSA 头信息（Siemens 特有）
if isfield(info, 'Private_0019_100a')
    flag = true;
end
end

function vol3d = demosaic_epi(mosaic, nSlices, sliceSize)
% 将 MOSAIC 2D 图像解码为 3D 体数据
% mosaic   - [M×M] MOSAIC 大图（M = tilesPerRow * sliceSize）
% nSlices  - 实际层数
% sliceSize - 每层体素尺寸（sliceSize × sliceSize）
%
% Siemens MOSAIC 排列规则:
%   从左上角开始，按行优先排列，每个子块对应一层
%   层序号: [1  2  3  ...
%             n  n+1 ...]

tilesPerRow = ceil(sqrt(nSlices));
vol3d = zeros(sliceSize, sliceSize, nSlices);

for s = 1:nSlices
    row_tile = floor((s-1) / tilesPerRow);  % 0-based 行块
    col_tile = mod( (s-1),  tilesPerRow);   % 0-based 列块

    row_start = row_tile * sliceSize + 1;
    row_end   = row_start + sliceSize - 1;
    col_start = col_tile * sliceSize + 1;
    col_end   = col_start + sliceSize - 1;

    % 防越界
    if row_end > size(mosaic,1) || col_end > size(mosaic,2)
        break;
    end

    % 注意: DICOM 行=Y，列=X；NIfTI 习惯 X×Y，需转置
    patch = mosaic(row_start:row_end, col_start:col_end);
    vol3d(:,:,s) = patch';  % 转置使 x=列, y=行
end
end

function affine = build_affine_from_dicom(info, dx, dy, dz, sliceSize, nSlices)
% 从 DICOM 头的 ImageOrientationPatient 和 ImagePositionPatient 构建仿射矩阵
% 若头信息缺失则使用默认值

F = eye(3);  % 方向余弦矩阵
if isfield(info, 'ImageOrientationPatient') && numel(info.ImageOrientationPatient) == 6
    iop = double(info.ImageOrientationPatient(:));
    F(:,1) = iop(1:3);  % 列方向（X体素步进）
    F(:,2) = iop(4:6);  % 行方向（Y体素步进）
    F(:,3) = cross(iop(1:3), iop(4:6));  % 层方向
end

% 第一层位置
T1 = [-sliceSize/2*dx; -sliceSize/2*dy; -nSlices/2*dz];
if isfield(info, 'ImagePositionPatient') && numel(info.ImagePositionPatient) >= 3
    T1 = double(info.ImagePositionPatient(:));
end

affine = [F(:,1)*dx, F(:,2)*dy, F(:,3)*dz, T1;
          0,         0,         0,          1 ];
end
