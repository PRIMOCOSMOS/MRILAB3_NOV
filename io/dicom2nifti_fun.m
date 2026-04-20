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
        vol3d = zeros(sliceSize, sliceSize, 1);
        patch = img2d(1:sliceSize, 1:sliceSize);
        % 与 MOSAIC 路径一致：先转置到 X×Y，再对 Y 轴翻转
        vol3d(:,:,1) = patch';
        vol3d(:,:,1) = vol3d(:, end:-1:1, 1);
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
affine = build_affine_from_dicom(info1, dx, dy, dz, sliceSize, sliceSize, nSlices, isMosaic);

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

    % 注意: DICOM 行=Y，列=X；NIfTI 习惯 X×Y，需转置。
    % 额外对 Y 轴翻转以与 dcm2niix/SPM 生成的体素排列一致。
    patch = mosaic(row_start:row_end, col_start:col_end);
    vol3d(:,:,s) = patch';
    vol3d(:,:,s) = vol3d(:, end:-1:1, s);
end
end

function affine = build_affine_from_dicom(info, dx, dy, dz, nx, ny, nSlices, isMosaic)
% 从 DICOM 方向余弦构建仿射矩阵（对齐 dcm2niix 的关键顺序）:
% 1) 在 LPS 中构建初始 Q44（row/col/slice + IPP）；
% 2) 若为 Siemens mosaic，先施加 mosaic 平移补偿；
% 3) LPS->RAS（翻转前两行）；
% 4) 与像素数据一致地做 Y 翻转后的头更新（nii_flipY 等价）。

rowDir = [1; 0; 0];
colDir = [0; 1; 0];
if isfield(info, 'ImageOrientationPatient') && numel(info.ImageOrientationPatient) == 6
    iop = double(info.ImageOrientationPatient(:));
    rowDir = iop(1:3);
    colDir = iop(4:6);
end
sliceDir = cross(rowDir, colDir);

ipp = [0; 0; 0];
if isfield(info, 'ImagePositionPatient') && numel(info.ImagePositionPatient) >= 3
    ipp = double(info.ImagePositionPatient(1:3));
end

% 先用 DICOM LPS 坐标构建未翻转的 Q44。
Q44 = [rowDir * dx, colDir * dy, sliceDir * dz, ipp;
       0, 0, 0, 1];

% Siemens mosaic：将起点从整幅 mosaic 左上角平移到有效切片网格中心。
if isMosaic
    nRowCol = ceil(sqrt(double(nSlices)));
    xdim = nx * nRowCol;
    ydim = ny * nRowCol;
    if isfield(info, 'Columns') && isfield(info, 'Rows')
        xdim = double(info.Columns);
        ydim = double(info.Rows);
    end
    lFactorX = (xdim - xdim / nRowCol) / 2.0;
    lFactorY = (ydim - ydim / nRowCol) / 2.0;
    Q44(1:3,4) = Q44(1:3,4) + Q44(1:3,1) * lFactorX + Q44(1:3,2) * lFactorY;
end

% LPS -> RAS
Q44(1:2,:) = -Q44(1:2,:);

% 数据写出时做了 Y 轴翻转，这里同步更新头（与 dcm2niix nii_flipY 一致）。
v = Q44 * [0; ny - 1; 0; 1];
R = Q44(1:3,1:3) * diag([1, -1, 1]);

affine = [R, v(1:3);
          0, 0, 0, 1];
end
