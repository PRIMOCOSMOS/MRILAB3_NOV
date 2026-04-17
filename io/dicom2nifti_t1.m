function outFile = dicom2nifti_t1(dicomDir, outDir)
% dicom2nifti_t1 - 将 T1 结构像 DICOM 序列转换为 3D NIfTI
% 每个 DICOM 文件对应一个 2D 切片，按层位置排序后重组为 3D 体
% 完全 standalone，使用 MATLAB 内置 dicominfo/dicomread
%
% 输入:
%   dicomDir - T1 DICOM 文件目录
%   outDir   - 输出目录
%
% 输出:
%   outFile  - 输出 3D NIfTI 文件路径 (t1.nii)

fprintf('[dicom2nifti_t1] 开始处理 T1 DICOM 目录: %s\n', dicomDir);

% -------- 收集 DICOM 文件 --------
dcmFiles = get_dcm_files(dicomDir);
if isempty(dcmFiles)
    error('[dicom2nifti_t1] 未找到 DICOM 文件: %s', dicomDir);
end
fprintf('[dicom2nifti_t1] 找到 %d 个切片文件\n', numel(dcmFiles));

nSlices = numel(dcmFiles);

% -------- 读取所有切片头信息，并按层位置排序 --------
slicePos = zeros(1, nSlices);
for i = 1:nSlices
    try
        inf_i = dicominfo(dcmFiles{i});
        if isfield(inf_i, 'ImagePositionPatient') && numel(inf_i.ImagePositionPatient) >= 3
            % 取层法向量方向的投影（通常为 Z 分量或3个分量的合值）
            slicePos(i) = inf_i.ImagePositionPatient(3);  % Z 坐标
        elseif isfield(inf_i, 'InstanceNumber')
            slicePos(i) = inf_i.InstanceNumber;
        else
            slicePos(i) = i;
        end
    catch
        slicePos(i) = i;
    end
end
[~, sortIdx] = sort(slicePos);
dcmFiles = dcmFiles(sortIdx);

% -------- 读取第一个切片获取图像尺寸 --------
info1 = dicominfo(dcmFiles{1});
img1  = double(dicomread(dcmFiles{1}));
[ny, nx] = size(img1);  % DICOM row=Y, col=X

% -------- 逐层读取并堆叠成3D --------
vol3d = zeros(nx, ny, nSlices, 'single');
for s = 1:nSlices
    img = double(dicomread(dcmFiles{s}));
    vol3d(:,:,s) = single(img');  % 转置: X×Y 顺序
end

% -------- 提取空间参数 --------
pixSpacing = [1 1];
sliceThick = 1;
if isfield(info1, 'PixelSpacing') && numel(info1.PixelSpacing) >= 2
    pixSpacing = double(info1.PixelSpacing(:)');
end
if isfield(info1, 'SliceThickness')
    sliceThick = double(info1.SliceThickness);
end
dx = pixSpacing(2);
dy = pixSpacing(1);
dz = sliceThick;

% -------- 构建仿射矩阵 --------
affine = build_t1_affine(info1, dx, dy, dz, nx, ny, nSlices);

% -------- 写出 NIfTI --------
ensure_dir(outDir);
hdr = nifti_default_hdr([nx ny nSlices 1], [dx dy dz 0]);
hdr.affine = affine;
hdr.srow_x = affine(1,:);
hdr.srow_y = affine(2,:);
hdr.srow_z = affine(3,:);
hdr.pixdim = single([1, dx, dy, dz, 0, 0, 0, 0]);
hdr.descrip = sprintf('T1 3D nx=%d ny=%d nz=%d', nx, ny, nSlices);

outFile = fullfile(outDir, 't1.nii');
nifti_write(outFile, vol3d, hdr);
fprintf('[dicom2nifti_t1] 已写出: %s  维度=[%d %d %d]\n', outFile, nx, ny, nSlices);
end

% ======================================================================
% 内部辅助函数
% ======================================================================

function files = get_dcm_files(dir_path)
d = dir(dir_path);
files = {};
for i = 1:numel(d)
    if d(i).isdir, continue; end
    [~,~,ext] = fileparts(d(i).name);
    if strcmpi(ext,'.dcm') || isempty(ext)
        files{end+1} = fullfile(dir_path, d(i).name); %#ok<AGROW>
    end
end
if isempty(files)
    for i = 1:numel(d)
        if d(i).isdir && ~strcmp(d(i).name,'.') && ~strcmp(d(i).name,'..')
            sub = get_dcm_files(fullfile(dir_path, d(i).name));
            files = [files, sub]; %#ok<AGROW>
        end
    end
end
end

function affine = build_t1_affine(info, dx, dy, dz, nx, ny, nSlices)
F = eye(3);
if isfield(info, 'ImageOrientationPatient') && numel(info.ImageOrientationPatient)==6
    iop = double(info.ImageOrientationPatient(:));
    F(:,1) = iop(1:3);
    F(:,2) = iop(4:6);
    F(:,3) = cross(iop(1:3), iop(4:6));
end
T1 = [-nx/2*dx; -ny/2*dy; -nSlices/2*dz];
if isfield(info, 'ImagePositionPatient') && numel(info.ImagePositionPatient)>=3
    T1 = double(info.ImagePositionPatient(:));
end
affine = [F(:,1)*dx, F(:,2)*dy, F(:,3)*dz, T1; 0, 0, 0, 1];
end
