function outFile = dicom2nifti_t1(dicomDir, outDir, cfg)
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

if nargin < 3
    cfg = struct();
end

if use_dpabi_dicom_backend(cfg)
    outFile = dicom2nifti_t1_dpabi(dicomDir, outDir, cfg);
    return;
end

fprintf('[dicom2nifti_t1] 开始处理 T1 DICOM 目录: %s\n', dicomDir);

% -------- 收集 DICOM 文件 --------
dcmFiles = get_dcm_files(dicomDir);
if isempty(dcmFiles)
    error('[dicom2nifti_t1] 未找到 DICOM 文件: %s', dicomDir);
end
fprintf('[dicom2nifti_t1] 找到 %d 个切片文件\n', numel(dcmFiles));

nSlices = numel(dcmFiles);

% -------- 读取所有切片头信息，并按层位置排序 --------
refInfo = dicominfo(dcmFiles{1});
rowDir = [1; 0; 0];
colDir = [0; 1; 0];
if isfield(refInfo, 'ImageOrientationPatient') && numel(refInfo.ImageOrientationPatient) == 6
    iop = double(refInfo.ImageOrientationPatient(:));
    rowDir = iop(1:3);
    colDir = iop(4:6);
end
sliceNormal = cross(rowDir, colDir);
if norm(sliceNormal) > 0
    sliceNormal = sliceNormal / norm(sliceNormal);
end

slicePos = zeros(1, nSlices);
for i = 1:nSlices
    try
        inf_i = dicominfo(dcmFiles{i});
        if isfield(inf_i, 'ImagePositionPatient') && numel(inf_i.ImagePositionPatient) >= 3
            ipp = double(inf_i.ImagePositionPatient(:));
            % 采用沿层法向量的投影进行排序，避免仅按 Z 坐标导致斜切数据乱序。
            if norm(sliceNormal) > 0
                slicePos(i) = dot(ipp, sliceNormal);
            else
                slicePos(i) = ipp(3);
            end
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
    vol3d(:,:,s) = single(img');
    vol3d(:,:,s) = vol3d(:, end:-1:1, s);  % 与功能像一致的 Y 轴翻转
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

% -------- 重排到近 RAS 轴顺序（与 dcm2niix -x 行为一致）--------
affineRaw = build_t1_affine(info1, dx, dy, dz, ny);
[vol3d, Mcanon] = canonical_reorient_3d(vol3d, affineRaw);
affine = affineRaw * Mcanon;

[nx, ny, nSlices] = size(vol3d);

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

function affine = build_t1_affine(info, dx, dy, dz, ny)
rowDir = [1; 0; 0];
colDir = [0; 1; 0];
if isfield(info, 'ImageOrientationPatient') && numel(info.ImageOrientationPatient)==6
    iop = double(info.ImageOrientationPatient(:));
    rowDir = iop(1:3);
    colDir = iop(4:6);
end
sliceDir = cross(rowDir, colDir);

ipp = [0; 0; 0];
if isfield(info, 'ImagePositionPatient') && numel(info.ImagePositionPatient) >= 3
    ipp = double(info.ImagePositionPatient(1:3));
end

Q44 = [rowDir * dx, colDir * dy, sliceDir * dz, ipp;
       0, 0, 0, 1];

% LPS -> RAS
Q44(1:2,:) = -Q44(1:2,:);

% 与数据写出中的 Y 轴翻转保持一致（nii_flipY 头更新等价）。
v = Q44 * [0; ny - 1; 0; 1];
R = Q44(1:3,1:3) * diag([1, -1, 1]);

affine = [R, v(1:3); 0, 0, 0, 1];
end

function [volOut, M] = canonical_reorient_3d(volIn, affineIn)
% 将 3D 体数据重排到最接近 RAS 轴的体素顺序，同时返回对应的 0-based 坐标变换 M。
R = affineIn(1:3,1:3);
sz = size(volIn);

permList = perms(1:3);
bestScore = -inf;
bestPerm = [1 2 3];
bestSign = [1 1 1];

for pi = 1:size(permList,1)
    p = permList(pi,:);
    for mask = 0:7
        s = [1 1 1];
        if bitget(mask,1), s(1) = -1; end
        if bitget(mask,2), s(2) = -1; end
        if bitget(mask,3), s(3) = -1; end

        Rt = R(:, p) * diag(s);
        score = sum(abs(diag(Rt)));
        if score > bestScore
            bestScore = score;
            bestPerm = p;
            bestSign = s;
        end
    end
end

volOut = permute(volIn, bestPerm);
for k = 1:3
    if bestSign(k) < 0
        volOut = flip(volOut, k);
    end
end

M = eye(4);
M(1:3,:) = 0;
for oldAxis = 1:3
    newAxis = find(bestPerm == oldAxis, 1, 'first');
    M(oldAxis, newAxis) = bestSign(newAxis);
    if bestSign(newAxis) < 0
        M(oldAxis, 4) = sz(oldAxis) - 1;
    end
end
end

    function tf = use_dpabi_dicom_backend(cfg)
    tf = isstruct(cfg) && isfield(cfg, 'dpabi') && ...
        isfield(cfg.dpabi, 'useDicomConvert') && logical(cfg.dpabi.useDicomConvert);
    end
