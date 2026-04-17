function outFile = normalize_apply(inFile, flowFile, outDir, mniCfg)
% normalize_apply - 将 4D 功能像应用形变场，标准化到 MNI 空间
%
% 背景:
%   空间标准化（Normalize）将个体脑图像"扭曲"变形到公共参考空间（MNI），
%   使不同被试的脑区在空间上对齐，从而支持群体级统计分析。
%
% 方法:
%   对 MNI 目标网格中的每个体素，通过形变场（位移场）
%   找到其在个体空间中的对应位置，然后三线性插值采样。
%   （逆变换插值，避免采样空洞）
%
% 输入:
%   inFile   - 输入 4D NIfTI 文件（个体空间）
%   flowFile - 位移场 NIfTI 文件（由 dartel_warp 生成, [nx ny nz 3]）
%   outDir   - 输出目录
%   mniCfg   - 目标网格配置（两种格式均支持）:
%              A) DPABI风格:
%                 mniCfg.boundingBox = [xmin ymin zmin; xmax ymax zmax]（mm）
%                 mniCfg.voxSize     = [vx vy vz]（mm）
%              B) 传统网格:
%                 mniCfg.dims    = [nx ny nz]
%                 mniCfg.voxSize = [vx vy vz]
%                 mniCfg.origin  = [ox oy oz]（1-based）
%
% 输出:
%   outFile - 标准化后的 4D NIfTI 文件路径（前缀 'w'）

fprintf('[normalize_apply] 读取功能像: %s\n', inFile);
[data, hdr] = nifti_read(inFile);
[nx_in, ny_in, nz_in, nt] = size(data);

fprintf('[normalize_apply] 读取位移场: %s\n', flowFile);
[flowData, ~] = nifti_read(flowFile);
% flowData: [nx ny nz 3]，单位为个体空间体素位移

% 位移场尺寸应与输入图像匹配
[nx_f, ny_f, nz_f, nd] = size(flowData);
if nd ~= 3
    error('[normalize_apply] 位移场第4维应为3（x/y/z方向）');
end

% -------- MNI 目标网格参数 --------
[mni_dims, mni_vox, mni_origin] = parse_target_grid(mniCfg);

mx = mni_dims(1); my = mni_dims(2); mz = mni_dims(3);

% MNI 体素 → 世界坐标（mm）
% MNI 仿射：以 origin 为零点
mni_affine = [mni_vox(1) 0 0 -mni_vox(1)*(mni_origin(1)-1);
              0 mni_vox(2) 0 -mni_vox(2)*(mni_origin(2)-1);
              0 0 mni_vox(3) -mni_vox(3)*(mni_origin(3)-1);
              0 0 0           1                             ];

% 输入图像仿射（个体空间世界坐标 → 体素坐标）
src_affine_inv = inv(hdr.affine);

% -------- 构建 MNI 网格 --------
[Xm, Ym, Zm] = ndgrid(1:mx, 1:my, 1:mz);
mni_coords = [Xm(:)'; Ym(:)'; Zm(:)'];  % [3 × nMNI]

% MNI 体素 → 世界坐标
mni0 = [mni_coords - 1; ones(1, mx*my*mz)];
mni_world = mni_affine * mni0;  % [4 × nMNI]

% 世界坐标 → 输入图像体素坐标（近似：忽略形变，先找个体空间中的体素）
src_vox = src_affine_inv * mni_world;  % [4 × nMNI]
src_vox = src_vox(1:3,:) + 1;  % 0-based → 1-based

% -------- 应用形变场（从 MNI 到个体空间的逆变形）--------
% 注意：flowData 存储的是 "个体 → MNI" 的位移，
%       需要通过迭代求逆（简化：直接取负值近似逆变形）
% 实际精确实现应使用 Powell's method 或固定点迭代求逆
% 此处采用简化的近似逆：将位移场在 MNI 网格中插值并取反
if nx_f == nx_in && ny_f == ny_in && nz_f == nz_in
    % 位移场与输入图像同尺寸：直接使用
    dx_mni = reshape(trilinear_interp(flowData(:,:,:,1), src_vox), mx, my, mz);
    dy_mni = reshape(trilinear_interp(flowData(:,:,:,2), src_vox), mx, my, mz);
    dz_mni = reshape(trilinear_interp(flowData(:,:,:,3), src_vox), mx, my, mz);
else
    % 尺寸不匹配：对位移场进行三线性重采样到与输入图像匹配的网格，再插值到 MNI 网格
    % 构建位移场体素坐标（目标：在 src_vox 处采样位移场）
    % 位移场尺寸为 [nx_f ny_f nz_f 3]，需将 src_vox 坐标缩放到位移场坐标系
    scale_x = nx_f / nx_in;  scale_y = ny_f / ny_in;  scale_z = nz_f / nz_in;
    src_vox_f = [src_vox(1,:) * scale_x; src_vox(2,:) * scale_y; src_vox(3,:) * scale_z];
    dx_mni = reshape(trilinear_interp(flowData(:,:,:,1), src_vox_f), mx, my, mz);
    dy_mni = reshape(trilinear_interp(flowData(:,:,:,2), src_vox_f), mx, my, mz);
    dz_mni = reshape(trilinear_interp(flowData(:,:,:,3), src_vox_f), mx, my, mz);
    fprintf('[normalize_apply] 位移场已重采样（%d×%d×%d → %d×%d×%d）\n', ...
        nx_f, ny_f, nz_f, nx_in, ny_in, nz_in);
end

% 最终采样坐标（个体空间体素坐标）
X_src = reshape(src_vox(1,:), mx, my, mz) - dx_mni;
Y_src = reshape(src_vox(2,:), mx, my, mz) - dy_mni;
Z_src = reshape(src_vox(3,:), mx, my, mz) - dz_mni;
sampCoords = [X_src(:)'; Y_src(:)'; Z_src(:)'];

% -------- 逐时间点插值 --------
fprintf('[normalize_apply] 开始对 %d 个时间点进行空间标准化...\n', nt);
dataOut = zeros(mx, my, mz, nt, 'single');

for t = 1:nt
    vol_t = double(data(:,:,:,t));
    vals  = trilinear_interp(vol_t, sampCoords);
    dataOut(:,:,:,t) = single(reshape(vals, mx, my, mz));

    if mod(t,50)==0
        fprintf('[normalize_apply]  已处理 %d/%d 时间点\n', t, nt);
    end
end

% -------- 构建输出 NIfTI 头 --------
hdrOut = nifti_default_hdr([mx my mz nt], [mni_vox(:)' hdr.pixdim(5)]);
hdrOut.affine = mni_affine;
hdrOut.srow_x = mni_affine(1,:);
hdrOut.srow_y = mni_affine(2,:);
hdrOut.srow_z = mni_affine(3,:);
hdrOut.descrip = 'NormalizedToMNI';

% -------- 写出 --------
ensure_dir(outDir);
[~, fname, ext] = fileparts(inFile);
outFile = fullfile(outDir, ['w' fname ext]);
nifti_write(outFile, dataOut, hdrOut);
fprintf('[normalize_apply] 已写出: %s  [%d %d %d %d]\n', outFile, mx, my, mz, nt);
end

function [dims, vox, origin] = parse_target_grid(cfgIn)
% 兼容两类标准化目标参数格式
if isfield(cfgIn, 'boundingBox')
    bbox = cfgIn.boundingBox;
    if ~isequal(size(bbox), [2 3])
        error('[normalize_apply] boundingBox 必须为 2x3 矩阵');
    end
    if ~isfield(cfgIn, 'voxSize')
        error('[normalize_apply] 使用 boundingBox 时必须提供 voxSize');
    end
    vox = cfgIn.voxSize(:)';
    if numel(vox) ~= 3 || any(vox <= 0)
        error('[normalize_apply] voxSize 必须为 1x3 正数向量');
    end
    bboxMin = bbox(1,:);
    bboxMax = bbox(2,:);
    dims = round((bboxMax - bboxMin) ./ vox) + 1;
    % 将 mm 坐标的 bbox 最小角转换为 1-based 体素原点：
    % 令体素坐标 i=1 时对应 bboxMin，可得
    % bboxMin_i = vox_i * (1 - origin_i)  =>  origin_i = 1 - bboxMin_i/vox_i
    origin = 1 - bboxMin ./ vox;
else
    if ~isfield(cfgIn, 'dims') || ~isfield(cfgIn, 'voxSize') || ~isfield(cfgIn, 'origin')
        error('[normalize_apply] 目标网格配置缺少字段，需提供 boundingBox+voxSize 或 dims+voxSize+origin');
    end
    dims = cfgIn.dims(:)';
    vox = cfgIn.voxSize(:)';
    origin = cfgIn.origin(:)';
end
end
