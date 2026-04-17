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
%   mniCfg   - MNI 目标网格配置:
%               mniCfg.dims    = [91 109 91]（2mm 各向同性）
%               mniCfg.voxSize = [2 2 2]（mm）
%               mniCfg.origin  = [46 64 37]（原点体素，1-based）
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
mni_dims   = mniCfg.dims(:)';
mni_vox    = mniCfg.voxSize(:)';
mni_origin = mniCfg.origin(:)';

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
    % 尺寸不匹配：位移场插值到源空间
    dx_mni = zeros(mx, my, mz);
    dy_mni = zeros(mx, my, mz);
    dz_mni = zeros(mx, my, mz);
    fprintf('[normalize_apply] 警告：位移场尺寸不匹配，跳过形变（使用线性映射）\n');
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
