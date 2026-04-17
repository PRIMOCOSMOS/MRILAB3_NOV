function outFile = smooth_3d(inFile, outDir, fwhm_mm)
% smooth_3d - 3D 高斯核空间平滑
%
% 背景:
%   空间平滑用于:
%   1. 提高信噪比（SNR）
%   2. 补偿个体间解剖差异（配准后残余误差）
%   3. 满足高斯随机场（GRF）理论的假设（用于多重比较校正）
%
% 实现原理:
%   3D 高斯平滑通过 3 次 1D 卷积（可分离性）高效实现:
%   G_3d(x,y,z) = G_1d(x) × G_1d(y) × G_1d(z)
%   FWHM 与标准差 σ 的关系: σ = FWHM / (2√(2ln2)) ≈ FWHM / 2.3548
%
% 输入:
%   inFile   - 输入 4D（或3D）NIfTI 文件路径
%   outDir   - 输出目录
%   fwhm_mm  - 高斯核 FWHM（mm），标量或 [fx fy fz] 向量
%
% 输出:
%   outFile - 平滑后的 NIfTI 文件路径（前缀 's'）

fprintf('[smooth_3d] 读取: %s\n', inFile);
[data, hdr] = nifti_read(inFile);
[nx, ny, nz, nt] = size(data);

% 体素尺寸（mm）
voxSz = sqrt(sum(hdr.affine(1:3,1:3).^2, 1));
dx = voxSz(1); dy = voxSz(2); dz = voxSz(3);

% FWHM 向量
if isscalar(fwhm_mm)
    fwhm_mm = fwhm_mm * [1 1 1];
end
fx = fwhm_mm(1); fy = fwhm_mm(2); fz = fwhm_mm(3);
fprintf('[smooth_3d] FWHM=[%.1f %.1f %.1f] mm, 体素尺寸=[%.2f %.2f %.2f] mm\n', ...
    fx, fy, fz, dx, dy, dz);

% σ（体素单位）
sig_x = (fx/dx) / (2*sqrt(2*log(2)));  % FWHM → σ（体素）
sig_y = (fy/dy) / (2*sqrt(2*log(2)));
sig_z = (fz/dz) / (2*sqrt(2*log(2)));

% -------- 构建 1D 高斯核 --------
kx = make_gauss_kernel(sig_x);
ky = make_gauss_kernel(sig_y);
kz = make_gauss_kernel(sig_z);

fprintf('[smooth_3d] 卷积核大小: [%d %d %d]\n', numel(kx), numel(ky), numel(kz));

% -------- 对每个时间点平滑 --------
dataOut = single(data);
for t = 1:nt
    vol = double(data(:,:,:,t));

    % 可分离 3D 卷积（3次1D卷积）
    vol = conv_axis(vol, kx, 1);  % 沿 X 轴
    vol = conv_axis(vol, ky, 2);  % 沿 Y 轴
    vol = conv_axis(vol, kz, 3);  % 沿 Z 轴

    dataOut(:,:,:,t) = single(vol);

    if mod(t,50)==0
        fprintf('[smooth_3d]  已平滑 %d/%d 时间点\n', t, nt);
    end
end

% -------- 写出 --------
hdr.descrip = sprintf('Smoothed_FWHM=[%.1f %.1f %.1f]mm', fx, fy, fz);
ensure_dir(outDir);
[~, fname, ext] = fileparts(inFile);
outFile = fullfile(outDir, ['s' fname ext]);
nifti_write(outFile, dataOut, hdr);
fprintf('[smooth_3d] 已写出: %s\n', outFile);
end

% ======================================================================
% 内部辅助函数
% ======================================================================

function k = make_gauss_kernel(sigma)
% 构建归一化 1D 高斯核
% 核半径为 3σ（截断），确保长度为奇数
if sigma <= 0
    k = 1;
    return;
end
r = ceil(3 * sigma);
x = -r:r;
k = exp(-x.^2 / (2*sigma^2));
k = k / sum(k);  % 归一化
end

function vol_out = conv_axis(vol, kernel, dim)
% 沿指定维度对3D数据进行1D卷积
% 使用边界反射填充（Replicate 边界条件）
sz = size(vol);
n  = sz(dim);
klen = numel(kernel);
pad  = (klen - 1) / 2;  % 每侧填充量（klen为奇数）
pad  = floor(pad);

% 重排，使目标维度在第1维
order = 1:3;
order(dim) = 1; order(1) = dim;
vol_p = permute(vol, order);  % [n, ...]
sp = size(vol_p);

% 展平后两维
vol_mat = reshape(vol_p, n, []);  % [n × nOther]

% 边界填充（Replicate: 用边界值复制）
vol_pad = [repmat(vol_mat(1,:), pad, 1); vol_mat; repmat(vol_mat(end,:), pad, 1)];

% 卷积（逐列）
nOther  = size(vol_mat, 2);
out_mat = zeros(n, nOther);
for c = 1:nOther
    out_mat(:, c) = conv(vol_pad(:,c), kernel(:), 'valid');
end

% 恢复形状
vol_out = ipermute(reshape(out_mat, sp), order);
end
