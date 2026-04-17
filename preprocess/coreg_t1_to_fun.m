function [outFile, finalParams] = coreg_t1_to_fun(t1File, funFile, outDir)
% coreg_t1_to_fun - 将 T1 结构像刚体配准到功能像（EPI）空间
%
% 背景:
%   T1 和 EPI 来自同一次扫描，但坐标系可能略有偏差。
%   通过刚体配准（6自由度）将 T1 对齐到 EPI，
%   以便后续分割结果能直接对应到功能像空间。
%
% 方法:
%   互信息（Mutual Information, MI）作为相似性度量，
%   Gauss-Newton 方法优化刚体参数。
%   注意：只修改 T1 的 NIfTI 仿射矩阵，不重采样（等同于 SPM coreg estimate）。
%
% 输入:
%   t1File  - T1 NIfTI 文件路径（移动图像）
%   funFile - EPI 功能像 NIfTI 文件路径（参考图像，取第1个时间点）
%   outDir  - 输出目录
%
% 输出:
%   outFile     - 仿射更新后的 T1 NIfTI（前缀 'coreg_'）
%   finalParams - 最终刚体参数 [tx ty tz rx ry rz]

fprintf('[coreg_t1_to_fun] T1=%s\n', t1File);
fprintf('[coreg_t1_to_fun] Fun=%s\n', funFile);

% -------- 读取图像 --------
[t1Data, t1Hdr] = nifti_read(t1File);
[funData, funHdr] = nifti_read(funFile);

% EPI 参考：取第1帧
if ndims(funData) == 4
    funRef = double(funData(:,:,:,1));
else
    funRef = double(funData);
end
t1Vol = double(t1Data(:,:,:,1));

% -------- 归一化强度（提高 MI 数值稳定性）--------
funRef = normalize_intensity(funRef, 64);
t1Vol  = normalize_intensity(t1Vol,  64);

% -------- Gauss-Newton 迭代（互信息最大化）--------
p = zeros(1,6);  % 初始参数
maxIter = 64;
tol     = 1e-6;

% 功能像体素坐标网格（下采样用于快速配准）
sep = 4;  % 采样间距（体素）
[nx, ny, nz] = size(funRef);
voxSzFun = sqrt(sum(funHdr.affine(1:3,1:3).^2, 1));
sep_vox  = max(1, round(sep ./ voxSzFun));

xs = 1:sep_vox(1):nx;
ys = 1:sep_vox(2):ny;
zs = 1:sep_vox(3):nz;
[Xf, Yf, Zf] = ndgrid(xs, ys, zs);
nPts = numel(Xf);

% 功能像世界坐标
funCoords_world = vox2world(funHdr.affine, [Xf(:)'; Yf(:)'; Zf(:)']);

for iter = 1:maxIter
    % 将功能像世界坐标变换到当前 T1 体素坐标
    % T1 当前仿射 = 原始仿射 × 参数变换矩阵（的逆）
    % （将 p 理解为 T1 相对于 Fun 的偏移）
    M  = rigid_mat(p);
    t1AffineNew = t1Hdr.affine * inv(M);  % 更新后的T1仿射

    t1Coords_vox = world2vox(t1AffineNew, funCoords_world);

    % 在 T1 中插值
    t1Vals = trilinear_interp(t1Vol, t1Coords_vox);

    % 参考像（功能像）在采样点的值
    funVals = trilinear_interp(funRef, [Xf(:)'; Yf(:)'; Zf(:)']);

    % 计算互信息梯度（使用有限差分近似）
    dMI_dp = zeros(1,6);
    eps_p   = [1e-3, 1e-3, 1e-3, 1e-5, 1e-5, 1e-5];
    mi0 = mutual_info(funVals, t1Vals);

    for k = 1:6
        dp_k    = zeros(1,6);
        dp_k(k) = eps_p(k);
        M_k     = rigid_mat(p + dp_k);
        t1Aff_k = t1Hdr.affine * inv(M_k);
        t1Vals_k = trilinear_interp(t1Vol, world2vox(t1Aff_k, funCoords_world));
        mi_k     = mutual_info(funVals, t1Vals_k);
        dMI_dp(k) = (mi_k - mi0) / eps_p(k);
    end

    % 梯度上升（最大化互信息）
    stepSize = 0.1;
    dp = stepSize * dMI_dp;
    p  = p + dp;

    if norm(dp) < tol, break; end
    if mod(iter,16)==0
        fprintf('[coreg_t1_to_fun]  iter=%d MI=%.4f dp_norm=%.2e\n', ...
            iter, mi0, norm(dp));
    end
end

finalParams = p;
fprintf('[coreg_t1_to_fun] 配准完成: [tx=%.2f ty=%.2f tz=%.2f rx=%.4f ry=%.4f rz=%.4f]\n', ...
    p(1),p(2),p(3),p(4),p(5),p(6));

% -------- 更新 T1 仿射矩阵 --------
M = rigid_mat(p);
t1HdrNew = t1Hdr;
t1HdrNew.affine = t1Hdr.affine * inv(M);
t1HdrNew.srow_x = t1HdrNew.affine(1,:);
t1HdrNew.srow_y = t1HdrNew.affine(2,:);
t1HdrNew.srow_z = t1HdrNew.affine(3,:);
t1HdrNew.descrip = sprintf('CoregedToFun p=[%.2f %.2f %.2f %.4f %.4f %.4f]', p);

% -------- 写出 --------
ensure_dir(outDir);
[~, fname, ext] = fileparts(t1File);
outFile = fullfile(outDir, ['coreg_' fname ext]);
nifti_write(outFile, single(t1Data), t1HdrNew);
fprintf('[coreg_t1_to_fun] 已写出: %s\n', outFile);
end

% ======================================================================
% 内部辅助函数
% ======================================================================

function worldCoords = vox2world(affine, voxCoords)
% 体素坐标（1-based）→ 世界坐标（mm）
vox0 = [voxCoords - 1; ones(1, size(voxCoords,2))];  % 0-based齐次
w = affine * vox0;
worldCoords = w(1:3,:);
end

function voxCoords = world2vox(affine, worldCoords)
% 世界坐标（mm）→ 体素坐标（1-based）
w = [worldCoords; ones(1, size(worldCoords,2))];
v = affine \ w;
voxCoords = v(1:3,:) + 1;  % 0-based → 1-based
end

function mi = mutual_info(a, b, nBins)
% 计算两个向量的互信息（基于联合直方图）
if nargin < 3, nBins = 32; end
a = a(:); b = b(:);
% 归一化到 [1, nBins]
a = min(max(round((a - min(a)) / (max(a)-min(a)+eps) * (nBins-1)) + 1, 1), nBins);
b = min(max(round((b - min(b)) / (max(b)-min(b)+eps) * (nBins-1)) + 1, 1), nBins);

% 联合直方图
joint = accumarray([a, b], 1, [nBins, nBins]);
joint = joint / sum(joint(:));

pa = sum(joint, 2);
pb = sum(joint, 1)';

% 互信息: MI = H(a) + H(b) - H(a,b)
ha = -sum(pa(pa>0) .* log(pa(pa>0)));
hb = -sum(pb(pb>0) .* log(pb(pb>0)));
hj = -sum(joint(joint>0) .* log(joint(joint>0)));
mi = ha + hb - hj;
end

function vol_norm = normalize_intensity(vol, nLevels)
% 将体数据强度归一化到 [0, nLevels-1] 的整数
vol_norm = vol(:);
lo = prctile(vol_norm, 1);
hi = prctile(vol_norm, 99);
vol_norm = (vol - lo) / (hi - lo + eps);
vol_norm = round(min(max(vol_norm, 0), 1) * (nLevels-1));
vol_norm = reshape(vol_norm, size(vol));
end
