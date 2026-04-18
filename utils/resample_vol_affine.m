function V_tgt = resample_vol_affine(V_src, src_affine, tgt_affine, tgt_dims)
% resample_vol_affine - 使用仿射矩阵将源体数据重采样到目标坐标空间
% 完全 standalone，不依赖 SPM 或其他工具箱
%
% 原理:
%   对目标空间的每个体素 t（1-based）:
%     1. 目标世界坐标: w = tgt_affine  * [t-1; 1]   （0-based 体素索引送入仿射）
%     2. 源体素坐标:   s = src_affine_inv * w + 1     （0-based → 1-based）
%     3. 三线性插值:   V_tgt(t) = trilinear_interp(V_src, s)
%
% 此函数正确处理两个图像在不同分辨率、不同视野（FOV）、不同朝向时的
% 空间对齐问题。超出源图像范围的体素赋值为 0。
%
% 输入:
%   V_src      - 源 3D 数组
%   src_affine - 源图像 4×4 仿射矩阵（NIfTI srow，体素0-based → 世界坐标 mm）
%   tgt_affine - 目标空间 4×4 仿射矩阵（同上）
%   tgt_dims   - 目标空间维度 [tx ty tz]（整数向量）
%
% 输出:
%   V_tgt - 重采样后的 3D 数组，尺寸为 tgt_dims
%
% 使用示例:
%   % 将个体 GM 图重采样到模板坐标空间
%   gm_in_tmpl = resample_vol_affine(gm, gm_hdr.affine, tmpl_hdr.affine, [91 109 91]);

tx = tgt_dims(1);
ty = tgt_dims(2);
tz = tgt_dims(3);

[Xt, Yt, Zt] = ndgrid(1:tx, 1:ty, 1:tz);
nTgt = tx * ty * tz;

% 目标体素坐标（0-based）→ 世界坐标
tgt_vox_0 = [Xt(:)'-1; Yt(:)'-1; Zt(:)'-1; ones(1, nTgt)];
tgt_world  = tgt_affine * tgt_vox_0;  % [4 × nTgt]

% 世界坐标 → 源体素坐标（1-based，用于 trilinear_interp）
src_vox_0 = src_affine \ tgt_world;   % [4 × nTgt]，0-based
src_vox_1 = src_vox_0(1:3,:) + 1;    % 0-based → 1-based

V_tgt = reshape(trilinear_interp(V_src, src_vox_1), tx, ty, tz);
end
