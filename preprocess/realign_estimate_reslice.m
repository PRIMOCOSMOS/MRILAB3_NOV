function [outFile, rpFile] = realign_estimate_reslice(inFile, outDir, rpDir, cfg)
% realign_estimate_reslice - 头动校正（刚体 6-DOF Realignment）
%
% 物理背景:
%   fMRI 扫描中被试不可避免地有微小头动，将引入系统性信号变化。
%   头动校正将所有时间点的脑图像刚体配准到参考体（均值像或第1帧）。
%
% 数学原理:
%   刚体变换有6个自由度: [tx ty tz rx ry rz]
%   使用高斯-牛顿（Gauss-Newton）方法 + 一阶泰勒展开迭代估计参数。
%   代价函数: 最小化参考图像与变换后图像的残差平方和（均方误差）。
%
%   更新规则:
%     J = ∂I/∂p  (图像相对参数的雅可比矩阵，近似为空间梯度 × 参数梯度)
%     Δp = (J'J)^{-1} J' r   (r = 参考 - 当前)
%     p  ← p + Δp
%
% 输入:
%   inFile  - 输入 4D NIfTI 文件
%   outDir  - 重采样输出目录（前缀 'r'）
%   rpDir   - 头动参数文件输出目录
%   cfg     - 配置结构体:
%               cfg.realign.maxIter (最大迭代次数, 默认24)
%               cfg.realign.tol     (收敛阈值, 默认1e-8)
%               cfg.realign.rtm     (1=配准到均值像, 0=配准到第1帧)
%               cfg.realign.sep     (采样间距mm, 默认4)
%               cfg.subID           (被试ID，用于命名rp文件)
%
% 输出:
%   outFile - 重采样后的 4D NIfTI（前缀 'r'）
%   rpFile  - 头动参数文件路径（rp_*.txt，格式: nT×6）

fprintf('[realign] 读取: %s\n', inFile);
[data, hdr] = nifti_read(inFile);
[nx, ny, nz, nt] = size(data);

maxIter = 24;
tol     = 1e-8;
rtm     = 1;
sep_mm  = 4;
if isfield(cfg, 'realign')
    if isfield(cfg.realign, 'maxIter'), maxIter = cfg.realign.maxIter; end
    if isfield(cfg.realign, 'tol'),     tol     = cfg.realign.tol;     end
    if isfield(cfg.realign, 'rtm'),     rtm     = cfg.realign.rtm;     end
    if isfield(cfg.realign, 'sep'),     sep_mm  = cfg.realign.sep;     end
end

fprintf('[realign] 维度=[%d %d %d %d], 最大迭代=%d\n', nx, ny, nz, nt, maxIter);

% -------- 体素尺寸（mm）--------
voxSz = sqrt(sum(hdr.affine(1:3,1:3).^2, 1));  % [dx dy dz]

% -------- 采样步长（体素）--------
sep_vox = max(1, round(sep_mm ./ voxSz));  % [sx sy sz]

% -------- 采样网格（稀疏，用于快速估计）--------
xs = 1:sep_vox(1):nx;
ys = 1:sep_vox(2):ny;
zs = 1:sep_vox(3):nz;
[Xs, Ys, Zs] = ndgrid(xs, ys, zs);
nPts = numel(Xs);

% -------- 确定参考体 --------
if rtm
    % 两遍：先对第1帧配准，再以均值像为参考
    refVol = double(data(:,:,:,1));
else
    refVol = double(data(:,:,:,1));
end

% -------- 第一遍：估计头动参数 --------
params = zeros(nt, 6);  % [tx ty tz rx ry rz]

for pass = 1:(1+rtm)  % 若rtm=1则跑两遍
    if pass == 2
        % 第二遍用均值像作为参考
        resliced = reslice_all(data, params, hdr);
        refVol = mean(resliced, 4);
    end

    % 预计算参考图像及其梯度
    [gx, gy, gz] = vol_gradient(refVol);

    for t = 1:nt
        if t == 1 && pass == 1
            params(t,:) = zeros(1,6);
            continue;
        end

        curVol = double(data(:,:,:,t));

        % 高斯-牛顿 迭代
        p = params(t,:);
        for iter = 1:maxIter
            M   = rigid_mat(p);  % 当前参数的变换矩阵

            % 将参考网格坐标（1-based 体素）变换到当前图像坐标（1-based 体素）
            % 构建 [4×nPts] 齐次坐标矩阵（1-based 体素坐标）
            coordsRef_h = [Xs(:)'; Ys(:)'; Zs(:)'; ones(1,nPts)];  % [4×nPts]，1-based
            coordsCur   = M \ coordsRef_h;  % [4×nPts]，1-based（M 在体素空间近似有效）

            % 在当前图像中插值
            Ival = trilinear_interp(curVol,  coordsCur(1:3,:));  % I(cur)
            Rval = trilinear_interp(refVol,  [Xs(:)'; Ys(:)'; Zs(:)']);  % Ref

            % 残差
            r = Rval(:) - Ival(:);

            % 雅可比矩阵: ∂I/∂p = [∂I/∂x, ∂I/∂y, ∂I/∂z] × ∂(Mx)/∂p
            Gx = trilinear_interp(gx, coordsCur(1:3,:));
            Gy = trilinear_interp(gy, coordsCur(1:3,:));
            Gz = trilinear_interp(gz, coordsCur(1:3,:));

            % ∂(Mx)/∂p for rigid body (6×3 per point -> 3×6 per point)
            % J = [G] × [dM/dp]   → J is [nPts × 6]
            J = compute_jacobian(Gx(:)', Gy(:)', Gz(:)', Xs(:)', Ys(:)', Zs(:)');

            % Gauss-Newton 更新
            JtJ = J' * J + 1e-6 * eye(6);  % 正则化
            Jtr = J' * r;
            dp  = JtJ \ Jtr;

            p = p + dp';
            if norm(dp) < tol, break; end
        end

        params(t,:) = p;
        if mod(t, 50)==0
            fprintf('[realign]  已处理 %d/%d 时间点\n', t, nt);
        end
    end
end

fprintf('[realign] 参数估计完成，最大平移=%.2fmm，最大旋转=%.4frad\n', ...
    max(abs(params(:,1:3)),[],'all'), max(abs(params(:,4:6)),[],'all'));

% -------- 写出头动参数文件（rp_*.txt）--------
ensure_dir(rpDir);
rpFile = fullfile(rpDir, sprintf('rp_%s.txt', cfg.subID));
dlmwrite(rpFile, params, 'delimiter', '\t', 'precision', 8);
fprintf('[realign] 头动参数已写出: %s\n', rpFile);

% -------- 重采样所有时间点到参考空间 --------
fprintf('[realign] 开始重采样所有时间点...\n');
dataR = reslice_all(data, params, hdr);

% -------- 写出重采样数据 --------
ensure_dir(outDir);
hdr.descrip = 'Realigned';
[~, fname, ext] = fileparts(inFile);
outFile = fullfile(outDir, ['r' fname ext]);
nifti_write(outFile, single(dataR), hdr);
fprintf('[realign] 重采样完成: %s\n', outFile);
end

% ======================================================================
% 内部辅助函数
% ======================================================================

function [gx, gy, gz] = vol_gradient(vol)
% 使用有限差分计算体数据的空间梯度
[gx, gy, gz] = gradient(vol);
end

function J = compute_jacobian(Gx, Gy, Gz, X, Y, Z)
% 计算刚体变换参数的雅可比矩阵
% 对 6 个参数 [tx ty tz rx ry rz] 求偏导
% 每行对应一个采样点，共 6 列
%
% 简化：使用小角近似（线性化）
%   ∂x/∂tx=1, ∂y/∂ty=1, ∂z/∂tz=1
%   ∂x/∂rx=0,      ∂y/∂rx=-z, ∂z/∂rx=y  (绕X旋转)
%   ∂x/∂ry=z,      ∂y/∂ry=0,  ∂z/∂ry=-x (绕Y旋转)
%   ∂x/∂rz=-y,     ∂y/∂rz=x,  ∂z/∂rz=0  (绕Z旋转)

nPts = numel(Gx);
J = zeros(nPts, 6);

% 平移分量（直接对应梯度）
J(:,1) = Gx(:);    % ∂I/∂tx
J(:,2) = Gy(:);    % ∂I/∂ty
J(:,3) = Gz(:);    % ∂I/∂tz

% 旋转分量（链式法则）
J(:,4) = Gy(:) .* (-Z(:)) + Gz(:) .* Y(:);   % ∂I/∂rx
J(:,5) = Gx(:) .* Z(:)    + Gz(:) .* (-X(:)); % ∂I/∂ry
J(:,6) = Gx(:) .* (-Y(:)) + Gy(:) .* X(:);   % ∂I/∂rz
end

function dataR = reslice_all(data, params, hdr)
% 将所有时间点重采样到参考空间
[nx, ny, nz, nt] = size(data);
dataR = zeros(nx, ny, nz, nt, 'single');

% 参考空间网格（1-based 体素坐标）
[Xr, Yr, Zr] = ndgrid(1:nx, 1:ny, 1:nz);
refCoords = [Xr(:)'; Yr(:)'; Zr(:)'; ones(1,nx*ny*nz)];

for t = 1:nt
    M = rigid_mat(params(t,:));

    % 参考坐标 → 当前图像坐标
    curCoords = M \ refCoords;

    vals = trilinear_interp(double(data(:,:,:,t)), curCoords(1:3,:));
    dataR(:,:,:,t) = single(reshape(vals, nx, ny, nz));
end
end
