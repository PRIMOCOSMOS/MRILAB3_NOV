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

if use_spm_functional_backend(cfg)
    [outFile, rpFile] = realign_estimate_reslice_spm(inFile, outDir, rpDir, cfg);
    return;
end

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
refCoords = [Xs(:)'; Ys(:)'; Zs(:)'];  % [3 x nPts], 1-based
nPts = size(refCoords, 2);

% 在图像中心坐标系下估计旋转可显著降低平移-旋转耦合
center = (double([nx; ny; nz]) + 1) / 2;  % 1-based 体素中心
refCentered = refCoords - center;

% -------- 确定参考体 --------
% 两遍策略时，第一遍也以第1帧为初始参考
refVol = double(data(:,:,:,1));

% -------- 估计头动参数（内部参数：体素平移 + 弧度旋转）--------
paramsVox = zeros(nt, 6);  % [tx_vox ty_vox tz_vox rx ry rz]

for pass = 1:(1+rtm)  % 若rtm=1则跑两遍
    if pass == 2
        % 第二遍用均值像作为参考
        resliced = reslice_all(data, paramsVox, center);
        refVol = mean(resliced, 4);
    end

    % 参考体在采样点上的强度（当前 pass 固定）
    refVals = trilinear_interp(refVol, refCoords);

    for t = 1:nt
        if t == 1 && pass == 1
            paramsVox(t,:) = zeros(1,6);
            continue;
        end

        curVol = double(data(:,:,:,t));
        [cgx, cgy, cgz] = vol_gradient(curVol);

        % 以相邻时间点为初值，符合 fMRI 头动的时间连续性
        if pass == 1 && t > 1
            p = paramsVox(t-1,:);
        else
            p = paramsVox(t,:);
        end

        bestP = p;
        bestCost = inf;
        stagnationCount = 0;

        % 阻尼高斯-牛顿迭代
        for iter = 1:maxIter
            coordsCur = transform_ref_to_cur(refCoords, p, center);
            validMask = coordsCur(1,:) >= 1 & coordsCur(1,:) <= nx & ...
                        coordsCur(2,:) >= 1 & coordsCur(2,:) <= ny & ...
                        coordsCur(3,:) >= 1 & coordsCur(3,:) <= nz;

            if nnz(validMask) < 256
                warning('[realign] t=%d 有效重叠采样点不足（%d/%d），回退最佳参数', t, nnz(validMask), nPts);
                p = bestP;
                break;
            end

            coordsCurValid = coordsCur(:, validMask);
            refValid = refVals(validMask);

            % 残差定义：r = I_ref - I_cur(T(x,p))
            Ival = trilinear_interp(curVol, coordsCurValid);
            r = refValid(:) - Ival(:);
            curCost = mean(r.^2);

            if curCost < bestCost
                bestCost = curCost;
                bestP = p;
                stagnationCount = 0;
            else
                stagnationCount = stagnationCount + 1;
            end

            % 使用“当前图像”梯度（符合 Gauss-Newton 线性化）
            Gx = trilinear_interp(cgx, coordsCurValid);
            Gy = trilinear_interp(cgy, coordsCurValid);
            Gz = trilinear_interp(cgz, coordsCurValid);

            J = compute_jacobian( ...
                Gx(:)', Gy(:)', Gz(:)', ...
                refCentered(1, validMask), ...
                refCentered(2, validMask), ...
                refCentered(3, validMask));

            JtJ = J' * J;
            damping = 1e-4 * (trace(JtJ) / 6 + eps);
            dp = (JtJ + damping * eye(6)) \ (J' * r);

            % 步长限制：避免一次更新过大导致发散
            dp = limit_update_step(dp);
            p = clamp_parameter_bounds(p + dp');

            if norm(dp) < tol || stagnationCount >= 4
                p = bestP;
                break;
            end
        end

        paramsVox(t,:) = bestP;
        if mod(t, 50)==0
            fprintf('[realign]  已处理 %d/%d 时间点\n', t, nt);
        end
    end
end

% rp 文件使用 mm/rad 约定（与 SPM 输出习惯一致）
params = paramsVox;
params(:,1) = paramsVox(:,1) * voxSz(1);
params(:,2) = paramsVox(:,2) * voxSz(2);
params(:,3) = paramsVox(:,3) * voxSz(3);

fprintf('[realign] 参数估计完成，最大平移=%.2fmm，最大旋转=%.4frad\n', ...
    max(abs(params(:,1:3)),[],'all'), max(abs(params(:,4:6)),[],'all'));

% -------- 写出头动参数文件（rp_*.txt）--------
ensure_dir(rpDir);
rpFile = fullfile(rpDir, sprintf('rp_%s.txt', cfg.subID));
if exist('writematrix', 'file')
    writematrix(params, rpFile, 'Delimiter', 'tab');
else
    dlmwrite(rpFile, params, 'delimiter', '\t', 'precision', 8); %#ok<DLMWRT>
end
fprintf('[realign] 头动参数已写出: %s\n', rpFile);

% -------- 重采样所有时间点到参考空间 --------
fprintf('[realign] 开始重采样所有时间点...\n');
dataR = reslice_all(data, paramsVox, center);

% -------- 写出重采样数据 --------
ensure_dir(outDir);
hdr.descrip = 'Realigned';
[~, fname, ext] = fileparts(inFile);
outFile = fullfile(outDir, ['r' fname ext]);
nifti_write(outFile, single(dataR), hdr);
fprintf('[realign] 重采样完成: %s\n', outFile);

% 与 SPM 行为对齐：额外写出 mean 图像，便于 QC 与后续流程对照
hdrMean = hdr;
hdrMean.nt = 1;
hdrMean.dim = int16([3, nx, ny, nz, 1, 1, 1, 1]);
hdrMean.descrip = 'MeanRealigned';
meanFile = fullfile(outDir, ['mean' fname ext]);
nifti_write(meanFile, single(mean(dataR, 4)), hdrMean);
fprintf('[realign] 均值像已写出: %s\n', meanFile);
end

% ======================================================================
% 内部辅助函数
% ======================================================================

function [gx, gy, gz] = vol_gradient(vol)
% 使用有限差分计算体数据的空间梯度
[gx, gy, gz] = gradient(vol);
end

function J = compute_jacobian(Gx, Gy, Gz, Xc, Yc, Zc)
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
J(:,4) = Gy(:) .* (-Zc(:)) + Gz(:) .* Yc(:);    % ∂I/∂rx
J(:,5) = Gx(:) .* Zc(:)    + Gz(:) .* (-Xc(:)); % ∂I/∂ry
J(:,6) = Gx(:) .* (-Yc(:)) + Gy(:) .* Xc(:);    % ∂I/∂rz
end

function dataR = reslice_all(data, paramsVox, center)
% 将所有时间点重采样到参考空间
[nx, ny, nz, nt] = size(data);
dataR = zeros(nx, ny, nz, nt, 'single');

% 参考空间网格（1-based 体素坐标）
[Xr, Yr, Zr] = ndgrid(1:nx, 1:ny, 1:nz);
refCoords = [Xr(:)'; Yr(:)'; Zr(:)'];

for t = 1:nt
    curCoords = transform_ref_to_cur(refCoords, paramsVox(t,:), center);

    vals = trilinear_interp(double(data(:,:,:,t)), curCoords);
    dataR(:,:,:,t) = single(reshape(vals, nx, ny, nz));
end

% 基本质量控制：防止重采样后出现“纯白/近常数/非有限值”异常
finiteRatio = nnz(isfinite(dataR)) / numel(dataR);
if finiteRatio < 1
    error('[realign] 重采样结果包含非有限值 (finiteRatio=%.6f)，请检查配准变换', finiteRatio);
end
volStd = squeeze(std(reshape(double(dataR), [], nt), 0, 1));
% 经验阈值：float32 数据在正常 fMRI 强度范围下，体数据标准差远大于 1e-6；
% 该阈值仅用于捕获“近乎常数图像（纯白/纯黑）”这类明显失败情形。
MIN_ACCEPTABLE_STD = 1e-6;
badIdx = find(volStd < MIN_ACCEPTABLE_STD, 1, 'first');
if ~isempty(badIdx)
    error('[realign] 重采样结果包含近似常数图像 (t=%d, std=%.3e)，请检查配准参数与坐标变换', ...
      badIdx, volStd(badIdx));
end
end

function coordsCur = transform_ref_to_cur(refCoords, p, center)
% 将参考坐标映射到当前图像坐标（均为1-based体素坐标）
R = rotation_matrix(p(4), p(5), p(6));
t = p(1:3).';
coordsCur = R * (refCoords - center) + center + t;
end

function R = rotation_matrix(rx, ry, rz)
% 与 rigid_mat 保持一致的旋转顺序：R = Rz * Ry * Rx
Rx = [1,      0,       0;
    0,  cos(rx), -sin(rx);
    0,  sin(rx),  cos(rx)];
Ry = [ cos(ry), 0, sin(ry);
        0,  1,      0;
    -sin(ry), 0, cos(ry)];
Rz = [cos(rz), -sin(rz), 0;
    sin(rz),  cos(rz), 0;
         0,        0,  1];
R = Rz * Ry * Rx;
end

function dp = limit_update_step(dp)
% 迭代步长限制，防止大步长导致优化发散
dp = dp(:);
MAX_TRANS_STEP_VOX = 0.5;
MAX_ROT_STEP_RAD = 0.01;
dp(1:3) = max(min(dp(1:3), MAX_TRANS_STEP_VOX), -MAX_TRANS_STEP_VOX);
dp(4:6) = max(min(dp(4:6), MAX_ROT_STEP_RAD), -MAX_ROT_STEP_RAD);
end

function p = clamp_parameter_bounds(p)
% 参数绝对范围限制（内部单位：voxel + rad）
MAX_TRANS_ABS_VOX = 20;
MAX_ROT_ABS_RAD = 0.35;
p(1:3) = max(min(p(1:3), MAX_TRANS_ABS_VOX), -MAX_TRANS_ABS_VOX);
p(4:6) = max(min(p(4:6), MAX_ROT_ABS_RAD), -MAX_ROT_ABS_RAD);
end

function tf = use_spm_functional_backend(cfg)
tf = isstruct(cfg) && isfield(cfg, 'spm') && ...
    isfield(cfg.spm, 'useFunctional') && logical(cfg.spm.useFunctional);
end
