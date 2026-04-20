function run_firstlevel_glm(smoothFile, rpFile, outDir, cfg)
% run_firstlevel_glm - 一阶 GLM 分析主函数
% 整合设计矩阵构建、OLS 估计、T-contrast 输出
%
% 输入:
%   smoothFile - 平滑后的 4D NIfTI 文件路径（标准化 + 平滑空间）
%   rpFile     - 头动参数文件（nT×6 txt，由 realign_estimate_reslice 生成）
%   outDir     - 一阶分析输出目录（如 Sub01_1stLevel）
%   cfg        - 配置结构体（见 config_sub01.m）

if use_spm_functional_backend(cfg)
    run_firstlevel_glm_spm(smoothFile, rpFile, outDir, cfg);
    return;
end

fprintf('[run_firstlevel_glm] 开始一阶 GLM 分析\n');
fprintf('[run_firstlevel_glm] 输入: %s\n', smoothFile);

% -------- 读取数据 --------
[data4d, hdr] = nifti_read(smoothFile);
[nx, ny, nz, nScans] = size(data4d);
fprintf('[run_firstlevel_glm] 4D 维度: [%d %d %d %d]\n', nx, ny, nz, nScans);

% -------- 读取头动参数 --------
rp = load(rpFile);  % [nScans × 6]
if size(rp,1) ~= nScans
    error('[run_firstlevel_glm] 头动参数行数 %d 与扫描数 %d 不匹配', ...
        size(rp,1), nScans);
end
fprintf('[run_firstlevel_glm] 已读取头动参数: %s\n', rpFile);

% -------- 构建设计矩阵 --------
fprintf('[run_firstlevel_glm] 构建设计矩阵...\n');
[X_base, condNames] = build_design_matrix(cfg, nScans);

% 合并运动参数（6列）：插入在条件列之后、常数项之前
nConds = numel(cfg.cond.names);
% X_base = [条件(nConds) | 常数(1) | 漂移(k)]
% 插入运动参数：[条件 | 运动(6) | 常数 | 漂移]
% 说明: 该顺序与 SPM 常用建模习惯一致，便于保持条件列索引稳定，
% 同时将运动参数作为 nuisance regressors 放在任务条件后、漂移项前。
X = [X_base(:, 1:nConds), rp, X_base(:, nConds+1:end)];

% 更新列名
spmCondNames = cell(1, nConds);
for i = 1:nConds
    spmCondNames{i} = sprintf('Sn(1) %s*bf(1)', condNames{i});
end
rpNames = arrayfun(@(k) sprintf('Sn(1) R%d', k), 1:6, 'UniformOutput', false);
baseTailNames = condNames(nConds+1:end);
if ~isempty(baseTailNames)
    baseTailNames{1} = 'Sn(1) constant';
    for k = 2:numel(baseTailNames)
        baseTailNames{k} = sprintf('Sn(1) %s', baseTailNames{k});
    end
end
allNames = [spmCondNames, rpNames, baseTailNames];

% -------- 设计矩阵质量控制与秩修复（仅移除 nuisance 共线列）--------
protectIdx = 1:nConds;  % 任务条件列必须保留
if isfield(cfg, 'tcons') && isfield(cfg.tcons, 'weight')
    c = cfg.tcons.weight(:);
    cIdxNonZero = find(abs(c) > 0);
    protectIdx = unique([protectIdx(:); cIdxNonZero(:)]);
    protectIdx = protectIdx(protectIdx >= 1 & protectIdx <= size(X,2))';
end
[X, allNames] = enforce_fullrank_design(X, allNames, protectIdx);

fprintf('[run_firstlevel_glm] 设计矩阵: [%d × %d]\n', size(X,1), size(X,2));
rankX = rank(X);
rcondXtX = rcond(X' * X + 1e-12 * eye(size(X,2)));
fprintf('[run_firstlevel_glm] 设计矩阵秩: %d / %d, rcond(X''X)=%.3e\n', ...
    rankX, size(X,2), rcondXtX);
if rankX < size(X,2)
    error('[run_firstlevel_glm] 设计矩阵仍秩不足 (%d/%d)，请检查任务设计/协变量', rankX, size(X,2));
end

% -------- 保存设计矩阵图像 --------
ensure_dir(outDir);
try
    fig = figure('Visible','off');
    imagesc(X);
    colormap(gray); colorbar;
    title(sprintf('设计矩阵 [%d × %d]', size(X,1), size(X,2)));
    xlabel('回归列'); ylabel('扫描时间点');
    saveas(fig, fullfile(outDir, 'design_matrix.png'));
    close(fig);
catch; end

% -------- 脑掩模（去除背景体素）--------
meanVol = mean(data4d, 4);
maskMethod = 'percentile';
if isfield(cfg, 'glm') && isfield(cfg.glm, 'maskMethod')
    maskMethod = cfg.glm.maskMethod;
end
if isstring(maskMethod)
    maskMethod = char(maskMethod);
end
if iscell(maskMethod)
    maskMethod = maskMethod{1};
end
maskMethod = lower(char(maskMethod));

switch maskMethod
    case 'percentile'
        p = 30;
        if isfield(cfg, 'glm') && isfield(cfg.glm, 'maskPercentile')
            p = cfg.glm.maskPercentile;
        end
        p = min(max(p, 0), 100);
        thresh = prctile(meanVol(:), p);
        baseMask = meanVol > thresh;
        fprintf('[run_firstlevel_glm] 掩模策略: percentile(%.1f), 阈值=%.4g\n', p, thresh);
    case 'globalfraction'
        frac = 0.5;
        if isfield(cfg, 'glm') && isfield(cfg.glm, 'maskGlobalFraction')
            frac = cfg.glm.maskGlobalFraction;
        end
        frac = min(max(frac, 0), 1);
        posVals = meanVol(meanVol > 0);
        if isempty(posVals)
            gMean = mean(meanVol(:));
        else
            gMean = mean(posVals(:));
        end
        thresh = frac * gMean;
        baseMask = meanVol > thresh;
        fprintf('[run_firstlevel_glm] 掩模策略: globalFraction(%.3f), 全局均值=%.4g, 阈值=%.4g\n', frac, gMean, thresh);
    otherwise
        error('[run_firstlevel_glm] 未知掩模策略: %s', maskMethod);
end

brainMask = baseMask;

% 优先叠加标准空间脑掩模（若已配置）
if isfield(cfg, 'templates') && isfield(cfg.templates, 'standard') && ...
   isfield(cfg.templates.standard, 'brainMaskNii') && ...
   exist(cfg.templates.standard.brainMaskNii, 'file')
    try
        stdMaskThreshold = 0.5;
        [maskData, maskHdr] = nifti_read(cfg.templates.standard.brainMaskNii);
        maskVol = double(maskData(:,:,:,1));
        if any(size(maskVol) ~= [nx ny nz])
            % 使用仿射矩阵对齐重采样，确保脑掩模与功能像在同一坐标系
            maskVol = resample_vol_affine(maskVol, maskHdr.affine, hdr.affine, [nx ny nz]);
            fprintf('[run_firstlevel_glm] 标准脑掩模已仿射重采样到 [%d %d %d]\n', nx, ny, nz);
        end
        stdMask = maskVol >= stdMaskThreshold;
        if any(stdMask(:))
            brainMask = brainMask & stdMask;
            fprintf('[run_firstlevel_glm] 已应用标准脑掩模: %s\n', cfg.templates.standard.brainMaskNii);
        else
            warning('[run_firstlevel_glm] 标准脑掩模为空，回退到强度阈值掩模');
        end
    catch ME
        warning('[run_firstlevel_glm] 读取/应用标准脑掩模失败，回退到强度阈值掩模: %s', ME.message);
    end
end

if ~any(brainMask(:))
    warning('[run_firstlevel_glm] 脑掩模为空，回退为基础掩模');
    brainMask = baseMask;
end
if ~any(brainMask(:))
    warning('[run_firstlevel_glm] 基础掩模为空，回退为 meanVol>0');
    brainMask = meanVol > 0;
end
nVox = sum(brainMask(:));
if ~isscalar(nx) || ~isscalar(ny) || ~isscalar(nz) || ~isscalar(nScans)
    error('[run_firstlevel_glm] 数据维度异常（非标量），请检查上游变量覆盖');
end
fprintf('[run_firstlevel_glm] 脑掩模体素数: %d / %d\n', nVox, nx*ny*nz);

% -------- 展平脑内体素: [nScans × nVox] --------
Y_flat = reshape(data4d, nx*ny*nz, nScans)';  % [nScans × nVox_all]
maskIdx = find(brainMask);
Y_brain = Y_flat(:, maskIdx);  % [nScans × nVox_brain]

% -------- 可选：SPM风格高通滤波（作用于 X 与 Y）--------
X_est = X;
Y_est = Y_brain;
K = [];
applyHPF = false;
if isfield(cfg, 'glm') && isfield(cfg.glm, 'applyHighPassFilter')
    applyHPF = logical(cfg.glm.applyHighPassFilter);
end
if applyHPF && isfield(cfg, 'hpf') && cfg.hpf > 0
    [X_est, Y_est, K] = apply_highpass_filter(X, Y_brain, cfg.TR, cfg.hpf);
    fprintf('[run_firstlevel_glm] 已应用高通滤波到估计矩阵: HParam=%.1fs\n', cfg.hpf);
end

% -------- 可选：AR(1) 预白化（SPM风格串行相关修正）--------
ar1Info = struct('enabled', false, 'rho', 0, 'nVoxUsed', 0);
useAR1 = true;
if isfield(cfg, 'glm') && isfield(cfg.glm, 'useAR1Whitening')
    useAR1 = logical(cfg.glm.useAR1Whitening);
end
if useAR1
    [X_est, Y_est, ar1Info] = apply_ar1_prewhiten(X_est, Y_est);
    if ar1Info.enabled
        fprintf('[run_firstlevel_glm] 已应用 AR(1) 预白化: rho=%.4f (vox=%d)\n', ar1Info.rho, ar1Info.nVoxUsed);
    else
        fprintf('[run_firstlevel_glm] AR(1) 预白化未启用（rho≈0 或估计失败）\n');
    end
end

% -------- OLS 估计 --------
fprintf('[run_firstlevel_glm] 开始 OLS 估计...\n');
[beta_brain, res_brain, sigma2_brain] = glm_ols(Y_est, X_est);

% -------- 将结果填回3D空间 --------
nColsX = size(X,2);
beta_all   = zeros(nColsX, nx*ny*nz);
sigma2_all = zeros(1, nx*ny*nz);
beta_all(:, maskIdx)   = beta_brain;
sigma2_all(maskIdx)    = sigma2_brain;

% -------- 重塑 beta 为4D（第4维=参数列）--------
beta4d = reshape(beta_all', nx, ny, nz, nColsX);

% 写出 beta 图（每列一个3D NIfTI）
for k = 1:nColsX
    hdr_b = hdr;
    hdr_b.nt = 1;
    hdr_b.dim = int16([3, nx, ny, nz, 1, 1, 1, 1]);
    hdr_b.descrip = sprintf('beta_%02d_%s', k, allNames{min(k,numel(allNames))});
    betaFile = fullfile(outDir, sprintf('beta_%04d.nii', k));
    nifti_write(betaFile, single(beta4d(:,:,:,k)), hdr_b);
end
fprintf('[run_firstlevel_glm] 已写出 %d 个 beta 图\n', nColsX);

% 写出方差图
hdr_s2 = hdr;
hdr_s2.nt = 1;
hdr_s2.dim = int16([3, nx, ny, nz, 1, 1, 1, 1]);
hdr_s2.descrip = 'ResMS (residual variance)';
nifti_write(fullfile(outDir, 'ResMS.nii'), ...
    single(reshape(sigma2_all, nx, ny, nz)), hdr_s2);

% -------- 保存 SPM.mat 兼容结构 --------
SPM.swd     = outDir;
SPM.xY.P    = smoothFile;
SPM.xY.RT   = cfg.TR;
SPM.xBF.RT  = cfg.TR;
SPM.xBF.UNITS = cfg.units;
if isfield(cfg, 'glm') && isfield(cfg.glm, 'microtimeResolution')
    SPM.xBF.T = cfg.glm.microtimeResolution;
else
    SPM.xBF.T = 16;
end
if isfield(cfg, 'glm') && isfield(cfg.glm, 'microtimeOnsetBin')
    SPM.xBF.T0 = cfg.glm.microtimeOnsetBin;
else
    SPM.xBF.T0 = 8;
end
SPM.xBF.name = 'hrf';
SPM.xX.X    = X;
SPM.xX.name = allNames;
if ~isempty(K)
    SPM.xX.K = K;
end
if ar1Info.enabled
    SPM.xVi.form = 'AR(1)';
    SPM.xVi.rho = ar1Info.rho;
end
SPM.nscan = nScans;
SPM.Sess(1).U = make_session_u(cfg);
SPM.beta    = beta4d;
SPM.VResMS.fname = fullfile(outDir, 'ResMS.nii');
SPM.xCon    = [];  % 对比结构（后续添加）
matFile = fullfile(outDir, 'SPM.mat');
save(matFile, 'SPM');
fprintf('[run_firstlevel_glm] SPM.mat 已保存: %s\n', matFile);

% -------- 计算 T-contrast --------
% 头信息（用于写出T图）
hdr3d = hdr;
hdr3d.nt = 1;
hdr3d.nx = nx; hdr3d.ny = ny; hdr3d.nz = nz;
hdr3d.dim = int16([3, nx, ny, nz, 1, 1, 1, 1]);

fprintf('[run_firstlevel_glm] 计算 T-contrast: %s\n', cfg.tcons.name);
[tMap, pMap, ~, contrastFiles] = compute_tcontrast(...
    beta_all, sigma2_all, X_est, ...
    cfg.tcons.weight(:), ...
    hdr3d, outDir, cfg.tcons.name);

% -------- 交互式3D激活可视化（现代 Renderer）--------
if isfield(cfg, 'visualization') && isfield(cfg.visualization, 'enable') && cfg.visualization.enable
    fprintf('[run_firstlevel_glm] 生成交互式3D激活图...\n');
    render_activation_3d(contrastFiles.tFile, cfg.visualization.brainTemplateNii, outDir, cfg.visualization);
end

function [Xf, Yf, K] = apply_highpass_filter(X, Y, TR, hpfSec)
% 单会话高通滤波（SPM风格）：R = I - X0*pinv(X0), X0 为 DCT 低频基（去掉常数列）
nScans = size(X, 1);
nBases = fix(2 * nScans * TR / hpfSec + 1);
if nBases <= 1
    Xf = X;
    Yf = Y;
    K = struct('HParam', hpfSec, 'X0', [], 'R', eye(nScans));
    return;
end

t = (0:nScans-1)';
X0 = zeros(nScans, nBases);
X0(:,1) = 1 / sqrt(nScans);
for k = 2:nBases
    X0(:,k) = sqrt(2 / nScans) * cos(pi * (2*t + 1) * (k-1) / (2*nScans));
end

% 去掉 DC 常数项，保留低频漂移基
X0 = X0(:,2:end);
if isempty(X0)
    R = eye(nScans);
else
    R = eye(nScans) - X0 * pinv(X0);
end

Xf = R * X;
Yf = R * Y;

K = struct();
K.HParam = hpfSec;
K.X0 = X0;
K.R = R;
end

function U = make_session_u(cfg)
% 生成与 SPM.Sess(1).U 结构近似的条件描述
nConds = numel(cfg.cond.names);
U = repmat(struct('name',{{''}}, 'ons',[], 'dur',[]), nConds, 1);
for i = 1:nConds
    U(i).name = {cfg.cond.names{i}};
    U(i).ons = cfg.cond.onsets{i}(:);
    U(i).dur = cfg.cond.durations{i}(:);
end
end

function [Xw, Yw, info] = apply_ar1_prewhiten(X, Y)
% 估计全局 AR(1) 系数并对白化 X/Y。
% 估计策略：先用 OLS 残差求每体素 lag-1 自相关，再取稳健中位数。
nScans = size(X, 1);
if nScans < 3
    Xw = X;
    Yw = Y;
    info = struct('enabled', false, 'rho', 0, 'nVoxUsed', 0);
    return;
end

% 初始残差（不打印额外日志，避免污染主流程输出）
pX = pinv(X);
res0 = Y - X * (pX * Y);

num = sum(res0(2:end,:) .* res0(1:end-1,:), 1);
den = sum(res0(1:end-1,:).^2, 1) + eps;
rhoVec = num ./ den;
rhoVec = rhoVec(isfinite(rhoVec));

if isempty(rhoVec)
    Xw = X;
    Yw = Y;
    info = struct('enabled', false, 'rho', 0, 'nVoxUsed', 0);
    return;
end

rho = median(rhoVec);
rho = min(max(rho, -0.5), 0.95);

if abs(rho) < 1e-4
    Xw = X;
    Yw = Y;
    info = struct('enabled', false, 'rho', rho, 'nVoxUsed', numel(rhoVec));
    return;
end

Xw = zeros(size(X), class(X));
Yw = zeros(size(Y), class(Y));

firstScale = sqrt(max(1 - rho^2, eps));
Xw(1,:) = firstScale * X(1,:);
Yw(1,:) = firstScale * Y(1,:);
Xw(2:end,:) = X(2:end,:) - rho * X(1:end-1,:);
Yw(2:end,:) = Y(2:end,:) - rho * Y(1:end-1,:);

info = struct('enabled', true, 'rho', rho, 'nVoxUsed', numel(rhoVec));
end

fprintf('[run_firstlevel_glm] === 一阶 GLM 分析完成 ===\n');
fprintf('[run_firstlevel_glm] 输出目录: %s\n', outDir);
end

function [X_out, names_out] = enforce_fullrank_design(X_in, names_in, protectIdx)
% 仅通过移除 nuisance 列修复秩缺陷；保护任务/对比涉及列
X_out = X_in;
names_out = names_in;
nCols = size(X_out, 2);
protectIdx = unique(protectIdx(:))';
protectIdx = protectIdx(protectIdx >= 1 & protectIdx <= nCols);

if rank(X_out) == nCols
    return;
end

keepMask = true(1, nCols);
nuisance = setdiff(1:nCols, protectIdx, 'stable');
if isempty(nuisance)
    error('[run_firstlevel_glm] 所有列都被保护，无法移除共线列修复秩缺陷');
end

while true
    curCols = find(keepMask);
    curRank = rank(X_out(:, curCols));
    if curRank == numel(curCols)
        break;
    end

    bestCol = 0;
    bestRank = -1;
    canDrop = nuisance(keepMask(nuisance));
    for i = 1:numel(canDrop)
        c = canDrop(i);
        km = keepMask;
        km(c) = false;
        r = rank(X_out(:, km));
        if r > bestRank
            bestRank = r;
            bestCol = c;
        end
    end

    if bestCol == 0 || bestRank <= curRank
        break;
    end
    keepMask(bestCol) = false;
end

keptCols = find(keepMask);
droppedCols = find(~keepMask);

if rank(X_out(:, keptCols)) < numel(keptCols)
    error('[run_firstlevel_glm] 无法在保留任务/对比列前提下获得满秩设计矩阵');
end

X_out = X_out(:, keptCols);
names_out = names_out(keptCols);

if ~isempty(droppedCols)
    droppedTxt = strjoin(arrayfun(@(i) sprintf('%d:%s', i, names_in{i}), droppedCols, 'UniformOutput', false), ', ');
    warning('[run_firstlevel_glm] 检测到共线 nuisance 列并已移除: %s', droppedTxt);
end
end

function tf = use_spm_functional_backend(cfg)
tf = isstruct(cfg) && isfield(cfg, 'spm') && ...
     isfield(cfg.spm, 'useFunctional') && logical(cfg.spm.useFunctional);
end
