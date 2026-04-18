function run_firstlevel_glm(smoothFile, rpFile, outDir, cfg)
% run_firstlevel_glm - 一阶 GLM 分析主函数
% 整合设计矩阵构建、OLS 估计、T-contrast 输出
%
% 输入:
%   smoothFile - 平滑后的 4D NIfTI 文件路径（标准化 + 平滑空间）
%   rpFile     - 头动参数文件（nT×6 txt，由 realign_estimate_reslice 生成）
%   outDir     - 一阶分析输出目录（如 Sub01_1stLevel）
%   cfg        - 配置结构体（见 config_sub01.m）

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
rpNames = {'rp_tx','rp_ty','rp_tz','rp_rx','rp_ry','rp_rz'};
allNames = [condNames(1:nConds), rpNames, condNames(nConds+1:end)];

% -------- 设计矩阵质量控制与秩修复（仅移除 nuisance 共线列）--------
protectIdx = 1:nConds;  % 任务条件列必须保留
if isfield(cfg, 'tcons') && isfield(cfg.tcons, 'weight')
    c = cfg.tcons.weight(:);
    nz = find(abs(c) > 0);
    protectIdx = unique([protectIdx(:); nz(:)]);
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
thresh  = prctile(meanVol(:), 30);  % 取第30百分位数作为阈值
brainMask = meanVol > thresh;

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
    warning('[run_firstlevel_glm] 脑掩模为空，回退为强度阈值掩模');
    brainMask = meanVol > thresh;
end
nVox = sum(brainMask(:));
fprintf('[run_firstlevel_glm] 脑掩模体素数: %d / %d\n', nVox, nx*ny*nz);

% -------- 展平脑内体素: [nScans × nVox] --------
Y_flat = reshape(data4d, nx*ny*nz, nScans)';  % [nScans × nVox_all]
maskIdx = find(brainMask);
Y_brain = Y_flat(:, maskIdx);  % [nScans × nVox_brain]

% -------- OLS 估计 --------
fprintf('[run_firstlevel_glm] 开始 OLS 估计...\n');
[beta_brain, res_brain, sigma2_brain] = glm_ols(Y_brain, X);

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
SPM.xBF.RT  = cfg.TR;
SPM.xX.X    = X;
SPM.xX.name = allNames;
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
    beta_all, sigma2_all, X, ...
    cfg.tcons.weight(:), ...
    hdr3d, outDir, cfg.tcons.name);

% -------- 交互式3D激活可视化（现代 Renderer）--------
if isfield(cfg, 'visualization') && isfield(cfg.visualization, 'enable') && cfg.visualization.enable
    fprintf('[run_firstlevel_glm] 生成交互式3D激活图...\n');
    render_activation_3d(contrastFiles.tFile, cfg.visualization.brainTemplateNii, outDir, cfg.visualization);
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
