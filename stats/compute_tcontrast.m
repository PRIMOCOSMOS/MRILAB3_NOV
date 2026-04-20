function [tMap, pMap, betaConMap, outFiles] = compute_tcontrast(beta, sigma2, X, contrast, hdr, outDir, contrastName)
% compute_tcontrast - 计算 T-contrast 统计图并写出结果
%
% T-统计量公式:
%   t = c' β̂ / √(σ̂² c' (X'X)^{-1} c)
%
%   其中:
%     c      = 对比向量（contrast vector），[nCols × 1]
%     β̂      = OLS 估计的参数，[nCols × nVox]
%     σ̂²     = 残差方差，[1 × nVox]
%     (X'X)^{-1} 是设计矩阵的协方差矩阵
%
% 输入:
%   beta         - [nCols × nVox] 参数估计（由 glm_ols 返回）
%   sigma2       - [1 × nVox] 体素噪声方差（由 glm_ols 返回）
%   X            - [nScans × nCols] 设计矩阵
%   contrast     - [nCols × 1] 对比向量（例如 [1 -1 0 0 0 0 0 0 0]）
%   hdr          - 输出 NIfTI 的头模板（空间信息）
%   outDir       - 输出目录
%   contrastName - 对比名称字符串（用于文件命名）
%
% 输出:
%   tMap      - [nx ny nz] T 统计图（NIfTI 文件已写出）
%   pMap      - [nx ny nz] 未校正 p 值图
%   betaConMap- [nx ny nz] 对比效应量图（c'β）

fprintf('[compute_tcontrast] 对比: %s\n', contrastName);

[nCols, nVox] = size(beta);
contrast = contrast(:);

% -------- 扩展对比向量（必要时补零）--------
if numel(contrast) < nCols
    contrast = [contrast; zeros(nCols - numel(contrast), 1)];
elseif numel(contrast) > nCols
    contrast = contrast(1:nCols);
end
contrast = reshape(contrast, [], 1);

% -------- 计算 β 协方差核（数值稳定）--------
% 对于 beta=pinv(X)*Y，有 Cov(beta)=sigma^2 * (pinv(X)*pinv(X)')
pX = pinv(X);                 % [nCols × nScans]
covBeta = pX * pX.';          % [nCols × nCols]

% 对比方差因子: c' Cov(beta)/sigma^2 c
cov_factor = contrast' * covBeta * contrast;  % 标量
cov_factor = max(cov_factor, 1e-12);

% -------- 计算对比效应量 c'β --------
conBeta = contrast' * beta;  % [1 × nVox]

% -------- T 统计量 --------
se   = sqrt(sigma2 * cov_factor);  % [1 × nVox] 标准误
se   = max(se, 1e-12);
tVec = conBeta ./ se;              % [1 × nVox]

% -------- 重塑为3D图 --------
nx = hdr.nx; ny = hdr.ny; nz = hdr.nz;
tMap      = reshape(tVec,      nx, ny, nz);
betaConMap= reshape(conBeta,   nx, ny, nz);

% -------- 计算 p 值（双尾，未校正）--------
df = size(X,1) - rank(X);
if df <= 0
    error('[compute_tcontrast] 自由度无效 (df=%d)，请检查设计矩阵秩与扫描数', df);
elseif df < 5
    warning('[compute_tcontrast] 自由度较低 (df=%d)，统计推断稳定性有限', df);
end
pMap = 2 * (1 - tcdf_approx(abs(tVec), df));
pMap = reshape(pMap, nx, ny, nz);
pMap(pMap < eps) = eps;

fprintf('[compute_tcontrast] T 统计完成: max(T)=%.2f, min(T)=%.2f\n', ...
    max(tVec), min(tVec));
fprintf('[compute_tcontrast] 显著体素 (|T|>3): %d / %d\n', ...
    sum(abs(tVec) > 3), nVox);

% -------- 写出 NIfTI --------
ensure_dir(outDir);

% 清理文件名（去除特殊字符）
safeName = regexprep(contrastName, '[^a-zA-Z0-9_]', '_');

% T 统计图
hdr_t = hdr;
hdr_t.nt = 1;
hdr_t.dim = int16([3, nx, ny, nz, 1, 1, 1, 1]);
hdr_t.datatype = 16;  % float32
hdr_t.scl_slope = 1; hdr_t.scl_inter = 0;
hdr_t.descrip = sprintf('T-contrast: %s', contrastName);
tFile = fullfile(outDir, sprintf('spmT_%s.nii', safeName));
nifti_write(tFile, single(tMap), hdr_t);
fprintf('[compute_tcontrast] T图: %s\n', tFile);

% p 值图（-log10(p) 形式，便于查看）
hdr_p = hdr_t;
hdr_p.descrip = sprintf('neg_log10p: %s', contrastName);
pFile = fullfile(outDir, sprintf('neg_log10p_%s.nii', safeName));
nifti_write(pFile, single(-log10(pMap)), hdr_p);
fprintf('[compute_tcontrast] -log10(p)图: %s\n', pFile);

% 对比效应量图
hdr_b = hdr_t;
hdr_b.descrip = sprintf('con_beta: %s', contrastName);
bFile = fullfile(outDir, sprintf('con_%s.nii', safeName));
nifti_write(bFile, single(betaConMap), hdr_b);
fprintf('[compute_tcontrast] 效应量图: %s\n', bFile);

% -------- 保存 SPM.mat 命名兼容的结果文件 --------
SPM_result.contrastName = contrastName;
SPM_result.contrast     = contrast;
SPM_result.df           = df;
SPM_result.tMap_file    = tFile;
SPM_result.pMap_file    = pFile;
SPM_result.con_file     = bFile;
SPM_result.max_T        = max(tVec);
SPM_result.n_sig_vox    = sum(abs(tVec) > 3);

matFile = fullfile(outDir, sprintf('SPMresult_%s.mat', safeName));
save(matFile, 'SPM_result');
fprintf('[compute_tcontrast] 结果已保存: %s\n', matFile);

outFiles = struct();
outFiles.tFile = tFile;
outFiles.pFile = pFile;
outFiles.conFile = bFile;
outFiles.matFile = matFile;
end

% ======================================================================
% 内部辅助函数
% ======================================================================

function p = tcdf_approx(t, df)
% t 分布 CDF 的近似（standalone，不需要统计工具箱）
% 使用正则化不完全 Beta 函数近似
% 对于大 df（>30），近似为正态分布

if df > 30
    % 正态近似
    p = 0.5 * (1 + erf(t / sqrt(2)));
else
    % 精确 t 分布 CDF（使用 MATLAB 内置的 betainc）
    x  = df ./ (df + t.^2);
    p  = 1 - 0.5 * betainc(x, df/2, 0.5);
    p(t < 0) = 1 - p(t < 0);
end
end
