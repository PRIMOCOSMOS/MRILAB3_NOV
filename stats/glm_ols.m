function [beta, res, sigma2] = glm_ols(Y, X)
% glm_ols - 一般线性模型（GLM）最小二乘（OLS）参数估计
%
% 模型: Y = X * β + ε，其中 ε ~ N(0, σ²I)
%
% OLS 解:
%   β̂ = (X'X)^{-1} X' Y  （当 X 满列秩时）
%   残差: r = Y - X*β̂
%   方差估计: σ̂² = r'r / (n - p)（每个体素独立估计）
%
% 输入:
%   Y - [nScans × nVox] 数据矩阵（nVox 可以是展平的脑体素数）
%   X - [nScans × nCols] 设计矩阵（已包含所有回归列）
%
% 输出:
%   beta   - [nCols × nVox] 参数估计（β̂）
%   res    - [nScans × nVox] 残差矩阵
%   sigma2 - [1 × nVox] 每体素的噪声方差估计
%
% 注意:
%   - 若 X 接近奇异，使用 pinv（伪逆）
%   - 大脑体素数通常很大（~100万），建议传入稀疏数据或分批处理

[nScans, nCols] = size(X);
nVox = size(Y, 2);

fprintf('[glm_ols] 开始 OLS 估计: Y=[%d × %d], X=[%d × %d]\n', ...
    nScans, nVox, nScans, nCols);

% -------- 检查矩阵秩 --------
rankX = rank(X);
if rankX < nCols
    warning('[glm_ols] 设计矩阵秩不足 (%d/%d)，使用伪逆（pinv）', rankX, nCols);
    pX = pinv(X);  % 伪逆
else
    % 高效 OLS: β̂ = (X'X)^{-1} X' Y
    % 预计算 (X'X)^{-1}X' （hat matrix 的一部分）
    XtX = X' * X;
    % 正则化（数值稳定性）
    XtX = XtX + 1e-10 * eye(nCols);
    pX  = XtX \ X';  % [(X'X)^{-1} X'] = [nCols × nScans]
end

% -------- 估计 β --------
beta = pX * Y;  % [nCols × nVox]

% -------- 计算残差 --------
res = Y - X * beta;  % [nScans × nVox]

% -------- 方差估计 --------
df      = nScans - rankX;  % 自由度
sigma2  = sum(res.^2, 1) / df;  % [1 × nVox]
sigma2  = max(sigma2, 1e-12);   % 防止零方差

fprintf('[glm_ols] 估计完成，自由度=%d\n', df);
end
