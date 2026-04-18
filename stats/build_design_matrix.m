function [X, condNames] = build_design_matrix(cfg, nScans)
% build_design_matrix - 构建一阶 GLM 设计矩阵
%
% 设计矩阵 X 的结构（列组成）:
%   [HRF卷积后的任务回归列 | 运动参数（6列）| 常数项（1列）| 可选漂移项]
%
% 原理:
%   对于每个任务条件 c:
%     1. 在时间轴上构建δ函数序列（neural signal model）:
%        n(t) = Σ_j δ(t - onset_j) * duration_j（或 boxcar 函数）
%     2. 与 HRF 卷积得到预期 BOLD 响应:
%        x(t) = n(t) * h(t)
%     3. 在 TR 时刻采样（下采样到 fMRI 采样率）
%
% 输入:
%   cfg    - 配置结构体:
%              cfg.cond.names     - 条件名称 cell array
%              cfg.cond.onsets    - 每个条件的 onset cell array（秒）
%              cfg.cond.durations - 每个条件的 duration cell array（秒）
%              cfg.TR             - 重复时间（秒）
%              cfg.hpf            - 高通滤波截止周期（秒），0=不过滤
%   nScans - 时间点数（扫描数）
%
% 输出:
%   X         - [nScans × nCols] 设计矩阵
%   condNames - 各列名称 cell array
%
% 注意: 运动参数需在 GLM 估计时额外传入（作为回归项）

nConds = numel(cfg.cond.names);
TR     = cfg.TR;
dt     = 0.1;  % 微时间分辨率（秒），用于 HRF 卷积

% 扫描时间轴（TR 采样时刻，从 0 开始）
scanTimes = (0:nScans-1) * TR;

% 高分辨率时间轴（用于 HRF 卷积）
totalTime = nScans * TR;
tFine     = 0:dt:totalTime;
nFine     = numel(tFine);

% -------- HRF --------
hrf_t  = 0:dt:32;  % HRF 时间向量（32s 覆盖完整响应）
hrf    = hrf_canonical(hrf_t);
hrf    = hrf(:);

% -------- 条件回归列 --------
X_cond = zeros(nScans, nConds);
condNames = cfg.cond.names;

for c = 1:nConds
    onsets    = cfg.cond.onsets{c}(:);
    durations = cfg.cond.durations{c}(:);

    % 构建神经信号（boxcar）
    neural = zeros(nFine, 1);
    for j = 1:numel(onsets)
        t_start = onsets(j);
        t_end   = t_start + durations(j);
        % 在高分辨率时间轴上设置 boxcar
        idx = tFine >= t_start & tFine < t_end;
        neural(idx) = 1;
    end

    % 与 HRF 卷积
    bold = conv(neural, hrf);
    bold = bold(1:nFine);  % 截断到原长度

    % 下采样到 TR 时刻
    for t = 1:nScans
        % 在扫描时刻附近插值
        [~, idx_fine] = min(abs(tFine - scanTimes(t)));
        X_cond(t, c)  = bold(idx_fine);
    end
end

% -------- 常数项 --------
X_const = ones(nScans, 1);
condNames{end+1} = 'Constant';

% -------- 漂移回归项（高通滤波的代替：DCT 基函数）--------
X_drift = [];
if isfield(cfg, 'hpf') && cfg.hpf > 0
    X_drift = build_dct_drift(nScans, TR, cfg.hpf);
    for k = 1:size(X_drift,2)
        condNames{end+1} = sprintf('Drift_%02d', k);
    end
end

% -------- 组合设计矩阵 --------
% 最终顺序: [条件 | 常数 | 漂移]
% 注意: 运动参数 (6列) 由调用者添加（在 glm_ols 中合并）
X = [X_cond, X_const, X_drift];

fprintf('[build_design_matrix] 设计矩阵构建完成: [%d × %d]\n', size(X,1), size(X,2));
fprintf('[build_design_matrix] 条件数: %d, 总列数: %d (不含运动参数)\n', nConds, size(X,2));

% -------- 可视化（保存为 PNG）--------
try
    fig = figure('Visible','off');
    imagesc(X);
    colormap(gray);
    colorbar;
    xlabel('回归列'); ylabel('扫描时间点（TR）');
    title('设计矩阵 X');
    % 添加列标签
    set(gca,'XTick',1:size(X,2));
    short_names = condNames;
    for i = 1:numel(short_names)
        if numel(short_names{i}) > 10
            short_names{i} = short_names{i}(1:10);
        end
    end
    set(gca,'XTickLabel', short_names, 'XTickLabelRotation', 45);
    saveas(fig, fullfile(fileparts(which('build_design_matrix')), ...
        '..', 'design_matrix.png'));
    close(fig);
catch
    % 保存图像不是核心功能，忽略错误
end
end

% ======================================================================
% 内部辅助函数
% ======================================================================

function X_dct = build_dct_drift(nScans, TR, hpf_sec)
% 构建 DCT（离散余弦变换）漂移基函数
% 截断频率: 1/hpf_sec Hz
% 保留周期 > hpf_sec 的成分作为漂移

% 参照 SPM(spm_filter): n = fix(2*nScans*TR/hpf_sec + 1)，并移除 DC 分量
nBases = fix(2 * nScans * TR / hpf_sec + 1);
if nBases <= 1
    X_dct = [];
    return;
end

% DCT-II 基函数（标准正交化）
% k=1 为 DC 常数项（与模型常数列共线），此处显式移除，仅保留 k=2..nBases
n_col = nBases - 1;
X_dct = zeros(nScans, n_col);
for k = 2:nBases
    norm_factor = sqrt(2/nScans);
    X_dct(:,k-1) = cos(pi * (0:nScans-1)' * (k-1) / nScans) * norm_factor;
end
end
