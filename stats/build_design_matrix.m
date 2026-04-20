function [X, condNames] = build_design_matrix(cfg, nScans)
% build_design_matrix - 构建一阶 GLM 设计矩阵
%
% 设计矩阵 X 的结构（列组成）:
%   [HRF卷积后的任务回归列 | 常数项（1列）| 可选漂移项]
%
% 实现说明（SPM风格）:
%   1) 在微时间分辨率 T 上构建神经活动 boxcar（units 支持 scans/secs）
%   2) 与 HRF 卷积（离散积分乘 dt）
%   3) 在每个 TR 的微时间 bin T0 处采样得到回归量
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
% 注意: 运动参数由 run_firstlevel_glm 在外部拼接

nConds = numel(cfg.cond.names);
TR     = cfg.TR;

% SPM风格微时间配置
microT = 16;
microT0 = 8;
if isfield(cfg, 'glm')
    if isfield(cfg.glm, 'microtimeResolution') && ~isempty(cfg.glm.microtimeResolution)
        microT = max(1, round(cfg.glm.microtimeResolution));
    end
    if isfield(cfg.glm, 'microtimeOnsetBin') && ~isempty(cfg.glm.microtimeOnsetBin)
        microT0 = max(1, round(cfg.glm.microtimeOnsetBin));
        microT0 = min(microT0, microT);
    end
end

dt = TR / microT;               % 微时间步长（秒）
nFine = nScans * microT;        % 总微时间点数
tFine = (0:nFine-1) * dt;       % 微时间轴（秒）

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

    if isfield(cfg, 'units') && strcmpi(cfg.units, 'scans')
        onsetsSec = onsets * TR;
        durationsSec = durations * TR;
    else
        onsetsSec = onsets;
        durationsSec = durations;
    end

    % 构建神经信号（boxcar）
    neural = zeros(nFine, 1);
    for j = 1:numel(onsetsSec)
        t_start = onsetsSec(j);
        t_end   = t_start + durationsSec(j);
        % 在高分辨率时间轴上设置 boxcar
        idx = tFine >= t_start & tFine < t_end;
        neural(idx) = neural(idx) + 1;
    end

    % 与 HRF 卷积，并进行离散积分缩放
    bold = conv(neural, hrf) * dt;
    bold = bold(1:nFine);  % 截断到原长度

    % 在每个 TR 的第 microT0 个微时间 bin 采样
    sampleIdx = (0:nScans-1) * microT + microT0;
    sampleIdx = max(1, min(nFine, sampleIdx));
    X_cond(:, c) = bold(sampleIdx);
end

% -------- 常数项 --------
X_const = ones(nScans, 1);
condNames{end+1} = 'Constant';

% -------- 漂移回归项（高通滤波的代替：DCT 基函数）--------
X_drift = [];
useExplicitDrift = true;
if isfield(cfg, 'glm') && isfield(cfg.glm, 'explicitDriftRegressors')
    useExplicitDrift = logical(cfg.glm.explicitDriftRegressors);
end
if useExplicitDrift && isfield(cfg, 'hpf') && cfg.hpf > 0
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

% -------- 可视化（可选保存为 PNG，默认关闭）--------
doSaveFig = false;
if isfield(cfg, 'debug') && isfield(cfg.debug, 'saveDesignMatrixPng')
    doSaveFig = logical(cfg.debug.saveDesignMatrixPng);
end
if doSaveFig
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
