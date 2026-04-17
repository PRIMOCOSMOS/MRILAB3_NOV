function validate_pipeline_config(cfg)
% validate_pipeline_config - 开箱即用配置校验（fail-fast）
% 在 pipeline 开始前检查关键目录、模板、参数合法性

% -------- 必填路径检查 --------
mustExistDirs = {cfg.funRawDir, cfg.t1RawDir};
for i = 1:numel(mustExistDirs)
    if ~exist(mustExistDirs{i}, 'dir')
        error('[validate_pipeline_config] 目录不存在: %s', mustExistDirs{i});
    end
end

% -------- 模板路径检查（应用级）--------
if ~isfield(cfg, 'templates') || ~isfield(cfg.templates, 'standard')
    error('[validate_pipeline_config] 缺少 cfg.templates.standard 配置');
end
if ~isfield(cfg.templates.standard, 'brainMaskNii') || ~isfield(cfg.templates.standard, 't1TemplateNii')
    error('[validate_pipeline_config] 缺少标准模板配置: brainMaskNii / t1TemplateNii');
end
requiredStandardTemplateFiles = {
    cfg.templates.standard.brainMaskNii, ...
    cfg.templates.standard.t1TemplateNii
};
for i = 1:numel(requiredStandardTemplateFiles)
    if ~exist(requiredStandardTemplateFiles{i}, 'file')
        error('[validate_pipeline_config] 模板文件不存在, 请在配置文件中填写正确绝对路径: %s', requiredStandardTemplateFiles{i});
    end
end

if ~isfield(cfg.templates, 'dartel')
    error('[validate_pipeline_config] 缺少 cfg.templates.dartel 配置');
end
hasDualFileDartel = isfield(cfg.templates.dartel, 'gmTemplateNii') && ...
                    isfield(cfg.templates.dartel, 'wmTemplateNii') && ...
                    ~isempty(cfg.templates.dartel.gmTemplateNii) && ...
                    ~isempty(cfg.templates.dartel.wmTemplateNii);
has4DDartel = isfield(cfg.templates.dartel, 'template4DNii') && ...
              ~isempty(cfg.templates.dartel.template4DNii);
if ~hasDualFileDartel && ~has4DDartel
    error(['[validate_pipeline_config] DARTEL 模板配置缺失: ', ...
           '请提供 (gmTemplateNii + wmTemplateNii) 或 (template4DNii + gmVolumeIndex + wmVolumeIndex)']);
end
if hasDualFileDartel
    if ~exist(cfg.templates.dartel.gmTemplateNii, 'file') || ~exist(cfg.templates.dartel.wmTemplateNii, 'file')
        error('[validate_pipeline_config] DARTEL 双文件模板不存在: GM=%s, WM=%s', ...
            cfg.templates.dartel.gmTemplateNii, cfg.templates.dartel.wmTemplateNii);
    end
end
if has4DDartel
    if ~exist(cfg.templates.dartel.template4DNii, 'file')
        error('[validate_pipeline_config] DARTEL 4D模板不存在: %s', cfg.templates.dartel.template4DNii);
    end
    if ~isfield(cfg.templates.dartel, 'gmVolumeIndex') || ~isfield(cfg.templates.dartel, 'wmVolumeIndex')
        error('[validate_pipeline_config] 使用 template4DNii 时必须配置 gmVolumeIndex / wmVolumeIndex');
    end
    gmIdx = cfg.templates.dartel.gmVolumeIndex;
    wmIdx = cfg.templates.dartel.wmVolumeIndex;
    if any([gmIdx, wmIdx] < 1) || any(mod([gmIdx, wmIdx], 1) ~= 0)
        error('[validate_pipeline_config] gmVolumeIndex / wmVolumeIndex 必须为正整数');
    end
    if gmIdx == wmIdx
        error('[validate_pipeline_config] gmVolumeIndex 与 wmVolumeIndex 不能相同');
    end
end

if cfg.visualization.enable && ~exist(cfg.visualization.brainTemplateNii, 'file')
    error('[validate_pipeline_config] visualization.brainTemplateNii 不存在: %s', cfg.visualization.brainTemplateNii);
end

% -------- 参数维度与范围 --------
if numel(cfg.sliceTimingMs) ~= cfg.nSlices
    error('[validate_pipeline_config] sliceTimingMs 长度(%d) != nSlices(%d)', ...
        numel(cfg.sliceTimingMs), cfg.nSlices);
end
if cfg.refSliceIdx < 1 || cfg.refSliceIdx > cfg.nSlices
    error('[validate_pipeline_config] refSliceIdx 超出范围: %d', cfg.refSliceIdx);
end
if cfg.TR <= 0
    error('[validate_pipeline_config] TR 必须 > 0');
end
if cfg.nDummy < 0 || floor(cfg.nDummy) ~= cfg.nDummy
    error('[validate_pipeline_config] nDummy 必须为非负整数');
end
if any(cfg.fwhm <= 0)
    error('[validate_pipeline_config] fwhm 必须均 > 0');
end
if isfield(cfg, 'bet') && isfield(cfg.bet, 'percentile')
    if cfg.bet.percentile < 0 || cfg.bet.percentile > 100
        error('[validate_pipeline_config] bet.percentile 必须在 [0,100]');
    end
end
if isfield(cfg, 'bet') && isfield(cfg.bet, 'smoothSigma')
    if cfg.bet.smoothSigma < 0
        error('[validate_pipeline_config] bet.smoothSigma 必须 >= 0');
    end
end

% -------- Normalize 配置（DPABI风格）--------
if ~isfield(cfg, 'normalize')
    error('[validate_pipeline_config] 缺少 normalize 配置');
end
if ~isfield(cfg.normalize, 'boundingBox') || ~isequal(size(cfg.normalize.boundingBox), [2 3])
    error('[validate_pipeline_config] normalize.boundingBox 必须为 2x3');
end
if ~isfield(cfg.normalize, 'voxSize') || numel(cfg.normalize.voxSize) ~= 3 || any(cfg.normalize.voxSize <= 0)
    error('[validate_pipeline_config] normalize.voxSize 必须为 1x3 正数向量');
end
if any(cfg.normalize.boundingBox(2,:) <= cfg.normalize.boundingBox(1,:))
    error('[validate_pipeline_config] normalize.boundingBox 上界必须大于下界');
end

% -------- 条件设计一致性 --------
nConds = numel(cfg.cond.names);
if numel(cfg.cond.onsets) ~= nConds || numel(cfg.cond.durations) ~= nConds
    error('[validate_pipeline_config] cond.names/onsets/durations 数量不一致');
end
for c = 1:nConds
    if numel(cfg.cond.onsets{c}) ~= numel(cfg.cond.durations{c})
        error('[validate_pipeline_config] 第 %d 个条件的 onset 与 duration 数量不一致', c);
    end
end

% -------- 权重长度（至少覆盖条件列）--------
minWeightLen = nConds + 6; % run_firstlevel_glm 默认拼接6列头动参数
if numel(cfg.tcons.weight) < minWeightLen
    error('[validate_pipeline_config] tcons.weight 长度(%d) < 最小要求(%d=条件列+6头动列)', ...
        numel(cfg.tcons.weight), minWeightLen);
end
end
