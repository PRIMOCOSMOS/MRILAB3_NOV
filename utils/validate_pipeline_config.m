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
requiredTemplateFiles = {
    cfg.templates.dartel.gmTemplateNii, ...
    cfg.templates.dartel.wmTemplateNii, ...
    cfg.templates.standard.brainMaskNii, ...
    cfg.templates.standard.t1TemplateNii
};
for i = 1:numel(requiredTemplateFiles)
    if ~exist(requiredTemplateFiles{i}, 'file')
        error('[validate_pipeline_config] 模板文件不存在, 请在配置文件中填写正确绝对路径: %s', requiredTemplateFiles{i});
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

% -------- 权重长度（至少覆盖条件列）--------
minWeightLen = numel(cfg.cond.names) + 6; % run_firstlevel_glm 默认拼接6列头动参数
if numel(cfg.tcons.weight) < minWeightLen
    error('[validate_pipeline_config] tcons.weight 长度(%d) < 最小要求(%d=条件列+6头动列)', ...
        numel(cfg.tcons.weight), minWeightLen);
end
end
