function cfg = resolve_pipeline_template_paths(cfg)
% resolve_pipeline_template_paths - 根据安装根目录自动推断并修正模板路径
% 目标：
% 1) 兼容 DPABI Templates 与 DPABI/Templates/SPMTemplates 两套目录结构
% 2) 兼容用户本机 SPM 安装目录（如 D:\spm）
% 3) 避免硬编码单一路径导致模板缺失

% ---- 安装根目录（可在 config 中覆盖）----
dpabiRoot = 'D:\DPABI_V9.0_250415';
spmRoot   = 'D:\spm';
if isfield(cfg, 'installPaths')
    if isfield(cfg.installPaths, 'dpabiRoot') && ~isempty(cfg.installPaths.dpabiRoot)
        dpabiRoot = cfg.installPaths.dpabiRoot;
    end
    if isfield(cfg.installPaths, 'spmRoot') && ~isempty(cfg.installPaths.spmRoot)
        spmRoot = cfg.installPaths.spmRoot;
    end
end

% ---- DARTEL 模板（4D）----
dartelCandidates = unique(nonempty({
    getfield_safe(cfg, {'templates','dartel','template4DNii'}), ...
    fullfile(spmRoot,   'toolbox', 'DARTEL', 'Template_6_IXI555_MNI152.nii'), ...
    fullfile(dpabiRoot, 'Templates', 'SPMTemplates', 'toolbox', 'DARTEL', 'Template_6_IXI555_MNI152.nii')
}));
[dartelPath, dartelFound] = first_existing_file(dartelCandidates);
if dartelFound
    cfg.templates.dartel.template4DNii = dartelPath;
end

% ---- 标准脑掩模（优先 DPABI 的 BrainMask_05_91x109x91）----
brainMaskCandidates = unique(nonempty({
    getfield_safe(cfg, {'templates','standard','brainMaskNii'}), ...
    fullfile(dpabiRoot, 'Templates', 'BrainMask_05_91x109x91.hdr'), ...
    fullfile(dpabiRoot, 'Templates', 'SPMTemplates', 'toolbox', 'OldNorm', 'mask_ICV.hdr')
}));
[brainMaskPath, brainMaskFound] = first_existing_file(brainMaskCandidates);
if brainMaskFound
    cfg.templates.standard.brainMaskNii = brainMaskPath;
end

% ---- T1 可视化模板（DPABI ch2.nii，回退 SPM canonical）----
t1TemplateCandidates = unique(nonempty({
    getfield_safe(cfg, {'templates','standard','t1TemplateNii'}), ...
    fullfile(dpabiRoot, 'Templates', 'ch2.nii'), ...
    fullfile(spmRoot,   'canonical', 'single_subj_T1.nii'), ...
    fullfile(dpabiRoot, 'Templates', 'SPMTemplates', 'canonical', 'single_subj_T1.nii')
}));
[t1TemplatePath, t1TemplateFound] = first_existing_file(t1TemplateCandidates);
if t1TemplateFound
    cfg.templates.standard.t1TemplateNii = t1TemplatePath;
    if isfield(cfg, 'visualization')
        cfg.visualization.brainTemplateNii = t1TemplatePath;
    end
end

% ---- SPM Renderer 模板 MAT（用于记录与兼容，不强依赖）----
renderMatCandidates = unique(nonempty({
    getfield_safe(cfg, {'visualization','spmRenderTemplateMat'}), ...
    fullfile(spmRoot,   'rend', 'render_single_subj.mat'), ...
    fullfile(dpabiRoot, 'Templates', 'SPMTemplates', 'rend', 'render_single_subj.mat')
}));
[renderMatPath, renderMatFound] = first_existing_file(renderMatCandidates);
if renderMatFound
    cfg.visualization.spmRenderTemplateMat = renderMatPath;
end

% ---- 记录推断结果 ----
cfg.templateResolution.dpabiRoot = dpabiRoot;
cfg.templateResolution.spmRoot   = spmRoot;
cfg.templateResolution.dartelTemplateFound  = dartelFound;
cfg.templateResolution.brainMaskFound       = brainMaskFound;
cfg.templateResolution.t1TemplateFound      = t1TemplateFound;
cfg.templateResolution.spmRenderMatFound    = renderMatFound;
cfg.templateResolution.dartelTemplate       = getfield_safe(cfg, {'templates','dartel','template4DNii'});
cfg.templateResolution.brainMask            = getfield_safe(cfg, {'templates','standard','brainMaskNii'});
cfg.templateResolution.t1Template           = getfield_safe(cfg, {'templates','standard','t1TemplateNii'});
cfg.templateResolution.spmRenderTemplateMat = getfield_safe(cfg, {'visualization','spmRenderTemplateMat'});

end

function [pathOut, found] = first_existing_file(candidates)
pathOut = '';
found = false;
for i = 1:numel(candidates)
    p = candidates{i};
    if exist(p, 'file')
        pathOut = p;
        found = true;
        return;
    end
end
if ~isempty(candidates)
    pathOut = candidates{1};
end
end

function out = nonempty(c)
out = c(~cellfun(@isempty, c));
end

function v = getfield_safe(s, fields)
v = [];
try
    v = s;
    for i = 1:numel(fields)
        if ~isstruct(v) || ~isfield(v, fields{i})
            v = [];
            return;
        end
        v = v.(fields{i});
    end
catch
    v = [];
end
if isempty(v), v = []; end
end
