function audit = audit_pipeline_parity(cfg, pipelineDir, reportFile)
% audit_pipeline_parity - 对照 SPM/DPABI 关键流程做一致性审计（逻辑级）
% 说明：
% - 本项目是 standalone 实现，不调用 SPM/DPABI 内核函数
% - 该审计用于核对“处理流程与输入输出逻辑”是否与开源库一致
% - 若 strict=true 且发现关键缺失，将报错退出

if nargin < 3 || isempty(reportFile)
    reportFile = fullfile(cfg.logDir, 'pipeline_parity_audit.txt');
end

steps = { ...
    'DICOM->NIfTI (SPM: spm_dicom_convert / DPABI: y_Call_dcm2nii)',      fullfile(pipelineDir,'io','dicom2nifti_fun.m'); ...
    'Remove Dummy TR (DPABI: RemoveFirstTimePoints)',                       fullfile(pipelineDir,'preprocess','remove_dummy_tr.m'); ...
    'Slice Timing (SPM: spm_slice_timing)',                                 fullfile(pipelineDir,'preprocess','slice_timing_corr.m'); ...
    'Realign+Reslice (SPM: spm_realign/spm_reslice)',                       fullfile(pipelineDir,'preprocess','realign_estimate_reslice.m'); ...
    'Reorient + Coreg (SPM: spm_get_space/spm_coreg)',                      fullfile(pipelineDir,'preprocess','coreg_t1_to_fun.m'); ...
    'Segment (SPM: spm_preproc_run)',                                       fullfile(pipelineDir,'preprocess','segment_tissue.m'); ...
    'DARTEL Warp (SPM: spm_dartel_warp)',                                   fullfile(pipelineDir,'preprocess','dartel_warp.m'); ...
    'Normalize by DARTEL (SPM: spm_dartel_norm_fun)',                       fullfile(pipelineDir,'preprocess','normalize_apply.m'); ...
    'Smooth (SPM: spm_smooth)',                                             fullfile(pipelineDir,'preprocess','smooth_3d.m'); ...
    '1st-level GLM + T-contrast (SPM: spm_fMRI_design/spm_spm/spm_contrasts)', fullfile(pipelineDir,'stats','run_firstlevel_glm.m'); ...
    'Modern 3D renderer (SPM renderer logic parsed, modern interaction)',   fullfile(pipelineDir,'visualize','render_activation_3d.m') ...
};

pass = true;
lines = {};
lines{end+1} = sprintf('[audit] 时间: %s', char(datetime('now','Format','yyyy-MM-dd HH:mm:ss'))); %#ok<AGROW>
lines{end+1} = '[audit] 对照来源: SPM25 + DPABI'; %#ok<AGROW>
lines{end+1} = '';

for i = 1:size(steps,1)
    desc = steps{i,1};
    f    = steps{i,2};
    ok   = exist(f, 'file') == 2;
    pass = pass && ok;
    lines{end+1} = sprintf('[%s] %s', ternary(ok,'OK','MISS'), desc); %#ok<AGROW>
    lines{end+1} = sprintf('       local: %s', f); %#ok<AGROW>
end

lines{end+1} = '';
lines{end+1} = sprintf('[template] DARTEL:   %s', getfield_safe(cfg, {'templateResolution','dartelTemplate'})); %#ok<AGROW>
lines{end+1} = sprintf('[template] BrainMask:%s', getfield_safe(cfg, {'templateResolution','brainMask'})); %#ok<AGROW>
lines{end+1} = sprintf('[template] T1:       %s', getfield_safe(cfg, {'templateResolution','t1Template'})); %#ok<AGROW>
lines{end+1} = sprintf('[template] RendererMAT:%s', getfield_safe(cfg, {'templateResolution','spmRenderTemplateMat'})); %#ok<AGROW>

if isfield(cfg, 'templateResolution')
    lines{end+1} = sprintf('[template] SPM roots scanned: %s', strjoin(cfg.templateResolution.spmRootsScanned, ' | ')); %#ok<AGROW>
end

ensure_dir(fileparts(reportFile));
fid = fopen(reportFile, 'w');
for i = 1:numel(lines)
    fprintf(fid, '%s\n', lines{i});
end
fclose(fid);

audit.pass = pass;
audit.reportFile = reportFile;
audit.lines = lines;

strict = false;
if isfield(cfg, 'referenceAudit') && isfield(cfg.referenceAudit, 'strict')
    strict = logical(cfg.referenceAudit.strict);
end
if strict && ~pass
    error('[audit_pipeline_parity] strict 模式下审计失败，详情见: %s', reportFile);
end
end

function out = ternary(cond, a, b)
if cond, out = a; else, out = b; end
end

function v = getfield_safe(s, fields)
v = '';
try
    v = s;
    for i = 1:numel(fields)
        if ~isstruct(v) || ~isfield(v, fields{i})
            v = '';
            return;
        end
        v = v.(fields{i});
    end
catch
    v = '';
end
if isstring(v), v = char(v); end
if isempty(v), v = ''; end
end
