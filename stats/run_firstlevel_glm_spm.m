function run_firstlevel_glm_spm(smoothFile, rpFile, outDir, cfg)
% run_firstlevel_glm_spm - 使用 SPM 标准流程执行一阶 GLM

fprintf('[run_firstlevel_glm_spm] 开始 SPM 一阶 GLM\n');

spmDir = resolve_spm_dir(cfg);
if ~exist(fullfile(spmDir, 'spm.m'), 'file')
    error('[run_firstlevel_glm_spm] 未找到 SPM: %s', spmDir);
end

addpath(spmDir);
spm('defaults', 'FMRI');
spm_jobman('initcfg');

ensure_dir(outDir);
cleanup_previous_outputs(outDir);

V = spm_vol(smoothFile);
nScans = numel(V);
if nScans < 2
    error('[run_firstlevel_glm_spm] 输入不是有效4D功能像: %s', smoothFile);
end

scans = cell(nScans, 1);
for i = 1:nScans
    scans{i} = sprintf('%s,%d', smoothFile, i);
end

matlabbatch = {};
matlabbatch{1}.spm.stats.fmri_spec.dir = {outDir};
matlabbatch{1}.spm.stats.fmri_spec.timing.units = cfg.units;
matlabbatch{1}.spm.stats.fmri_spec.timing.RT = cfg.TR;
matlabbatch{1}.spm.stats.fmri_spec.timing.fmri_t = get_glm_field(cfg, 'microtimeResolution', 16);
matlabbatch{1}.spm.stats.fmri_spec.timing.fmri_t0 = get_glm_field(cfg, 'microtimeOnsetBin', 8);

matlabbatch{1}.spm.stats.fmri_spec.sess.scans = scans;
for c = 1:numel(cfg.cond.names)
    matlabbatch{1}.spm.stats.fmri_spec.sess.cond(c).name = cfg.cond.names{c};
    matlabbatch{1}.spm.stats.fmri_spec.sess.cond(c).onset = double(cfg.cond.onsets{c}(:)');
    matlabbatch{1}.spm.stats.fmri_spec.sess.cond(c).duration = double(cfg.cond.durations{c}(:)');
    matlabbatch{1}.spm.stats.fmri_spec.sess.cond(c).tmod = 0;
    matlabbatch{1}.spm.stats.fmri_spec.sess.cond(c).pmod = struct('name', {}, 'param', {}, 'poly', {});
    matlabbatch{1}.spm.stats.fmri_spec.sess.cond(c).orth = 1;
end

matlabbatch{1}.spm.stats.fmri_spec.sess.multi = {''};
matlabbatch{1}.spm.stats.fmri_spec.sess.regress = struct('name', {}, 'val', {});
matlabbatch{1}.spm.stats.fmri_spec.sess.multi_reg = {rpFile};
matlabbatch{1}.spm.stats.fmri_spec.sess.hpf = cfg.hpf;

matlabbatch{1}.spm.stats.fmri_spec.fact = struct('name', {}, 'levels', {});
matlabbatch{1}.spm.stats.fmri_spec.bases.hrf.derivs = [0 0];
matlabbatch{1}.spm.stats.fmri_spec.volt = 1;
matlabbatch{1}.spm.stats.fmri_spec.global = 'None';
matlabbatch{1}.spm.stats.fmri_spec.mthresh = 0.8;
matlabbatch{1}.spm.stats.fmri_spec.mask = {''};
matlabbatch{1}.spm.stats.fmri_spec.cvi = 'AR(1)';

spmMat = fullfile(outDir, 'SPM.mat');
matlabbatch{2}.spm.stats.fmri_est.spmmat = {spmMat};
matlabbatch{2}.spm.stats.fmri_est.write_residuals = 0;
matlabbatch{2}.spm.stats.fmri_est.method.Classical = 1;

fprintf('[run_firstlevel_glm_spm] 运行 SPM 设计与估计\n');
spm_jobman('run', matlabbatch);

if ~exist(spmMat, 'file')
    error('[run_firstlevel_glm_spm] 未生成 SPM.mat: %s', spmMat);
end

S = load(spmMat, 'SPM');
SPM = S.SPM;

contrastVec = build_contrast_vector(SPM, cfg);
SPM.xCon = spm_FcUtil('Set', cfg.tcons.name, 'T', 'c', contrastVec', SPM.xX.xKXs);
SPM = spm_contrasts(SPM, 1:numel(SPM.xCon));
save(spmMat, 'SPM');

copy_named_contrast_outputs(outDir, cfg.tcons.name);

if isfield(cfg, 'visualization') && isfield(cfg.visualization, 'enable') && cfg.visualization.enable
    safeName = sanitize_name(cfg.tcons.name);
    tFile = fullfile(outDir, ['spmT_' safeName '.nii']);
    if exist(tFile, 'file')
        try
            render_activation_3d(tFile, cfg.visualization.brainTemplateNii, outDir, cfg.visualization);
        catch ME
            warning('[run_firstlevel_glm_spm] 可视化渲染失败: %s', ME.message);
        end
    end
end

fprintf('[run_firstlevel_glm_spm] 完成，输出目录: %s\n', outDir);
end

function contrastVec = build_contrast_vector(SPM, cfg)
nCols = size(SPM.xX.X, 2);
contrastVec = zeros(1, nCols);
weights = cfg.tcons.weight(:)';
nConds = numel(cfg.cond.names);
nUse = min(numel(weights), nConds);

for i = 1:nUse
    targetName = sprintf('Sn(1) %s*bf(1)', cfg.cond.names{i});
    idx = find(strcmp(SPM.xX.name, targetName), 1);
    if isempty(idx)
        idx = find(contains(SPM.xX.name, cfg.cond.names{i}), 1);
    end
    if ~isempty(idx)
        contrastVec(idx) = weights(i);
    end
end

if all(abs(contrastVec) < eps)
    error('[run_firstlevel_glm_spm] 未匹配到任何条件列，无法构建对比向量');
end
end

function copy_named_contrast_outputs(outDir, contrastName)
safeName = sanitize_name(contrastName);

spmTsrc = fullfile(outDir, 'spmT_0001.nii');
spmTdst = fullfile(outDir, ['spmT_' safeName '.nii']);
if exist(spmTsrc, 'file')
    copyfile(spmTsrc, spmTdst);
end

conSrc = fullfile(outDir, 'con_0001.nii');
conDst = fullfile(outDir, ['con_' safeName '.nii']);
if exist(conSrc, 'file')
    copyfile(conSrc, conDst);
end
end

function cleanup_previous_outputs(outDir)
patterns = { ...
    'SPM.mat', 'beta_*.nii', 'ResMS.nii', 'mask.nii', ...
    'RPV.nii', 'con_*.nii', 'spmT_*.nii', 'spmF_*.nii', ...
    'ess_*.nii', 'ResI_*.nii', 'design_matrix.png'};

for i = 1:numel(patterns)
    files = dir(fullfile(outDir, patterns{i}));
    for j = 1:numel(files)
        f = fullfile(files(j).folder, files(j).name);
        if exist(f, 'file')
            delete(f);
        end
    end
end
end

function value = get_glm_field(cfg, fieldName, defaultValue)
value = defaultValue;
if isstruct(cfg) && isfield(cfg, 'glm') && isfield(cfg.glm, fieldName)
    value = cfg.glm.(fieldName);
end
end

function safeName = sanitize_name(name)
safeName = regexprep(char(name), '[^A-Za-z0-9_]', '_');
if isempty(safeName)
    safeName = 'contrast';
end
end

function spmDir = resolve_spm_dir(cfg)
if isstruct(cfg) && isfield(cfg, 'spm') && isfield(cfg.spm, 'dir') && ~isempty(cfg.spm.dir)
    spmDir = cfg.spm.dir;
elseif isstruct(cfg) && isfield(cfg, 'installPaths') && isfield(cfg.installPaths, 'spmRoot')
    spmDir = cfg.installPaths.spmRoot;
else
    spmDir = 'D:/spm';
end
end
