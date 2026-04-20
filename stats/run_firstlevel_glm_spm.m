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

goldRef = load_gold_spm_reference(cfg);
spec = build_firstlevel_spec(cfg, goldRef);

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
matlabbatch{1}.spm.stats.fmri_spec.timing.units = spec.units;
matlabbatch{1}.spm.stats.fmri_spec.timing.RT = spec.RT;
matlabbatch{1}.spm.stats.fmri_spec.timing.fmri_t = spec.fmri_t;
matlabbatch{1}.spm.stats.fmri_spec.timing.fmri_t0 = spec.fmri_t0;

matlabbatch{1}.spm.stats.fmri_spec.sess.scans = scans;
for c = 1:numel(spec.condNames)
    matlabbatch{1}.spm.stats.fmri_spec.sess.cond(c).name = spec.condNames{c};
    matlabbatch{1}.spm.stats.fmri_spec.sess.cond(c).onset = double(spec.condOnsets{c}(:)');
    matlabbatch{1}.spm.stats.fmri_spec.sess.cond(c).duration = double(spec.condDurations{c}(:)');
    matlabbatch{1}.spm.stats.fmri_spec.sess.cond(c).tmod = 0;
    matlabbatch{1}.spm.stats.fmri_spec.sess.cond(c).pmod = struct('name', {}, 'param', {}, 'poly', {});
    matlabbatch{1}.spm.stats.fmri_spec.sess.cond(c).orth = spec.condOrth(c);
end

matlabbatch{1}.spm.stats.fmri_spec.sess.multi = {''};
matlabbatch{1}.spm.stats.fmri_spec.sess.regress = struct('name', {}, 'val', {});
matlabbatch{1}.spm.stats.fmri_spec.sess.multi_reg = {rpFile};
matlabbatch{1}.spm.stats.fmri_spec.sess.hpf = spec.hpf;

matlabbatch{1}.spm.stats.fmri_spec.fact = struct('name', {}, 'levels', {});
matlabbatch{1}.spm.stats.fmri_spec.bases.hrf.derivs = [0 0];
matlabbatch{1}.spm.stats.fmri_spec.volt = 1;
matlabbatch{1}.spm.stats.fmri_spec.global = spec.global;
matlabbatch{1}.spm.stats.fmri_spec.mthresh = spec.mthresh;
matlabbatch{1}.spm.stats.fmri_spec.mask = {''};
matlabbatch{1}.spm.stats.fmri_spec.cvi = spec.cvi;

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

contrastDef = build_contrast_definition(SPM, cfg, goldRef, spec.condNames);
SPM.xCon = spm_FcUtil('Set', contrastDef.name, 'T', 'c', contrastDef.vec', SPM.xX.xKXs);
SPM = spm_contrasts(SPM, 1:numel(SPM.xCon));
save(spmMat, 'SPM');

copy_named_contrast_outputs(outDir, contrastDef.name);
generate_multiple_comparison_report(outDir, 1, contrastDef.name, cfg, spmDir);

if isfield(cfg, 'visualization') && isfield(cfg.visualization, 'enable') && cfg.visualization.enable
    safeName = sanitize_name(contrastDef.name);
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

function contrastDef = build_contrast_definition(SPM, cfg, goldRef, condNames)
contrastDef = struct();

if isstruct(goldRef) && isfield(goldRef, 'available') && goldRef.available && ...
        isfield(goldRef, 'contrastVec') && ~isempty(goldRef.contrastVec)
    c = double(goldRef.contrastVec(:)');
    if numel(c) == size(SPM.xX.X, 2)
        contrastDef.vec = c;
        contrastDef.name = goldRef.contrastName;
        fprintf('[run_firstlevel_glm_spm] 使用 Gold SPM.mat 对比: %s\n', contrastDef.name);
        return;
    end
end

contrastDef.vec = build_contrast_vector(SPM, cfg, condNames);
contrastDef.name = cfg.tcons.name;
end

function contrastVec = build_contrast_vector(SPM, cfg, condNames)
nCols = size(SPM.xX.X, 2);
contrastVec = zeros(1, nCols);
weights = cfg.tcons.weight(:)';
nConds = numel(condNames);
nUse = min(numel(weights), nConds);

for i = 1:nUse
    targetName = sprintf('Sn(1) %s*bf(1)', condNames{i});
    idx = find(strcmp(SPM.xX.name, targetName), 1);
    if isempty(idx)
        idx = find(contains(SPM.xX.name, condNames{i}), 1);
    end
    if ~isempty(idx)
        contrastVec(idx) = weights(i);
    end
end

if all(abs(contrastVec) < eps)
    error('[run_firstlevel_glm_spm] 未匹配到任何条件列，无法构建对比向量');
end
end

function spec = build_firstlevel_spec(cfg, goldRef)
spec = struct();

spec.units = get_cfg_field(cfg, 'units', 'scans');
spec.RT = get_cfg_field(cfg, 'TR', 2);
spec.fmri_t = get_glm_field(cfg, 'microtimeResolution', 16);
spec.fmri_t0 = get_glm_field(cfg, 'microtimeOnsetBin', 8);
spec.hpf = get_cfg_field(cfg, 'hpf', 128);
spec.global = 'None';
spec.mthresh = 0.8;
spec.cvi = 'AR(1)';

spec.condNames = cfg.cond.names;
spec.condOnsets = cfg.cond.onsets;
spec.condDurations = cfg.cond.durations;
spec.condOrth = ones(1, numel(spec.condNames));

if ~isstruct(goldRef) || ~isfield(goldRef, 'available') || ~goldRef.available
    return;
end

spec.units = goldRef.units;
spec.RT = goldRef.RT;
spec.fmri_t = goldRef.fmri_t;
spec.fmri_t0 = goldRef.fmri_t0;
spec.hpf = goldRef.hpf;
spec.global = goldRef.global;
spec.cvi = goldRef.cvi;

if ~isempty(goldRef.condNames)
    spec.condNames = goldRef.condNames;
    spec.condOnsets = goldRef.condOnsets;
    spec.condDurations = goldRef.condDurations;
    spec.condOrth = goldRef.condOrth;
end

fprintf('[run_firstlevel_glm_spm] 已对齐 Gold SPM.mat: %s\n', goldRef.file);
end

function goldRef = load_gold_spm_reference(cfg)
goldRef = struct('available', false, 'file', '', 'units', 'scans', 'RT', 2, ...
    'fmri_t', 16, 'fmri_t0', 8, 'hpf', 128, 'global', 'None', 'cvi', 'AR(1)', ...
    'condNames', {{}}, 'condOnsets', {{}}, 'condDurations', {{}}, 'condOrth', [], ...
    'contrastName', '', 'contrastVec', []);

spmFile = resolve_gold_spm_mat(cfg);
if isempty(spmFile) || ~exist(spmFile, 'file')
    return;
end

S = load(spmFile, 'SPM');
if ~isfield(S, 'SPM') || ~isstruct(S.SPM)
    return;
end
SPM = S.SPM;

goldRef.available = true;
goldRef.file = spmFile;

if isfield(SPM, 'xBF') && isfield(SPM.xBF, 'UNITS') && ~isempty(SPM.xBF.UNITS)
    goldRef.units = char(SPM.xBF.UNITS);
end
if isfield(SPM, 'xY') && isfield(SPM.xY, 'RT') && ~isempty(SPM.xY.RT)
    goldRef.RT = double(SPM.xY.RT);
end
if isfield(SPM, 'xBF') && isfield(SPM.xBF, 'T') && ~isempty(SPM.xBF.T)
    goldRef.fmri_t = double(SPM.xBF.T);
end
if isfield(SPM, 'xBF') && isfield(SPM.xBF, 'T0') && ~isempty(SPM.xBF.T0)
    goldRef.fmri_t0 = double(SPM.xBF.T0);
end
if isfield(SPM, 'xX') && isfield(SPM.xX, 'K') && ~isempty(SPM.xX.K) && isfield(SPM.xX.K(1), 'HParam')
    goldRef.hpf = double(SPM.xX.K(1).HParam);
end
if isfield(SPM, 'xGX') && isfield(SPM.xGX, 'iGXcalc') && ~isempty(SPM.xGX.iGXcalc)
    goldRef.global = char(SPM.xGX.iGXcalc);
end
if isfield(SPM, 'xVi') && isfield(SPM.xVi, 'form') && ~isempty(SPM.xVi.form)
    goldRef.cvi = char(SPM.xVi.form);
end

if isfield(SPM, 'Sess') && ~isempty(SPM.Sess) && isfield(SPM.Sess(1), 'U')
    U = SPM.Sess(1).U;
    nU = numel(U);
    names = cell(1, nU);
    onsets = cell(1, nU);
    durations = cell(1, nU);
    orth = ones(1, nU);
    for i = 1:nU
        nm = U(i).name;
        if iscell(nm) && ~isempty(nm)
            nm = nm{1};
        end
        names{i} = char(nm);
        onsets{i} = double(U(i).ons(:)');
        durations{i} = double(U(i).dur(:)');
        if isfield(U(i), 'orth') && ~isempty(U(i).orth)
            orth(i) = double(U(i).orth(1));
        end
    end
    goldRef.condNames = names;
    goldRef.condOnsets = onsets;
    goldRef.condDurations = durations;
    goldRef.condOrth = orth;
end

if isfield(SPM, 'xCon') && ~isempty(SPM.xCon)
    idxT = find(arrayfun(@(x) isfield(x, 'STAT') && strcmpi(x.STAT, 'T'), SPM.xCon), 1);
    if isempty(idxT)
        idxT = 1;
    end
    goldRef.contrastName = char(SPM.xCon(idxT).name);
    goldRef.contrastVec = full(SPM.xCon(idxT).c(:)');
end
end

function spmFile = resolve_gold_spm_mat(cfg)
spmFile = '';

if isstruct(cfg) && isfield(cfg, 'gold') && isfield(cfg.gold, 'firstLevelSPMMat') && ...
        ischar(cfg.gold.firstLevelSPMMat) && ~isempty(cfg.gold.firstLevelSPMMat)
    if exist(cfg.gold.firstLevelSPMMat, 'file')
        spmFile = cfg.gold.firstLevelSPMMat;
        return;
    end
end

if ~(isstruct(cfg) && isfield(cfg, 'gold') && isfield(cfg.gold, 'baseDir') && ~isempty(cfg.gold.baseDir))
    return;
end

base = cfg.gold.baseDir;
subID = get_cfg_field(cfg, 'subID', 'Sub_01');
subNoUnderscore = regexprep(char(subID), '_', '');
cand = { ...
    fullfile(base, [subNoUnderscore '_1stLevel'], 'SPM.mat'), ...
    fullfile(base, [char(subID) '_1stLevel'], 'SPM.mat'), ...
    fullfile(base, 'Sub01_1stLevel', 'SPM.mat') ...
    };

for i = 1:numel(cand)
    if exist(cand{i}, 'file')
        spmFile = cand{i};
        return;
    end
end
end

function generate_multiple_comparison_report(outDir, contrastIndex, contrastName, cfg, spmDir)
p = 0.05;
if isstruct(cfg) && isfield(cfg, 'glm') && isfield(cfg.glm, 'significanceP') && ~isempty(cfg.glm.significanceP)
    p = double(cfg.glm.significanceP);
end
p = min(max(p, eps), 1-eps);

methods = {'FWE', 'FDR', 'none'};
if isstruct(cfg) && isfield(cfg, 'glm') && isfield(cfg.glm, 'multipleComparisonMethods') && ~isempty(cfg.glm.multipleComparisonMethods)
    raw = cfg.glm.multipleComparisonMethods;
    if ischar(raw) || isstring(raw)
        raw = {char(raw)};
    end
    if iscell(raw)
        methods = raw;
    end
end

lines = {};
lines{end+1} = sprintf('Contrast=%s', contrastName); %#ok<AGROW>
lines{end+1} = sprintf('p=%.6f', p); %#ok<AGROW>
lines{end+1} = sprintf('SPM=%s', spmDir); %#ok<AGROW>
lines{end+1} = ''; %#ok<AGROW>

for i = 1:numel(methods)
    m = normalize_method(methods{i});
    if isempty(m)
        continue;
    end
    try
        res = run_threshold_once(outDir, contrastIndex, contrastName, p, m);
        lines{end+1} = sprintf('[%s] u=%.6f vox=%d clusters=%d peakT=%.6f desc=%s', ...
            m, res.u, res.nVox, res.nClusters, res.peakT, res.desc); %#ok<AGROW>
    catch ME
        lines{end+1} = sprintf('[%s] ERROR: %s', m, ME.message); %#ok<AGROW>
    end
end

reportFile = fullfile(outDir, ['threshold_report_' sanitize_name(contrastName) '.txt']);
fid = fopen(reportFile, 'w');
for i = 1:numel(lines)
    fprintf(fid, '%s\n', lines{i});
end
fclose(fid);
fprintf('[run_firstlevel_glm_spm] 多重比较报告: %s\n', reportFile);
end

function res = run_threshold_once(outDir, contrastIndex, contrastName, p, method)
if strcmpi(method, 'FDR')
    spm_get_defaults('stats.topoFDR', 0);
end

xSPM = struct();
xSPM.swd = outDir;
xSPM.Ic = contrastIndex;
xSPM.Im = [];
xSPM.pm = [];
xSPM.Ex = [];
xSPM.title = contrastName;
xSPM.k = 0;
xSPM.u = p;
xSPM.thresDesc = method;
xSPM.n = 1;

[~, out] = spm_getSPM(xSPM);

res = struct();
res.u = out.u;
if isfield(out, 'Z') && ~isempty(out.Z)
    res.nVox = numel(out.Z);
    res.peakT = max(double(out.Z));
else
    res.nVox = 0;
    res.peakT = NaN;
end

if isfield(out, 'XYZ') && ~isempty(out.XYZ)
    c = spm_clusters(out.XYZ);
    if isempty(c)
        res.nClusters = 0;
    else
        res.nClusters = numel(unique(c(:)));
    end
else
    res.nClusters = 0;
end

if isfield(out, 'thresDesc')
    res.desc = char(out.thresDesc);
else
    res.desc = method;
end
end

function method = normalize_method(v)
method = '';
if isstring(v)
    v = char(v);
end
if ~ischar(v) || isempty(v)
    return;
end
t = lower(strtrim(v));
switch t
    case 'fwe'
        method = 'FWE';
    case 'fdr'
        method = 'FDR';
    case {'none', 'unc', 'uncorrected'}
        method = 'none';
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

function value = get_cfg_field(cfg, fieldName, defaultValue)
value = defaultValue;
if isstruct(cfg) && isfield(cfg, fieldName)
    value = cfg.(fieldName);
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
