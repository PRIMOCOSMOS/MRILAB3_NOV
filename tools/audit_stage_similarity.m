function results = audit_stage_similarity()
% audit_stage_similarity - 全流程逐步对照 Gold Standard 的相似度审计
% 说明:
% 1) 覆盖核心步骤产物（含 Step04 头动参数与 Step09 DARTEL 流场）
% 2) 对于 Gold 无可比文件的步骤，显式给出 unavailable/proxy 原因
% 3) 输出 FIRST_FAIL_DIRECT 仅统计 direct 可比步骤

addpath(fileparts(fileparts(mfilename('fullpath'))));
addpath(fullfile(fileparts(fileparts(mfilename('fullpath'))), 'utils'));

thrCorr = 0.90;
thrRpCorr = 0.90;

oursBase = 'D:/MRI_PRO/MRILAB3/BOLDCODE/BOLDDATA';
goldBase = resolve_gold_base_dir(oursBase);

fprintf('=== Stage Similarity Audit (Full Step Coverage) ===\n');
fprintf('THR_CORR=%.2f, THR_RP_CORR=%.2f\n', thrCorr, thrRpCorr);
fprintf('OURS_BASE=%s\n', oursBase);
fprintf('GOLD_BASE=%s\n', goldBase);

results = struct('name', {}, 'status', {}, 'metric', {}, 'value', {}, ...
    'threshold', {}, 'pass', {}, 'gating', {}, 'note', {}, 'ours', {}, 'gold', {});

% Step01: DICOM->NIfTI (Fun 原始阶段在 Gold 中无 NIfTI，对齐缺失)
results(end+1) = make_unavailable_result( ...
    'Step01_FunImg_raw', ...
    fullfile(oursBase, 'FunImg', 'Sub_01', 'bold_4d.nii'), ...
    fullfile(goldBase, 'FunImg', 'Sub_01', 'Sub_01_BOLD_lefthand_20190710154350_10.nii'), ...
    'Gold FunImg 阶段仅有 JSON，无可比 NIfTI');

results(end+1) = evaluate_nifti_corr_step( ...
    'Step01_T1Img_raw', ...
    fullfile(oursBase, 'T1Img', 'Sub_01', 't1.nii'), ...
    fullfile(goldBase, 'T1Img', 'Sub_01', 'Sub_01_t1_mprage_sag_p2_iso_20190710154350_5a.nii'), ...
    thrCorr, true, 'direct');

% Step02: Remove Dummy（Gold 无单独目录，使用 FunImgA 作为 proxy）
results(end+1) = evaluate_nifti_corr_step( ...
    'Step02_RemoveDummy_proxyA', ...
    fullfile(oursBase, 'FunImgA', 'Sub_01', 'abold_4d.nii'), ...
    fullfile(goldBase, 'FunImgA', 'Sub_01', 'aSub_01_BOLD_lefthand_20190710154350_10.nii'), ...
    thrCorr, false, 'proxy: Gold 无独立 RemoveDummy 产物');

% Step03: Slice Timing
results(end+1) = evaluate_nifti_corr_step( ...
    'Step03_SliceTiming', ...
    fullfile(oursBase, 'FunImgA', 'Sub_01', 'stabold_4d.nii'), ...
    fullfile(goldBase, 'FunImgA', 'Sub_01', 'aSub_01_BOLD_lefthand_20190710154350_10.nii'), ...
    thrCorr, true, 'direct');

% Step04: Realign 图像 + rp 参数
results(end+1) = evaluate_nifti_corr_step( ...
    'Step04_Realign_FunAR', ...
    fullfile(oursBase, 'FunImgAR', 'Sub_01', 'rstabold_4d.nii'), ...
    fullfile(goldBase, 'FunImgAR', 'Sub_01', 'raSub_01_BOLD_lefthand_20190710154350_10.nii'), ...
    thrCorr, false, 'proxy: Gold AR 产物可能已包含 reorient 几何');

results(end+1) = evaluate_rp_step( ...
    'Step04_Realign_RP', ...
    fullfile(oursBase, 'RealignParameter', 'Sub_01', 'rp_Sub_01.txt'), ...
    fullfile(goldBase, 'RealignParameter', 'Sub_01', 'rp_aSub_01_BOLD_lefthand_20190710154350_10.txt'), ...
    thrRpCorr, true);

% Step05: Reorient
results(end+1) = evaluate_nifti_corr_step( ...
    'Step05_Reorient_Fun', ...
    fullfile(oursBase, 'FunImgAR', 'Sub_01', 'reorient_rstabold_4d.nii'), ...
    fullfile(goldBase, 'FunImgAR', 'Sub_01', 'raSub_01_BOLD_lefthand_20190710154350_10.nii'), ...
    thrCorr, true, 'direct');

results(end+1) = evaluate_nifti_corr_step( ...
    'Step05_Reorient_T1_proxy', ...
    fullfile(oursBase, 'T1Img', 'Sub_01', 'reorient_t1.nii'), ...
    fullfile(goldBase, 'T1Img', 'Sub_01', 'Sub_01_t1_mprage_sag_p2_iso_20190710154350_5a.nii'), ...
    thrCorr, false, 'proxy: Gold 无显式 reorient 文件');

% Step06: BET
results(end+1) = evaluate_bet_step( ...
    'Step06_BET', ...
    fullfile(oursBase, 'T1ImgBet', 'Sub_01', 'bet_reorient_t1.nii'), ...
    fullfile(goldBase, 'T1ImgBet', 'Sub_01', 'Bet_Sub_01_t1_mprage_sag_p2_iso_20190710154350_5_Crop_1.nii'), ...
    thrCorr, true);

% Step07: Coreg T1->Fun
results(end+1) = evaluate_nifti_corr_step( ...
    'Step07_Coreg', ...
    fullfile(oursBase, 'T1ImgCoreg', 'Sub_01', 'coreg_reorient_t1.nii'), ...
    fullfile(goldBase, 'T1ImgCoreg', 'Sub_01', 'Sub_01_t1_mprage_sag_p2_iso_20190710154350_5a_Crop_1.nii'), ...
    thrCorr, false, 'proxy: coreg.mode=identity 下仅头信息透传');

% Step08: New Segment (c1/c2/c3)
results(end+1) = evaluate_nifti_corr_step_with_samegrid( ...
    'Step08_Segment_c1', ...
    fullfile(oursBase, 'T1ImgNewSegment', 'Sub_01', 'c1t1.nii'), ...
    fullfile(goldBase, 'T1ImgNewSegment', 'Sub_01', 'c1Sub_01_t1_mprage_sag_p2_iso_20190710154350_5_Crop_1.nii'), ...
    thrCorr, true, 'direct');

results(end+1) = evaluate_nifti_corr_step_with_samegrid( ...
    'Step08_Segment_c2', ...
    fullfile(oursBase, 'T1ImgNewSegment', 'Sub_01', 'c2t1.nii'), ...
    fullfile(goldBase, 'T1ImgNewSegment', 'Sub_01', 'c2Sub_01_t1_mprage_sag_p2_iso_20190710154350_5_Crop_1.nii'), ...
    thrCorr, true, 'direct');

results(end+1) = evaluate_nifti_corr_step_with_samegrid( ...
    'Step08_Segment_c3', ...
    fullfile(oursBase, 'T1ImgNewSegment', 'Sub_01', 'c3t1.nii'), ...
    fullfile(goldBase, 'T1ImgNewSegment', 'Sub_01', 'c3Sub_01_t1_mprage_sag_p2_iso_20190710154350_5_Crop_1.nii'), ...
    thrCorr, true, 'direct');

% Step09: DARTEL flow (5D u_rc1*)
results(end+1) = evaluate_flowfield_step( ...
    'Step09_DARTEL_Flow', ...
    fullfile(oursBase, 'T1ImgNewSegment', 'Sub_01', 'u_rc1t1_Template.nii'), ...
    fullfile(goldBase, 'T1ImgNewSegment', 'Sub_01', 'u_rc1Sub_01_t1_mprage_sag_p2_iso_20190710154350_5_Crop_1_Template.nii'), ...
    thrCorr, true);

% Step10: Normalize ARW
results(end+1) = evaluate_nifti_corr_step( ...
    'Step10_Normalize_ARW', ...
    fullfile(oursBase, 'FunImgARW', 'Sub_01', 'wreorient_rstabold_4d.nii'), ...
    fullfile(goldBase, 'FunImgARW', 'Sub_01', 'wraSub_01_BOLD_lefthand_20190710154350_10.nii'), ...
    thrCorr, true, 'direct');

% Step11: Smooth ARWS
results(end+1) = evaluate_nifti_corr_step( ...
    'Step11_Smooth_ARWS', ...
    fullfile(oursBase, 'FunImgARWS', 'Sub_01', 'swreorient_rstabold_4d.nii'), ...
    fullfile(goldBase, 'FunImgARWS', 'Sub_01', 'swraSub_01_BOLD_lefthand_20190710154350_10.nii'), ...
    thrCorr, true, 'direct');

% Step12: First-level T-contrast
results(end+1) = evaluate_nifti_corr_step_multiours( ...
    'Step12_GLM_T', ...
    { ...
        fullfile(oursBase, 'Sub_01_1stLevel', 'spmT_0001.nii'), ...
        fullfile(oursBase, 'Sub_01_1stLevel', 'spmT_t_contrast.nii'), ...
        fullfile(oursBase, 'Sub_01_1stLevel', 'spmT_Task_gt_Baseline.nii') ...
    }, ...
    fullfile(goldBase, 'Sub01_1stLevel', 'spmT_0001.nii'), ...
    thrCorr, true, 'direct');

print_results(results);
firstFail = find_first_fail_direct(results);
if isempty(firstFail)
    fprintf('FIRST_FAIL_DIRECT none\n');
else
    fprintf('FIRST_FAIL_DIRECT %s\n', firstFail);
end

end

function result = evaluate_nifti_corr_step(name, oursFile, goldFile, thr, gating, modeNote)
result = init_result(name, oursFile, goldFile, gating);
result.metric = 'corr';
result.threshold = thr;
result.note = modeNote;

if ~exist(oursFile, 'file')
    result.status = 'missing_ours';
    result.pass = 0;
    return;
end
if ~exist(goldFile, 'file')
    result.status = 'missing_gold';
    result.pass = 0;
    return;
end

c = corr_resampled(oursFile, goldFile);
result.value = c;
result.status = 'ok';
result.pass = isfinite(c) && (c >= thr);
end

function result = evaluate_nifti_corr_step_multiours(name, oursFiles, goldFile, thr, gating, modeNote)
if ~iscell(oursFiles)
    oursFiles = {oursFiles};
end

oursFile = '';
for i = 1:numel(oursFiles)
    if ischar(oursFiles{i}) && exist(oursFiles{i}, 'file')
        oursFile = oursFiles{i};
        break;
    end
end

if isempty(oursFile)
    oursFile = oursFiles{1};
end

result = evaluate_nifti_corr_step(name, oursFile, goldFile, thr, gating, modeNote);
[~, oursName, oursExt] = fileparts(oursFile);
result.note = sprintf('%s ours=%s%s', modeNote, oursName, oursExt);
end

function result = evaluate_nifti_corr_step_with_samegrid(name, oursFile, goldFile, thr, gating, modeNote)
result = init_result(name, oursFile, goldFile, gating);
result.metric = 'corr_resampled';
result.threshold = thr;
result.note = modeNote;

if ~exist(oursFile, 'file')
    result.status = 'missing_ours';
    result.pass = 0;
    return;
end
if ~exist(goldFile, 'file')
    result.status = 'missing_gold';
    result.pass = 0;
    return;
end

cResampled = corr_resampled(oursFile, goldFile);
cVoxel = corr_voxel_same_grid(oursFile, goldFile);
if isfinite(cVoxel)
    cEff = max(cResampled, cVoxel);
    result.note = sprintf('%s voxel=%.6f', modeNote, cVoxel);
else
    cEff = cResampled;
end

result.value = cEff;
result.status = 'ok';
result.pass = isfinite(cEff) && (cEff >= thr);
end

function result = evaluate_bet_step(name, oursFile, goldFile, thr, gating)
result = init_result(name, oursFile, goldFile, gating);
result.metric = 'corr_resampled';
result.threshold = thr;

if ~exist(oursFile, 'file')
    result.status = 'missing_ours';
    result.pass = 0;
    return;
end
if ~exist(goldFile, 'file')
    result.status = 'missing_gold';
    result.pass = 0;
    return;
end

cResampled = corr_resampled(oursFile, goldFile);
cVoxel = corr_voxel_same_grid(oursFile, goldFile);
if isfinite(cVoxel)
    cEff = max(cResampled, cVoxel);
    result.note = sprintf('voxel=%.6f', cVoxel);
else
    cEff = cResampled;
    result.note = 'voxel=nan';
end

result.value = cEff;
result.status = 'ok';
result.pass = isfinite(cEff) && (cEff >= thr);
end

function result = evaluate_rp_step(name, oursFile, goldFile, thr, gating)
result = init_result(name, oursFile, goldFile, gating);
result.metric = 'corr_mean';
result.threshold = thr;

if ~exist(oursFile, 'file')
    result.status = 'missing_ours';
    result.pass = 0;
    return;
end
if ~exist(goldFile, 'file')
    result.status = 'missing_gold';
    result.pass = 0;
    return;
end

A = load_numeric_matrix(oursFile);
B = load_numeric_matrix(goldFile);
if size(A,2) < 6 || size(B,2) < 6
    result.status = 'invalid_shape';
    result.note = sprintf('ours=[%d,%d], gold=[%d,%d]', size(A,1), size(A,2), size(B,1), size(B,2));
    result.pass = 0;
    return;
end

n = min(size(A,1), size(B,1));
A = double(A(1:n,1:6));
B = double(B(1:n,1:6));

colCorr = nan(1, 6);
for k = 1:6
    a = A(:,k);
    b = B(:,k);
    if std(a) < eps && std(b) < eps
        colCorr(k) = 1;
    else
        cc = corrcoef(a, b);
        colCorr(k) = cc(1,2);
    end
end

corrMean = mean(colCorr, 'omitnan');
mae = mean(abs(A - B), 'all');

result.value = corrMean;
result.status = 'ok';
result.pass = isfinite(corrMean) && (corrMean >= thr);
result.note = sprintf('mae=%.6g colCorr=[%.4f %.4f %.4f %.4f %.4f %.4f]', ...
    mae, colCorr(1), colCorr(2), colCorr(3), colCorr(4), colCorr(5), colCorr(6));
end

function result = evaluate_flowfield_step(name, oursFile, goldFile, thr, gating)
result = init_result(name, oursFile, goldFile, gating);
result.metric = 'corr_flow_focus';
result.threshold = thr;

if ~exist(oursFile, 'file')
    result.status = 'missing_ours';
    result.pass = 0;
    return;
end
if ~exist(goldFile, 'file')
    result.status = 'missing_gold';
    result.pass = 0;
    return;
end

try
    spmDir = resolve_spm_dir_for_audit();
    if ~exist(fullfile(spmDir, 'spm.m'), 'file')
        error('未找到 SPM: %s', spmDir);
    end
    addpath(spmDir);

    N1 = nifti(oursFile);
    N2 = nifti(goldFile);

    d1 = double(N1.dat.dim);
    d2 = double(N2.dat.dim);
    if numel(d1) < 5 || numel(d2) < 5 || d1(5) < 3 || d2(5) < 3
        error('流场维度异常: ours=%s, gold=%s', mat2str(d1), mat2str(d2));
    end

    compA = cell(1,3);
    compB = cell(1,3);
    for c = 1:3
        A = squeeze(double(N1.dat(:,:,:,:,c)));
        B = squeeze(double(N2.dat(:,:,:,:,c)));

        if ~isequal(size(A), size(B)) || any(abs(N1.mat(:) - N2.mat(:)) > 1e-6)
            A = resample_vol_affine(A, N1.mat, N2.mat, size(B));
        end

        compA{c} = A;
        compB{c} = B;
    end

    magA = sqrt(compA{1}.^2 + compA{2}.^2 + compA{3}.^2);
    magB = sqrt(compB{1}.^2 + compB{2}.^2 + compB{3}.^2);

    baseMask = isfinite(magA) & isfinite(magB) & ((magA ~= 0) | (magB ~= 0));
    if ~any(baseMask(:))
        result.status = 'invalid_mask';
        result.pass = 0;
        result.note = 'flow mask is empty';
        return;
    end

    focusThr = prctile(magB(baseMask), 70);
    focusMask = baseMask & ((magA >= focusThr) | (magB >= focusThr));
    if nnz(focusMask) < 1000
        focusMask = baseMask;
    end

    compCorr = nan(1,3);
    compMae = nan(1,3);
    for c = 1:3
        A = compA{c};
        B = compB{c};
        m = focusMask & isfinite(A) & isfinite(B);
        if any(m(:))
            cc = corrcoef(A(m), B(m));
            compCorr(c) = cc(1,2);
            compMae(c) = mean(abs(A(m) - B(m)), 'omitnan');
        end
    end

    flowCorr = mean(compCorr, 'omitnan');
    flowMae = mean(compMae, 'omitnan');
    result.value = flowCorr;
    result.status = 'ok';
    result.pass = isfinite(flowCorr) && (flowCorr >= thr);
    result.note = sprintf('focus=%d/%d mae=%.6g compCorr=[%.4f %.4f %.4f] compMae=[%.6g %.6g %.6g]', ...
        nnz(focusMask), nnz(baseMask), flowMae, compCorr(1), compCorr(2), compCorr(3), ...
        compMae(1), compMae(2), compMae(3));
catch ME
    result.status = 'error';
    result.pass = 0;
    result.note = ME.message;
end
end

function result = make_unavailable_result(name, oursFile, goldFile, note)
result = init_result(name, oursFile, goldFile, false);
result.status = 'unavailable';
result.metric = 'none';
result.value = NaN;
result.threshold = NaN;
result.pass = NaN;
result.note = note;
end

function result = init_result(name, oursFile, goldFile, gating)
result = struct();
result.name = name;
result.status = 'init';
result.metric = 'none';
result.value = NaN;
result.threshold = NaN;
result.pass = NaN;
result.gating = logical(gating);
result.note = '';
result.ours = oursFile;
result.gold = goldFile;
end

function print_results(results)
fprintf('--- Stage Results ---\n');
for i = 1:numel(results)
    r = results(i);
    passTxt = pass_to_text(r.pass);
    fprintf(['AUDIT_STEP %-24s status=%-12s metric=%-14s value=%10.6f ', ...
             'thr=%7.3f pass=%-3s gating=%d note="%s"\n'], ...
            r.name, r.status, r.metric, r.value, r.threshold, passTxt, r.gating, r.note);
end
end

function t = pass_to_text(v)
if isnan(v)
    t = 'na';
elseif v
    t = 'yes';
else
    t = 'no';
end
end

function name = find_first_fail_direct(results)
name = '';
for i = 1:numel(results)
    r = results(i);
    if r.gating && ~isnan(r.pass) && ~logical(r.pass)
        name = r.name;
        return;
    end
end
end

function M = load_numeric_matrix(filePath)
try
    M = load(filePath);
catch
    M = dlmread(filePath);
end
if isstruct(M)
    f = fieldnames(M);
    M = M.(f{1});
end
end

function c = corr_resampled(oursFile, goldFile)
[A, Ah] = nifti_read(oursFile);
[B, Bh] = nifti_read(goldFile);

if ndims(A) == 3
    A = reshape(A, size(A,1), size(A,2), size(A,3), 1);
end
if ndims(B) == 3
    B = reshape(B, size(B,1), size(B,2), size(B,3), 1);
end

nt = min(size(A,4), size(B,4));
Ar = zeros([size(B,1), size(B,2), size(B,3), nt], 'double');
for t = 1:nt
    Ar(:,:,:,t) = resample_vol_affine(double(A(:,:,:,t)), Ah.affine, Bh.affine, size(B(:,:,:,1)));
end
Bb = double(B(:,:,:,1:nt));

m = isfinite(Ar) & isfinite(Bb) & ((Ar ~= 0) | (Bb ~= 0));
if ~any(m(:))
    c = NaN;
    return;
end
cc = corrcoef(Ar(m), Bb(m));
c = cc(1,2);
end

function c = corr_voxel_same_grid(oursFile, goldFile)
[A, ~] = nifti_read(oursFile);
[B, ~] = nifti_read(goldFile);

if ndims(A) == 3
    A = reshape(A, size(A,1), size(A,2), size(A,3), 1);
end
if ndims(B) == 3
    B = reshape(B, size(B,1), size(B,2), size(B,3), 1);
end

if any(size(A,1:3) ~= size(B,1:3))
    c = NaN;
    return;
end

nt = min(size(A,4), size(B,4));
A = double(A(:,:,:,1:nt));
B = double(B(:,:,:,1:nt));

m = isfinite(A) & isfinite(B) & ((A ~= 0) | (B ~= 0));
if ~any(m(:))
    c = NaN;
    return;
end

cc = corrcoef(A(m), B(m));
c = cc(1,2);
end

function goldBase = resolve_gold_base_dir(oursBase)
candidates = { ...
    oursBase, ...
    'D:/MRI_PRO/MRILAB3/BOLDDATA'};

refRel = fullfile('FunImgARW', 'Sub_01', 'wraSub_01_BOLD_lefthand_20190710154350_10.nii');
goldBase = '';
for i = 1:numel(candidates)
    cand = candidates{i};
    if exist(fullfile(cand, refRel), 'file')
        goldBase = cand;
        return;
    end
end

error('未找到 Gold 数据根目录，请检查候选路径: %s | %s', candidates{1}, candidates{2});
end

function spmDir = resolve_spm_dir_for_audit()
candidates = {'D:/spm', 'D:/Coding2', 'D:/spm25'};
spmDir = '';
for i = 1:numel(candidates)
    cand = candidates{i};
    if exist(fullfile(cand, 'spm.m'), 'file')
        spmDir = cand;
        return;
    end
end
error('未找到可用 SPM 路径，请检查 D:/spm 或 D:/Coding2');
end
