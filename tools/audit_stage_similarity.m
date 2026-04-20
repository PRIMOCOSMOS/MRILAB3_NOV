function audit_stage_similarity()
% audit_stage_similarity - Stage-by-stage similarity audit against GOLD.
% Prints direct correlation and best flip-only correlation for each stage.

addpath(fileparts(fileparts(mfilename('fullpath'))));
addpath(fullfile(fileparts(fileparts(mfilename('fullpath'))), 'utils'));

thr = 0.90;

oursBase = 'D:/MRI_PRO/MRILAB3/BOLDCODE/BOLDDATA';
goldBase = resolve_gold_base_dir(oursBase);

steps = {
    'A', ...
    fullfile(oursBase,'FunImgA','Sub_01','abold_4d.nii'), ...
    fullfile(goldBase,'FunImgA','Sub_01','aSub_01_BOLD_lefthand_20190710154350_10.nii');
    'AR', ...
    fullfile(oursBase,'FunImgAR','Sub_01','reorient_rstabold_4d.nii'), ...
    fullfile(goldBase,'FunImgAR','Sub_01','raSub_01_BOLD_lefthand_20190710154350_10.nii');
    'T1_raw', ...
    fullfile(oursBase,'T1Img','Sub_01','reorient_t1.nii'), ...
    fullfile(goldBase,'T1Img','Sub_01','Sub_01_t1_mprage_sag_p2_iso_20190710154350_5a.nii');
    'T1_crop_ref', ...
    fullfile(oursBase,'T1Img','Sub_01','reorient_t1.nii'), ...
    fullfile(goldBase,'T1Img','Sub_01','Sub_01_t1_mprage_sag_p2_iso_20190710154350_5a_Crop_1.nii');
    'BET', ...
    fullfile(oursBase,'T1ImgBet','Sub_01','bet_reorient_t1.nii'), ...
    fullfile(goldBase,'T1ImgBet','Sub_01','Bet_Sub_01_t1_mprage_sag_p2_iso_20190710154350_5_Crop_1.nii');
    'Coreg', ...
    fullfile(oursBase,'T1ImgCoreg','Sub_01','coreg_reorient_t1.nii'), ...
    fullfile(goldBase,'T1ImgCoreg','Sub_01','Sub_01_t1_mprage_sag_p2_iso_20190710154350_5a_Crop_1.nii');
    'ARW', ...
    fullfile(oursBase,'FunImgARW','Sub_01','wreorient_rstabold_4d.nii'), ...
    fullfile(goldBase,'FunImgARW','Sub_01','wraSub_01_BOLD_lefthand_20190710154350_10.nii');
    'ARWS', ...
    fullfile(oursBase,'FunImgARWS','Sub_01','swreorient_rstabold_4d.nii'), ...
    fullfile(goldBase,'FunImgARWS','Sub_01','swraSub_01_BOLD_lefthand_20190710154350_10.nii');
    'T', ...
    fullfile(oursBase,'Sub_01_1stLevel','spmT_Task_gt_Baseline.nii'), ...
    fullfile(goldBase,'Sub01_1stLevel','spmT_0001.nii')};

fprintf('=== Stage Similarity Audit (threshold=%.2f) ===\n', thr);
fprintf('OURS_BASE=%s\n', oursBase);
fprintf('GOLD_BASE=%s\n', goldBase);
firstFail = '';
for i = 1:size(steps, 1)
    name = steps{i, 1};
    ours = steps{i, 2};
    gold = steps{i, 3};

    if ~exist(ours, 'file') || ~exist(gold, 'file')
        fprintf('AUDIT_%s missing\n', name);
        if isempty(firstFail)
            firstFail = [name ' (missing)'];
        end
        continue;
    end

    cDirect = corr_resampled(ours, gold);
    if strcmpi(name, 'BET')
        cVoxel = corr_voxel_same_grid(ours, gold);
        cEff = cDirect;
        if isfinite(cVoxel) && (~isfinite(cEff) || cVoxel > cEff)
            cEff = cVoxel;
        end
        pass = cEff >= thr;
        fprintf('AUDIT_%-12s direct=%.6f voxel=%.6f effective=%.6f pass=%d\n', name, cDirect, cVoxel, cEff, pass);
    else
        cBestFlip = best_flip_corr_firstvol(ours, gold);
        pass = cDirect >= thr;
        fprintf('AUDIT_%-12s direct=%.6f bestFlip=%.6f pass=%d\n', name, cDirect, cBestFlip, pass);
    end

    if ~pass && isempty(firstFail)
        firstFail = name;
    end
end

fprintf('--- Tissue Segmentation Matrix (ours c1/c2/c3 vs gold c1/c2/c3) ---\n');
oursSeg = {
    fullfile(oursBase,'T1ImgNewSegment','Sub_01','c1t1.nii'), ...
    fullfile(oursBase,'T1ImgNewSegment','Sub_01','c2t1.nii'), ...
    fullfile(oursBase,'T1ImgNewSegment','Sub_01','c3t1.nii')};
goldSeg = {
    fullfile(goldBase,'T1ImgNewSegment','Sub_01','c1Sub_01_t1_mprage_sag_p2_iso_20190710154350_5_Crop_1.nii'), ...
    fullfile(goldBase,'T1ImgNewSegment','Sub_01','c2Sub_01_t1_mprage_sag_p2_iso_20190710154350_5_Crop_1.nii'), ...
    fullfile(goldBase,'T1ImgNewSegment','Sub_01','c3Sub_01_t1_mprage_sag_p2_iso_20190710154350_5_Crop_1.nii')};

C = nan(3, 3);
for i = 1:3
    for j = 1:3
        if exist(oursSeg{i}, 'file') && exist(goldSeg{j}, 'file')
            C(i, j) = corr_resampled(oursSeg{i}, goldSeg{j});
        end
    end
end
disp(C);

if isempty(firstFail)
    fprintf('FIRST_FAIL none\n');
else
    fprintf('FIRST_FAIL %s\n', firstFail);
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

function c = best_flip_corr_firstvol(oursFile, goldFile)
[A, Ah] = nifti_read(oursFile);
[B, Bh] = nifti_read(goldFile);

if ndims(A) == 4
    A = A(:,:,:,1);
end
if ndims(B) == 4
    B = B(:,:,:,1);
end

Ar = resample_vol_affine(double(A), Ah.affine, Bh.affine, size(B));
B = double(B);

best = -Inf;
for fx = 0:1
    for fy = 0:1
        for fz = 0:1
            T = Ar;
            if fx
                T = T(end:-1:1,:,:);
            end
            if fy
                T = T(:,end:-1:1,:);
            end
            if fz
                T = T(:,:,end:-1:1);
            end

            m = isfinite(T) & isfinite(B) & ((T ~= 0) | (B ~= 0));
            if ~any(m(:))
                continue;
            end
            cc = corrcoef(T(m), B(m));
            best = max(best, cc(1,2));
        end
    end
end

c = best;
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
