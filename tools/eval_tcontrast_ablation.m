function results = eval_tcontrast_ablation()
% eval_tcontrast_ablation - Quantify T-map error attribution by ARWS/RP cross-combinations.

rootDir = fileparts(fileparts(mfilename('fullpath')));
addpath(rootDir);
addpath(fullfile(rootDir, 'utils'));
addpath(fullfile(rootDir, 'stats'));

cfg = config_sub01();
if isfield(cfg, 'visualization')
    cfg.visualization.enable = false;
end

goldTFile = 'D:/MRI_PRO/MRILAB3/BOLDDATA/Sub01_1stLevel/spmT_0001.nii';
if ~exist(goldTFile, 'file')
    error('Gold T-map not found: %s', goldTFile);
end

oursARWS = 'D:/MRI_PRO/MRILAB3/BOLDCODE/BOLDDATA/FunImgARWS/Sub_01/swreorient_rstabold_4d.nii';
goldARWS = 'D:/MRI_PRO/MRILAB3/BOLDDATA/FunImgARWS/Sub_01/swraSub_01_BOLD_lefthand_20190710154350_10.nii';
oursRP   = 'D:/MRI_PRO/MRILAB3/BOLDCODE/BOLDDATA/RealignParameter/Sub_01/rp_Sub_01.txt';
goldRP   = 'D:/MRI_PRO/MRILAB3/BOLDDATA/RealignParameter/Sub_01/rp_aSub_01_BOLD_lefthand_20190710154350_10.txt';

for p = {oursARWS, goldARWS, oursRP, goldRP}
    if ~exist(p{1}, 'file')
        error('Required file missing: %s', p{1});
    end
end

combos = {
    'oursARWS_oursRP', oursARWS, oursRP;
    'oursARWS_goldRP', oursARWS, goldRP;
    'goldARWS_oursRP', goldARWS, oursRP;
    'goldARWS_goldRP', goldARWS, goldRP
};

results = struct('name', {}, 'corrResampled', {}, 'corrVoxel', {}, 'tFile', {});

baseOut = 'D:/MRI_PRO/MRILAB3/BOLDCODE/BOLDDATA/Sub_01_1stLevel_ablation';
ensure_dir(baseOut);

fprintf('=== T-contrast Ablation ===\n');
for i = 1:size(combos,1)
    name = combos{i,1};
    funFile = combos{i,2};
    rpFile = combos{i,3};
    outDir = fullfile(baseOut, name);

    if exist(outDir, 'dir')
        try
            rmdir(outDir, 's');
        catch
        end
    end
    ensure_dir(outDir);

    run_firstlevel_glm(funFile, rpFile, outDir, cfg);

    tFile = fullfile(outDir, 'spmT_Task_gt_Baseline.nii');
    if ~exist(tFile, 'file')
        error('T-map missing for combo %s: %s', name, tFile);
    end

    corrResampled = corr_resampled(tFile, goldTFile);
    corrVoxel = corr_voxel_if_same_grid(tFile, goldTFile);

    fprintf('%-18s resampled=%.6f voxel=%.6f\n', name, corrResampled, corrVoxel);

    results(i).name = name; %#ok<AGROW>
    results(i).corrResampled = corrResampled; %#ok<AGROW>
    results(i).corrVoxel = corrVoxel; %#ok<AGROW>
    results(i).tFile = tFile; %#ok<AGROW>
end

save(fullfile(baseOut, 'tcontrast_ablation_results.mat'), 'results');
end

function c = corr_resampled(oursFile, goldFile)
[A, Ah] = nifti_read(oursFile);
[B, Bh] = nifti_read(goldFile);
A = double(A(:,:,:,1));
B = double(B(:,:,:,1));
Ar = resample_vol_affine(A, Ah.affine, Bh.affine, size(B));
m = isfinite(Ar) & isfinite(B) & ((Ar ~= 0) | (B ~= 0));
if ~any(m(:))
    c = NaN;
    return;
end
cc = corrcoef(Ar(m), B(m));
c = cc(1,2);
end

function c = corr_voxel_if_same_grid(oursFile, goldFile)
[A, ~] = nifti_read(oursFile);
[B, ~] = nifti_read(goldFile);
A = double(A(:,:,:,1));
B = double(B(:,:,:,1));
if any(size(A) ~= size(B))
    c = NaN;
    return;
end
m = isfinite(A) & isfinite(B) & ((A ~= 0) | (B ~= 0));
if ~any(m(:))
    c = NaN;
    return;
end
cc = corrcoef(A(m), B(m));
c = cc(1,2);
end
