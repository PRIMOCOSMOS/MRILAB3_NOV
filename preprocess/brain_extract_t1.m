function [outFile, maskFile] = brain_extract_t1(inFile, outDir, cfg)
% brain_extract_t1 - 简易 T1 脑提取（BET）
%
% 方法:
%   1) 使用强度百分位阈值去除背景
%   2) 可选平滑掩模边界（若配置提供 smoothSigma）
%   3) 输出 bet_*.nii（脑提取后 T1）与 betmask_*.nii（二值掩模）

fprintf('[brain_extract_t1] 读取: %s\n', inFile);
[vol, hdr] = nifti_read(inFile);
vol = single(vol(:,:,:,1));

if use_dpabi_bet_method(cfg)
    [volBet, mask] = build_bet_dpabi(inFile, cfg);
    methodDesc = sprintf('DPABI-BET option="%s"', get_bet_field(cfg, 'dpabiOption', ''));
elseif use_spm_prob_method(cfg)
    [volBet, mask] = build_bet_spm_prob(inFile, cfg);
    methodDesc = sprintf('SPMProb gamma=%.3f th=%.3f scale=%.1f', ...
        get_bet_field(cfg, 'spmProbGamma', 0.2), ...
        get_bet_field(cfg, 'spmProbThreshold', 0.01), ...
        get_bet_field(cfg, 'spmProbScale', 200));
else
    nonZeroVals = vol(vol > 0);
    if isempty(nonZeroVals)
        error('[brain_extract_t1] 输入图像全零，无法进行 BET');
    end

    if isfield(cfg, 'bet') && isfield(cfg.bet, 'percentile')
        p = cfg.bet.percentile;
    else
        p = 20;
    end

    thr = prctile(double(nonZeroVals), p);
    mask = vol > thr;

    if isfield(cfg, 'bet') && isfield(cfg.bet, 'smoothSigma') && cfg.bet.smoothSigma > 0
        try
            ksz = max(3, 2*ceil(2*cfg.bet.smoothSigma)+1);
            x = -(ksz-1)/2:(ksz-1)/2;
            g = exp(-(x.^2)/(2*cfg.bet.smoothSigma^2));
            g = g / sum(g);
            maskF = convn(single(mask), reshape(g,[],1,1), 'same');
            maskF = convn(maskF, reshape(g,1,[],1), 'same');
            maskF = convn(maskF, reshape(g,1,1,[]), 'same');
            mask = maskF > 0.5;
        catch
        end
    end

    volBet = vol;
    volBet(~mask) = 0;
    methodDesc = sprintf('Percentile p=%g', p);
end

ensure_dir(outDir);
[~, fname, ext] = fileparts(inFile);
outFile  = fullfile(outDir, ['bet_' fname ext]);
maskFile = fullfile(outDir, ['betmask_' fname ext]);

hdrVol = hdr;
hdrVol.descrip = ['BET ' methodDesc];
nifti_write(outFile, single(volBet), hdrVol);

hdrMask = hdr;
hdrMask.descrip = ['BETMask ' methodDesc];
nifti_write(maskFile, uint8(mask), hdrMask);

fprintf('[brain_extract_t1] 已写出: %s\n', outFile);
fprintf('[brain_extract_t1] 掩模: %s\n', maskFile);
end

function tf = use_dpabi_bet_method(cfg)
tf = isfield(cfg, 'bet') && isfield(cfg.bet, 'method') && strcmpi(cfg.bet.method, 'dpabi_bet');
end

function tf = use_spm_prob_method(cfg)
tf = isfield(cfg, 'bet') && isfield(cfg.bet, 'method') && strcmpi(cfg.bet.method, 'spm_prob');
end

function val = get_bet_field(cfg, fieldName, defaultVal)
val = defaultVal;
if isfield(cfg, 'bet') && isfield(cfg.bet, fieldName) && ~isempty(cfg.bet.(fieldName))
    val = cfg.bet.(fieldName);
end
end

function [volBet, mask] = build_bet_dpabi(inFile, cfg)
dpabiRoot = resolve_dpabi_root_for_bet(cfg);
if ~exist(fullfile(dpabiRoot, 'DPABI.m'), 'file')
    error('[brain_extract_t1] DPABI 路径无效: %s', dpabiRoot);
end

addpath(dpabiRoot);
addpath(fullfile(dpabiRoot, 'DPARSF'));
addpath(genpath(fullfile(dpabiRoot, 'DPARSF', 'Subfunctions')));

if exist('y_Call_bet', 'file') ~= 2
    error('[brain_extract_t1] 未找到 y_Call_bet，请检查 DPABI 安装路径');
end

tmpDir = tempname;
mkdir(tmpDir);
oldDir = pwd;
cleanupObj = onCleanup(@() restore_dir_and_cleanup(oldDir, tmpDir));

outFile = fullfile(tmpDir, 'bet_dpabi_out.nii');
option = get_bet_field(cfg, 'dpabiOption', '');
if isfield(cfg, 'baseDir') && exist(cfg.baseDir, 'dir')
    workDir = cfg.baseDir;
else
    workDir = tmpDir;
end

y_Call_bet(inFile, outFile, option, workDir);

if ~exist(outFile, 'file')
    error('[brain_extract_t1] y_Call_bet 未生成输出: %s', outFile);
end

[betVol, ~] = nifti_read(outFile);
volBet = single(betVol(:,:,:,1));
mask = volBet > 0;

clear cleanupObj;
restore_dir_and_cleanup(oldDir, tmpDir);
end

function [volBet, mask] = build_bet_spm_prob(inFile, cfg)
spmDir = resolve_spm_dir_for_bet(cfg);
if ~exist(fullfile(spmDir, 'spm.m'), 'file')
    error('[brain_extract_t1] SPM 路径无效: %s', spmDir);
end

addpath(spmDir);
spm('defaults', 'FMRI');
spm_jobman('initcfg');

tmpDir = tempname;
mkdir(tmpDir);
cleanupObj = onCleanup(@() cleanup_tmp_dir(tmpDir));

segInput = fullfile(tmpDir, 't1.nii');
copyfile(inFile, segInput);

tpmFile = fullfile(spmDir, 'tpm', 'TPM.nii');
if ~exist(tpmFile, 'file')
    error('[brain_extract_t1] 未找到 TPM 模板: %s', tpmFile);
end

job = struct();
job.channel.vols    = {segInput};
job.channel.biasreg = 0.001;
job.channel.biasfwhm= 60;
job.channel.write   = [0 0];

ngaus = [1 1 2 3 4 2];
for k = 1:6
    job.tissue(k).tpm   = {sprintf('%s,%d', tpmFile, k)};
    job.tissue(k).ngaus = ngaus(k);
    if k <= 3
        job.tissue(k).native = [1 0];
    else
        job.tissue(k).native = [0 0];
    end
    job.tissue(k).warped = [0 0];
end

job.warp.mrf     = 1;
job.warp.cleanup = 1;
job.warp.reg     = [0 0.001 0.5 0.05 0.2];
job.warp.affreg  = 'mni';
job.warp.fwhm    = 0;
job.warp.samp    = 3;
job.warp.write   = [0 0];
job.warp.bb      = NaN(2,3);
job.warp.vox     = NaN(1,3);
job.niterations  = 1;
job.alpha        = 12;

spm_preproc_run(job, 'run');

c1File = fullfile(tmpDir, 'c1t1.nii');
c2File = fullfile(tmpDir, 'c2t1.nii');
c3File = fullfile(tmpDir, 'c3t1.nii');
if ~exist(c1File, 'file') || ~exist(c2File, 'file') || ~exist(c3File, 'file')
    error('[brain_extract_t1] SPM 概率图生成失败');
end

[p1, ~] = nifti_read(c1File);
[p2, ~] = nifti_read(c2File);
[p3, ~] = nifti_read(c3File);

p1 = double(p1(:,:,:,1));
p2 = double(p2(:,:,:,1));
p3 = double(p3(:,:,:,1));

prob = max(cat(4, p1, p2, p3), [], 4);
prob = max(prob, 0) .^ get_bet_field(cfg, 'spmProbGamma', 0.2);
thr  = get_bet_field(cfg, 'spmProbThreshold', 0.01);
prob(prob < thr) = 0;

scale = get_bet_field(cfg, 'spmProbScale', 200);
volBet = single(scale * prob);
mask = prob > 0;

clear cleanupObj;
cleanup_tmp_dir(tmpDir);
end

function dpabiRoot = resolve_dpabi_root_for_bet(cfg)
if isfield(cfg, 'installPaths') && isfield(cfg.installPaths, 'dpabiRoot') && ~isempty(cfg.installPaths.dpabiRoot)
    dpabiRoot = cfg.installPaths.dpabiRoot;
else
    dpabiRoot = 'D:/DPABI_V9.0_250415';
end
end

function spmDir = resolve_spm_dir_for_bet(cfg)
if isfield(cfg, 'spm') && isfield(cfg.spm, 'dir') && ~isempty(cfg.spm.dir)
    spmDir = cfg.spm.dir;
elseif isfield(cfg, 'installPaths') && isfield(cfg.installPaths, 'spmRoot')
    spmDir = cfg.installPaths.spmRoot;
else
    spmDir = 'D:/spm';
end
end

function cleanup_tmp_dir(tmpDir)
if exist(tmpDir, 'dir')
    try
        rmdir(tmpDir, 's');
    catch
    end
end
end

function restore_dir_and_cleanup(oldDir, tmpDir)
if exist(oldDir, 'dir')
    try
        cd(oldDir);
    catch
    end
end
cleanup_tmp_dir(tmpDir);
end
