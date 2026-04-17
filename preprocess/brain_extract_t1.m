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

% 可选：平滑掩模边界（尽量不依赖额外工具箱）
if isfield(cfg, 'bet') && isfield(cfg.bet, 'smoothSigma') && cfg.bet.smoothSigma > 0
    try
        % 取约 ±2σ 范围的奇数长度核（2*ceil(2σ)+1，最小 3），兼顾平滑效果与计算量
        ksz = max(3, 2*ceil(2*cfg.bet.smoothSigma)+1);
        x = -(ksz-1)/2:(ksz-1)/2;
        g = exp(-(x.^2)/(2*cfg.bet.smoothSigma^2));
        g = g / sum(g);
        maskF = convn(single(mask), reshape(g,[],1,1), 'same');
        maskF = convn(maskF, reshape(g,1,[],1), 'same');
        maskF = convn(maskF, reshape(g,1,1,[]), 'same');
        mask = maskF > 0.5;
    catch
        % 忽略平滑失败，保持基础阈值掩模
    end
end

volBet = vol;
volBet(~mask) = 0;

ensure_dir(outDir);
[~, fname, ext] = fileparts(inFile);
outFile  = fullfile(outDir, ['bet_' fname ext]);
maskFile = fullfile(outDir, ['betmask_' fname ext]);

hdrVol = hdr;
hdrVol.descrip = sprintf('BET percentile=%g', p);
nifti_write(outFile, single(volBet), hdrVol);

hdrMask = hdr;
hdrMask.descrip = sprintf('BETMask percentile=%g', p);
nifti_write(maskFile, uint8(mask), hdrMask);

fprintf('[brain_extract_t1] 已写出: %s\n', outFile);
fprintf('[brain_extract_t1] 掩模: %s\n', maskFile);
end
