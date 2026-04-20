function tests = test_realign_regression
% test_realign_regression - 回归测试：避免 realign 发散导致全零体积
tests = functiontests(localfunctions);
end

function setupOnce(testCase)
rootDir = fileparts(fileparts(mfilename('fullpath')));
addpath(rootDir);
addpath(fullfile(rootDir, 'utils'));
addpath(fullfile(rootDir, 'preprocess'));

tmpDir = tempname;
mkdir(tmpDir);

testCase.TestData.rootDir = rootDir;
testCase.TestData.tmpDir = tmpDir;
end

function teardownOnce(testCase)
tmpDir = testCase.TestData.tmpDir;
if exist(tmpDir, 'dir')
    rmdir(tmpDir, 's');
end
end

function test_realign_no_allzero_volumes_on_real_data(testCase)
cfg = config_sub01();
inFile = fullfile(cfg.funImgADir, 'stabold_4d.nii');
assumeTrue(testCase, exist(inFile, 'file') == 2, ...
    sprintf('跳过：未找到真实数据文件 %s', inFile));

% 仅截取前20个时间点，加速回归测试
[data, hdr] = nifti_read(inFile);
ntUse = min(20, size(data, 4));
subset = single(data(:,:,:,1:ntUse));

subInDir = fullfile(testCase.TestData.tmpDir, 'in');
outDir = fullfile(testCase.TestData.tmpDir, 'out');
rpDir = fullfile(testCase.TestData.tmpDir, 'rp');
ensure_dir(subInDir);
ensure_dir(outDir);
ensure_dir(rpDir);

hdrSubset = hdr;
hdrSubset.dim = int16([4, size(subset,1), size(subset,2), size(subset,3), ntUse, 1, 1, 1]);
hdrSubset.nt = ntUse;
hdrSubset.datatype = 16;
hdrSubset.bitpix = 32;
hdrSubset.scl_slope = 1;
hdrSubset.scl_inter = 0;

subsetFile = fullfile(subInDir, 'stabold_subset.nii');
nifti_write(subsetFile, subset, hdrSubset);

cfgRun = cfg;
cfgRun.subID = 'UT_REALIGN';
cfgRun.realign.maxIter = 24;
cfgRun.realign.tol = 1e-8;
cfgRun.realign.rtm = 1;
cfgRun.realign.sep = 4;

[outFile, rpFile] = realign_estimate_reslice(subsetFile, outDir, rpDir, cfgRun);

[outData, ~] = nifti_read(outFile);
volStd = squeeze(std(reshape(outData, [], size(outData,4)), 0, 1));
volNnz = squeeze(sum(reshape(outData ~= 0, [], size(outData,4)), 1));

verifyTrue(testCase, all(isfinite(outData(:))), 'realign 输出包含非有限值');
verifyTrue(testCase, all(volStd > 1e-3), '存在近乎常数体积（疑似重采样失败）');
verifyTrue(testCase, all(volNnz > 0), '存在全零体积（疑似映射到边界外）');

P = readmatrix(rpFile);
maxTransMm = max(abs(P(:,1:3)), [], 'all');
maxRotRad = max(abs(P(:,4:6)), [], 'all');

% 保守上界：正常 realign 不应出现极端发散值
verifyLessThan(testCase, maxTransMm, 80, '平移参数异常发散');
verifyLessThan(testCase, maxRotRad, 1.2, '旋转参数异常发散');
end

function test_realign_synthetic_stability(testCase)
[nx, ny, nz, nt] = deal(48, 48, 24, 12);
[X, Y, Z] = ndgrid(1:nx, 1:ny, 1:nz);
base = exp(-((X-24).^2 + (Y-24).^2 + (Z-12).^2) / (2*7^2));
base = base + 0.1 * exp(-((X-16).^2 + (Y-34).^2 + (Z-9).^2) / (2*4^2));

rng(7);
vol4d = zeros(nx, ny, nz, nt, 'single');
for t = 1:nt
    vol4d(:,:,:,t) = single(base + 0.01 * randn(nx, ny, nz));
end

hdr = nifti_default_hdr([nx ny nz nt], [3 3 3 2]);
hdr.datatype = 16;
hdr.bitpix = 32;

inDir = fullfile(testCase.TestData.tmpDir, 'syn_in');
outDir = fullfile(testCase.TestData.tmpDir, 'syn_out');
rpDir = fullfile(testCase.TestData.tmpDir, 'syn_rp');
ensure_dir(inDir);
ensure_dir(outDir);
ensure_dir(rpDir);

inFile = fullfile(inDir, 'syn4d.nii');
nifti_write(inFile, vol4d, hdr);

cfgRun = struct();
cfgRun.subID = 'UT_SYN';
cfgRun.realign = struct('maxIter', 20, 'tol', 1e-8, 'rtm', 1, 'sep', 4);

[outFile, rpFile] = realign_estimate_reslice(inFile, outDir, rpDir, cfgRun);
[outData, ~] = nifti_read(outFile);

verifyTrue(testCase, all(isfinite(outData(:))), '合成数据输出包含非有限值');
verifyGreaterThan(testCase, std(outData(:)), 1e-3, '合成数据输出异常平坦');

P = readmatrix(rpFile);
verifyLessThan(testCase, max(abs(P(:,1:3)), [], 'all'), 10, '合成数据平移估计过大');
verifyLessThan(testCase, max(abs(P(:,4:6)), [], 'all'), 0.2, '合成数据旋转估计过大');
end
