function tests = test_compute_tcontrast_regression
% test_compute_tcontrast_regression - 回归测试：短对比向量应可自动补零并稳定计算
tests = functiontests(localfunctions);
end

function setupOnce(testCase)
rootDir = fileparts(fileparts(mfilename('fullpath')));
addpath(rootDir);
addpath(fullfile(rootDir, 'utils'));
addpath(fullfile(rootDir, 'stats'));

tmpDir = tempname;
mkdir(tmpDir);

testCase.TestData.tmpDir = tmpDir;
end

function teardownOnce(testCase)
tmpDir = testCase.TestData.tmpDir;
if exist(tmpDir, 'dir')
    rmdir(tmpDir, 's');
end
end

function test_short_contrast_vector_autopad(testCase)
rng(42);

nScans = 40;
nCols = 5;
nx = 2; ny = 2; nz = 1;
nVox = nx * ny * nz;

X = randn(nScans, nCols);
X(:,end) = 1;  % 常数项
beta = randn(nCols, nVox);
sigma2 = abs(randn(1, nVox)) + 0.1;

hdr = nifti_default_hdr([nx ny nz 1], [2 2 2 2]);
hdr.nt = 1;
hdr.dim = int16([3 nx ny nz 1 1 1 1]);

contrastShort = 1;  % 仅给任务列，剩余列应自动补零

[tMap, pMap, betaConMap, outFiles] = compute_tcontrast(...
    beta, sigma2, X, contrastShort, hdr, testCase.TestData.tmpDir, 'Task_gt_Baseline');

verifySize(testCase, tMap, [nx ny nz]);
verifySize(testCase, pMap, [nx ny nz]);
verifySize(testCase, betaConMap, [nx ny nz]);
verifyTrue(testCase, all(isfinite(tMap(:))));
verifyTrue(testCase, all(isfinite(pMap(:))));
verifyTrue(testCase, exist(outFiles.tFile, 'file') == 2);
verifyTrue(testCase, exist(outFiles.conFile, 'file') == 2);
verifyTrue(testCase, exist(outFiles.pFile, 'file') == 2);
end
