function outFile = slice_timing_corr_spm(inFile, outDir, sliceTimingMs, refSliceIdx, TR, cfg)
% slice_timing_corr_spm - 使用 SPM Slice Timing 执行切片时序校正

if nargin < 6
    cfg = struct();
end

if numel(sliceTimingMs) < refSliceIdx || refSliceIdx < 1
    error('[slice_timing_corr_spm] refSliceIdx 超出范围: %d', refSliceIdx);
end

spmDir = resolve_spm_dir(cfg);
if ~exist(fullfile(spmDir, 'spm.m'), 'file')
    error('[slice_timing_corr_spm] 未找到 SPM: %s', spmDir);
end

addpath(spmDir);
spm('defaults', 'FMRI');
spm_jobman('initcfg');

ensure_dir(outDir);
[~, inBase, inExt] = fileparts(inFile);
workInFile = fullfile(outDir, [inBase inExt]);
if ~is_same_path(inFile, workInFile)
    copyfile(inFile, workInFile);
else
    workInFile = inFile;
end

V = spm_vol(workInFile);
nVol = numel(V);
if nVol < 2
    error('[slice_timing_corr_spm] 输入不是有效的4D时间序列: %s', workInFile);
end

scans = cell(nVol, 1);
for k = 1:nVol
    scans{k} = sprintf('%s,%d', workInFile, k);
end

sliceOrderOrTimes = double(sliceTimingMs(:)');
refSliceTimeMs = double(sliceTimingMs(refSliceIdx));
timing = [0, TR];

outFile = fullfile(outDir, ['st' inBase inExt]);
if exist(outFile, 'file')
    delete(outFile);
end

fprintf('[slice_timing_corr_spm] 使用 SPM Slice Timing: nVol=%d, refSlice=%d\n', nVol, refSliceIdx);
spm_slice_timing(char(scans), sliceOrderOrTimes, refSliceTimeMs, timing, 'st');

if ~exist(outFile, 'file')
    error('[slice_timing_corr_spm] 未生成输出文件: %s', outFile);
end

fprintf('[slice_timing_corr_spm] 已写出: %s\n', outFile);
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

function tf = is_same_path(a, b)
na = lower(strrep(char(a), '\', '/'));
nb = lower(strrep(char(b), '\', '/'));
tf = strcmp(na, nb);
end
