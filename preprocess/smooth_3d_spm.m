function outFile = smooth_3d_spm(inFile, outDir, fwhm_mm, cfg)
% smooth_3d_spm - 使用 SPM Spatial Smooth 进行空间平滑

if nargin < 4
    cfg = struct();
end

if isscalar(fwhm_mm)
    fwhm_mm = [fwhm_mm, fwhm_mm, fwhm_mm];
end

spmDir = resolve_spm_dir(cfg);
if ~exist(fullfile(spmDir, 'spm.m'), 'file')
    error('[smooth_3d_spm] 未找到 SPM: %s', spmDir);
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
if nVol < 1
    error('[smooth_3d_spm] 无法读取输入影像: %s', workInFile);
end

scans = cell(nVol, 1);
for k = 1:nVol
    scans{k} = sprintf('%s,%d', workInFile, k);
end

matlabbatch = {};
matlabbatch{1}.spm.spatial.smooth.data = scans;
matlabbatch{1}.spm.spatial.smooth.fwhm = double(fwhm_mm(:)');
matlabbatch{1}.spm.spatial.smooth.dtype = 0;
matlabbatch{1}.spm.spatial.smooth.im = 0;
matlabbatch{1}.spm.spatial.smooth.prefix = 's';

fprintf('[smooth_3d_spm] 使用 SPM Smooth: nVol=%d, FWHM=[%.1f %.1f %.1f]\n', ...
    nVol, fwhm_mm(1), fwhm_mm(2), fwhm_mm(3));
spm_jobman('run', matlabbatch);

outFile = fullfile(outDir, ['s' inBase inExt]);
if ~exist(outFile, 'file')
    error('[smooth_3d_spm] 未生成平滑输出: %s', outFile);
end

fprintf('[smooth_3d_spm] 已写出: %s\n', outFile);
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
