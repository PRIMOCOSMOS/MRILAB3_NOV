function outFile = remove_dummy_tr_spm(inFile, outDir, nDummy, cfg)
% remove_dummy_tr_spm - 使用 SPM 对 4D 功能像移除起始 Dummy TR

if nargin < 4
    cfg = struct();
end

spmDir = resolve_spm_dir(cfg);
if ~exist(fullfile(spmDir, 'spm.m'), 'file')
    error('[remove_dummy_tr_spm] 未找到 SPM: %s', spmDir);
end

addpath(spmDir);
spm('defaults', 'FMRI');

V = spm_vol(inFile);
nt = numel(V);
if nt <= nDummy
    error('[remove_dummy_tr_spm] 时间点数 %d <= 去除量 %d', nt, nDummy);
end

idx = (nDummy + 1):nt;
scans = cell(numel(idx), 1);
for i = 1:numel(idx)
    scans{i} = sprintf('%s,%d', inFile, idx(i));
end

ensure_dir(outDir);
[~, fname, ext] = fileparts(inFile);
outFile = fullfile(outDir, ['a' fname ext]);

TR = 0;
if isfield(cfg, 'TR') && ~isempty(cfg.TR)
    TR = cfg.TR;
end

spm_file_merge(char(scans), outFile, TR);
fprintf('[remove_dummy_tr_spm] 输出=%s\n', outFile);
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
