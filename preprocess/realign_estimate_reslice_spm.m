function [outFile, rpFile] = realign_estimate_reslice_spm(inFile, outDir, rpDir, cfg)
% realign_estimate_reslice_spm - 使用 SPM Realign(Estimate+Reslice) 头动校正

if nargin < 4
    cfg = struct();
end

spmDir = resolve_spm_dir(cfg);
if ~exist(fullfile(spmDir, 'spm.m'), 'file')
    error('[realign_estimate_reslice_spm] 未找到 SPM: %s', spmDir);
end

addpath(spmDir);
spm('defaults', 'FMRI');
spm_jobman('initcfg');

ensure_dir(outDir);
ensure_dir(rpDir);

[~, inBase, inExt] = fileparts(inFile);
workInFile = fullfile(outDir, [inBase inExt]);
if ~is_same_path(inFile, workInFile)
    copyfile(inFile, workInFile);
else
    workInFile = inFile;
end

V = spm_vol(workInFile);
if numel(V) < 2
    error('[realign_estimate_reslice_spm] 输入不是有效的4D时间序列: %s', workInFile);
end

eflags = struct();
if isstruct(cfg) && isfield(cfg, 'realign')
    if isfield(cfg.realign, 'quality'), eflags.quality = cfg.realign.quality; end
    if isfield(cfg.realign, 'sep'),     eflags.sep = cfg.realign.sep; end
    if isfield(cfg.realign, 'fwhm'),    eflags.fwhm = cfg.realign.fwhm; end
    if isfield(cfg.realign, 'rtm'),     eflags.rtm = cfg.realign.rtm; end
end
if ~isfield(eflags, 'interp'), eflags.interp = 2; end
if ~isfield(eflags, 'wrap'),   eflags.wrap = [0 0 0]; end
if ~isfield(eflags, 'PW'),     eflags.PW = ''; end

fprintf('[realign_estimate_reslice_spm] 使用 SPM Realign Estimate\n');
spm_realign(V, eflags);

rflags = struct();
rflags.interp = 4;
rflags.wrap = [0 0 0];
rflags.mask = 1;
rflags.which = [2 1];
rflags.prefix = 'r';

fprintf('[realign_estimate_reslice_spm] 使用 SPM Reslice\n');
spm_reslice(V, rflags);

outFile = fullfile(outDir, ['r' inBase inExt]);
if ~exist(outFile, 'file')
    error('[realign_estimate_reslice_spm] 未生成重采样输出: %s', outFile);
end

rpSrc = find_latest_rp_file(outDir, inBase);
subID = get_subid(cfg);
rpFile = fullfile(rpDir, sprintf('rp_%s.txt', subID));
if strcmpi(rpSrc, rpFile)
    % already in place
else
    copyfile(rpSrc, rpFile);
end

fprintf('[realign_estimate_reslice_spm] 重采样输出: %s\n', outFile);
fprintf('[realign_estimate_reslice_spm] 头动参数: %s\n', rpFile);
end

function rpSrc = find_latest_rp_file(searchDir, inBase)
patternPrefer = fullfile(searchDir, sprintf('rp_%s*.txt', inBase));
cand = dir(patternPrefer);
if isempty(cand)
    cand = dir(fullfile(searchDir, 'rp_*.txt'));
end
if isempty(cand)
    error('[realign_estimate_reslice_spm] 未找到 SPM 生成的 rp_*.txt: %s', searchDir);
end
[~, idx] = max([cand.datenum]);
rpSrc = fullfile(cand(idx).folder, cand(idx).name);
end

function subID = get_subid(cfg)
subID = 'Sub_01';
if isstruct(cfg) && isfield(cfg, 'subID') && ~isempty(cfg.subID)
    subID = cfg.subID;
end
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
