function outFile = reorient_set_origin_spm(inFile, outDir, cfg, modality)
% reorient_set_origin_spm - 使用 SPM 仿射写回方式应用重定位矩阵

if nargin < 3
    cfg = struct();
end
if nargin < 4
    modality = '';
end

spmDir = resolve_spm_dir(cfg);
if ~exist(fullfile(spmDir, 'spm.m'), 'file')
    error('[reorient_set_origin_spm] 未找到 SPM: %s', spmDir);
end

addpath(spmDir);
spm('defaults', 'FMRI');

ensure_dir(outDir);
[~, inBase, inExt] = fileparts(inFile);
outFile = fullfile(outDir, ['reorient_' inBase inExt]);
outMatFile = fullfile(outDir, ['reorient_' inBase '.mat']);
if exist(outFile, 'file')
    delete(outFile);
end
if exist(outMatFile, 'file')
    delete(outMatFile);
end

if ~is_same_path(inFile, outFile)
    copyfile(inFile, outFile);
else
    error('[reorient_set_origin_spm] 输入与输出不能相同: %s', inFile);
end

[M, matFile] = load_reorient_matrix(cfg, modality);
V = spm_vol(outFile);
for i = 1:numel(V)
    oldM = spm_get_space(sprintf('%s,%d', outFile, i));
    spm_get_space(sprintf('%s,%d', outFile, i), M * oldM);
end

fprintf('[reorient_set_origin_spm] matFile=%s\n', matFile);
fprintf('[reorient_set_origin_spm] modality=%s 矩阵来源已应用\n', modality);
fprintf('[reorient_set_origin_spm] 输出=%s\n', outFile);
end

function [M, matFile] = load_reorient_matrix(cfg, modality)
matFile = '';
if isstruct(cfg) && isfield(cfg, 'reorient') && isstruct(cfg.reorient)
    if strcmpi(modality, 'fun') && isfield(cfg.reorient, 'funMatFile')
        matFile = cfg.reorient.funMatFile;
    elseif strcmpi(modality, 't1') && isfield(cfg.reorient, 't1MatFile')
        matFile = cfg.reorient.t1MatFile;
    end
end

if isempty(matFile) || ~exist(matFile, 'file')
    error('[reorient_set_origin_spm] 未找到重定位矩阵文件（modality=%s）', modality);
end

S = load(matFile);
if ~isfield(S, 'mat') || ~isequal(size(S.mat), [4 4])
    error('[reorient_set_origin_spm] 矩阵文件格式错误: %s', matFile);
end
M = double(S.mat);
end

function tf = is_same_path(a, b)
na = lower(strrep(char(a), '\\', '/'));
nb = lower(strrep(char(b), '\\', '/'));
tf = strcmp(na, nb);
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
