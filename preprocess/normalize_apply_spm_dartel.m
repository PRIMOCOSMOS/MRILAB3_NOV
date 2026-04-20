function outFile = normalize_apply_spm_dartel(inFile, flowFile, templateFile, outDir, mniCfg, cfg)
% normalize_apply_spm_dartel - 使用 SPM DARTEL 流场将功能像标准化到 MNI
%
% 输入:
%   inFile       - 输入 4D 功能像（个体空间）
%   flowFile     - SPM DARTEL 流场（u_rc1*.nii）
%   templateFile - DARTEL 最终模板（通常 Template_6.nii）
%   outDir       - 输出目录
%   mniCfg       - 目标网格配置，需含 boundingBox 与 voxSize
%   cfg          - 全局配置（用于解析 SPM 路径）
%
% 输出:
%   outFile - 标准化输出（w*.nii）

if nargin < 6
    cfg = struct();
end

spmDir = resolve_spm_dir(cfg);
if ~exist(fullfile(spmDir, 'spm.m'), 'file')
    error('[normalize_apply_spm_dartel] 未找到 SPM: %s', spmDir);
end
addpath(spmDir);
spm('defaults', 'FMRI');
spm_jobman('initcfg');

if nargin < 3 || isempty(templateFile) || ~exist(templateFile, 'file')
    [flowDir, ~, ~] = fileparts(flowFile);
    templateFile = fullfile(flowDir, 'Template_6.nii');
end
if ~exist(templateFile, 'file')
    error('[normalize_apply_spm_dartel] 未找到 DARTEL 模板: %s', templateFile);
end

if ~isfield(mniCfg, 'boundingBox') || ~isequal(size(mniCfg.boundingBox), [2 3])
    error('[normalize_apply_spm_dartel] mniCfg.boundingBox 必须为 2x3');
end
if ~isfield(mniCfg, 'voxSize') || numel(mniCfg.voxSize) ~= 3
    error('[normalize_apply_spm_dartel] mniCfg.voxSize 必须为 1x3');
end

ensure_dir(outDir);

job = struct();
job.template = {templateFile};
job.data.subj.flowfield = {flowFile};
job.data.subj.images = {inFile};
job.vox = mniCfg.voxSize(:)';
job.bb = mniCfg.boundingBox;
job.preserve = 0;
job.fwhm = [0 0 0];
job.output.option = 'allin';
job.output.outDir = outDir;

spm_dartel_norm_fun(job);

[~, inBase, inExt] = fileparts(inFile);
outFile = fullfile(outDir, ['w' inBase inExt]);
if ~exist(outFile, 'file')
    error('[normalize_apply_spm_dartel] 未生成标准化输出: %s', outFile);
end

fprintf('[normalize_apply_spm_dartel] 已写出: %s\n', outFile);
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
