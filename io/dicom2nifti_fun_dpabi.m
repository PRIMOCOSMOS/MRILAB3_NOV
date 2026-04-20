function outFile = dicom2nifti_fun_dpabi(dicomDir, outDir, cfg)
% dicom2nifti_fun_dpabi - 使用 DPABI 自带 dcm2niix 转换功能像 DICOM

if nargin < 3
    cfg = struct();
end

exePath = resolve_dcm2niix_exe(cfg);
if isempty(exePath) || ~exist(exePath, 'file')
    error('[dicom2nifti_fun_dpabi] 未找到 dcm2niix 可执行文件');
end
if ~exist(dicomDir, 'dir')
    error('[dicom2nifti_fun_dpabi] DICOM 目录不存在: %s', dicomDir);
end

ensure_dir(outDir);
tmpDir = tempname;
mkdir(tmpDir);
cleanupObj = onCleanup(@() safe_rmdir(tmpDir));

opt = '-b y -x y -z n';
if isfield(cfg, 'dpabi') && isfield(cfg.dpabi, 'dcm2niixOptionFun') && ~isempty(cfg.dpabi.dcm2niixOptionFun)
    opt = cfg.dpabi.dcm2niixOptionFun;
end

cmd = sprintf('"%s" %s -f "%%f" -o "%s" "%s"', exePath, opt, tmpDir, dicomDir);
status = system(cmd);
if status ~= 0
    error('[dicom2nifti_fun_dpabi] dcm2niix 执行失败，status=%d', status);
end

niiFiles = dir(fullfile(tmpDir, '*.nii'));
if isempty(niiFiles)
    error('[dicom2nifti_fun_dpabi] 未生成 NIfTI 文件');
end

bestFile = '';
bestNt = -1;
bestSize = -1;
for i = 1:numel(niiFiles)
    p = fullfile(niiFiles(i).folder, niiFiles(i).name);
    [d, ~] = nifti_read(p);
    nt = 1;
    if ndims(d) == 4
        nt = size(d, 4);
    end
    if nt > bestNt || (nt == bestNt && niiFiles(i).bytes > bestSize)
        bestNt = nt;
        bestSize = niiFiles(i).bytes;
        bestFile = p;
    end
end

if isempty(bestFile)
    error('[dicom2nifti_fun_dpabi] 无法选择输出 NIfTI');
end

outFile = fullfile(outDir, 'bold_4d.nii');
copyfile(bestFile, outFile);

[pth, name, ~] = fileparts(bestFile);
jsonSrc = fullfile(pth, [name '.json']);
if exist(jsonSrc, 'file')
    copyfile(jsonSrc, fullfile(outDir, 'bold_4d.json'));
end

fprintf('[dicom2nifti_fun_dpabi] dcm2niix=%s\n', exePath);
fprintf('[dicom2nifti_fun_dpabi] 输出=%s\n', outFile);

clear cleanupObj;
safe_rmdir(tmpDir);
end

function exePath = resolve_dcm2niix_exe(cfg)
exePath = '';
cands = {};

if isstruct(cfg) && isfield(cfg, 'dpabi') && isfield(cfg.dpabi, 'dcm2niixExeCandidates')
    cands = cfg.dpabi.dcm2niixExeCandidates;
elseif isstruct(cfg) && isfield(cfg, 't1') && isfield(cfg.t1, 'dcm2niixExeCandidates')
    cands = cfg.t1.dcm2niixExeCandidates;
end

if ~iscell(cands)
    cands = {cands};
end

for i = 1:numel(cands)
    p = cands{i};
    if ischar(p) && ~isempty(p) && exist(p, 'file')
        exePath = p;
        return;
    end
end
end

function safe_rmdir(p)
if exist(p, 'dir')
    try
        rmdir(p, 's');
    catch
    end
end
end
