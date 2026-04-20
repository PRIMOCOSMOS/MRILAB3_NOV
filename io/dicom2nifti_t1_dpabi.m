function outFile = dicom2nifti_t1_dpabi(dicomDir, outDir, cfg)
% dicom2nifti_t1_dpabi - 使用 DPABI 自带 dcm2niix 转换 T1 DICOM

if nargin < 3
    cfg = struct();
end

exePath = resolve_dcm2niix_exe(cfg);
if isempty(exePath) || ~exist(exePath, 'file')
    error('[dicom2nifti_t1_dpabi] 未找到 dcm2niix 可执行文件');
end
if ~exist(dicomDir, 'dir')
    error('[dicom2nifti_t1_dpabi] DICOM 目录不存在: %s', dicomDir);
end

ensure_dir(outDir);
tmpDir = tempname;
mkdir(tmpDir);
cleanupObj = onCleanup(@() safe_rmdir(tmpDir));

opt = '-b y -x y -z n';
if isfield(cfg, 'dpabi') && isfield(cfg.dpabi, 'dcm2niixOptionT1') && ~isempty(cfg.dpabi.dcm2niixOptionT1)
    opt = cfg.dpabi.dcm2niixOptionT1;
end

cmd = sprintf('"%s" %s -f "%%f" -o "%s" "%s"', exePath, opt, tmpDir, dicomDir);
status = system(cmd);
if status ~= 0
    error('[dicom2nifti_t1_dpabi] dcm2niix 执行失败，status=%d', status);
end

niiFiles = dir(fullfile(tmpDir, '*.nii'));
if isempty(niiFiles)
    error('[dicom2nifti_t1_dpabi] 未生成 NIfTI 文件');
end

bestFile = '';
bestVol = -1;
bestSize = -1;
for i = 1:numel(niiFiles)
    p = fullfile(niiFiles(i).folder, niiFiles(i).name);
    [d, ~] = nifti_read(p);
    dims = size(d);
    if numel(dims) < 3
        continue;
    end
    vol = double(dims(1)) * double(dims(2)) * double(dims(3));
    if vol > bestVol || (vol == bestVol && niiFiles(i).bytes > bestSize)
        bestVol = vol;
        bestSize = niiFiles(i).bytes;
        bestFile = p;
    end
end

if isempty(bestFile)
    error('[dicom2nifti_t1_dpabi] 无法选择输出 NIfTI');
end

outFile = fullfile(outDir, 't1.nii');
copyfile(bestFile, outFile);

[pth, name, ~] = fileparts(bestFile);
jsonSrc = fullfile(pth, [name '.json']);
if exist(jsonSrc, 'file')
    copyfile(jsonSrc, fullfile(outDir, 't1.json'));
end

fprintf('[dicom2nifti_t1_dpabi] dcm2niix=%s\n', exePath);
fprintf('[dicom2nifti_t1_dpabi] 输出=%s\n', outFile);

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
