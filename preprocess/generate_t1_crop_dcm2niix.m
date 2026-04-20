function cropFile = generate_t1_crop_dcm2niix(t1RawDir, outDir, cfg)
% generate_t1_crop_dcm2niix - 使用 dcm2niix 生成 T1 Crop_1 文件
%
% 输出:
%   cropFile - 标准化命名的裁剪文件路径（reorient_t1_Crop_1.nii）

if nargin < 3
    cfg = struct();
end

exePath = resolve_dcm2niix_exe(cfg);
if isempty(exePath) || ~exist(exePath, 'file')
    error('[generate_t1_crop_dcm2niix] 未找到 dcm2niix 可执行文件');
end

if ~exist(t1RawDir, 'dir')
    error('[generate_t1_crop_dcm2niix] T1 原始目录不存在: %s', t1RawDir);
end

ensure_dir(outDir);

cmd = sprintf('"%s" -b y -x y -z n -f "%%f" -o "%s" "%s"', exePath, outDir, t1RawDir);
status = system(cmd);
if status ~= 0
    error('[generate_t1_crop_dcm2niix] dcm2niix 执行失败，status=%d', status);
end

cropCandidates = dir(fullfile(outDir, '*_Crop_1.nii'));
if isempty(cropCandidates)
    error('[generate_t1_crop_dcm2niix] 未生成 *_Crop_1.nii: %s', outDir);
end

[~, idx] = max([cropCandidates.datenum]);
cropSrc = fullfile(cropCandidates(idx).folder, cropCandidates(idx).name);

cropFile = fullfile(outDir, 'reorient_t1_Crop_1.nii');
if ~strcmpi(cropSrc, cropFile)
    copyfile(cropSrc, cropFile);
end

fprintf('[generate_t1_crop_dcm2niix] dcm2niix=%s\n', exePath);
fprintf('[generate_t1_crop_dcm2niix] crop=%s\n', cropFile);
end

function exePath = resolve_dcm2niix_exe(cfg)
exePath = '';

cands = {};
if isstruct(cfg) && isfield(cfg, 't1') && isfield(cfg.t1, 'dcm2niixExeCandidates')
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
