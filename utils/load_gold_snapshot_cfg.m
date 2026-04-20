function cfg = load_gold_snapshot_cfg(cfg)
% load_gold_snapshot_cfg - 从 Gold 目录的 DPARSFA_AutoSave*.mat 读取参数并覆盖 cfg

if nargin < 1
    cfg = struct();
end

if ~isfield(cfg, 'gold') || ~isfield(cfg.gold, 'autoLoadSnapshot') || ~logical(cfg.gold.autoLoadSnapshot)
    return;
end

goldBase = get_gold_base(cfg);
if isempty(goldBase) || ~exist(goldBase, 'dir')
    warning('[load_gold_snapshot_cfg] Gold 根目录不存在，跳过快照加载: %s', goldBase);
    return;
end

snapshotFile = '';
if isfield(cfg.gold, 'snapshotFile') && ischar(cfg.gold.snapshotFile) && exist(cfg.gold.snapshotFile, 'file')
    snapshotFile = cfg.gold.snapshotFile;
else
    cand = dir(fullfile(goldBase, 'DPARSFA_AutoSave_*.mat'));
    if isempty(cand)
        warning('[load_gold_snapshot_cfg] 未找到 DPARSFA_AutoSave_*.mat: %s', goldBase);
        return;
    end
    [~, idx] = max([cand.datenum]);
    snapshotFile = fullfile(cand(idx).folder, cand(idx).name);
end

S = load(snapshotFile);
if ~isfield(S, 'Cfg') || ~isstruct(S.Cfg)
    warning('[load_gold_snapshot_cfg] 快照不包含 Cfg 结构: %s', snapshotFile);
    return;
end

G = S.Cfg;

if isfield(G, 'TR') && ~isempty(G.TR)
    cfg.TR = double(G.TR);
end
if isfield(G, 'RemoveFirstTimePoints') && ~isempty(G.RemoveFirstTimePoints)
    cfg.nDummy = double(G.RemoveFirstTimePoints);
end
if isfield(G, 'TimePoints') && ~isempty(G.TimePoints)
    cfg.gold.timePoints = double(G.TimePoints);
end

if isfield(G, 'SliceTiming') && isstruct(G.SliceTiming)
    st = G.SliceTiming;
    if isfield(st, 'SliceNumber') && ~isempty(st.SliceNumber)
        cfg.nSlices = double(st.SliceNumber);
    end
    if isfield(st, 'SliceOrder') && ~isempty(st.SliceOrder)
        cfg.sliceTimingMs = double(st.SliceOrder(:)');
    end
    if isfield(st, 'ReferenceSlice') && ~isempty(st.ReferenceSlice) && isfield(cfg, 'sliceTimingMs') && ~isempty(cfg.sliceTimingMs)
        cfg.refSliceIdx = resolve_reference_slice_idx(cfg.sliceTimingMs, double(st.ReferenceSlice));
    end
end

if isfield(G, 'Smooth') && isstruct(G.Smooth) && isfield(G.Smooth, 'FWHM') && ~isempty(G.Smooth.FWHM)
    cfg.fwhm = double(G.Smooth.FWHM(:)');
end

if isfield(G, 'Normalize') && isstruct(G.Normalize)
    if isfield(G.Normalize, 'BoundingBox') && isequal(size(G.Normalize.BoundingBox), [2 3])
        cfg.normalize.boundingBox = double(G.Normalize.BoundingBox);
    end
    if isfield(G.Normalize, 'VoxSize') && numel(G.Normalize.VoxSize) == 3
        cfg.normalize.voxSize = double(G.Normalize.VoxSize(:)');
    end
end

if isfield(G, 'Segment') && isstruct(G.Segment) && ...
   isfield(G.Segment, 'AffineRegularisationInSegmentation') && ~isempty(G.Segment.AffineRegularisationInSegmentation)
    cfg.seg.affreg = lower(char(G.Segment.AffineRegularisationInSegmentation));
end

cfg = refresh_normalize_compat_fields(cfg);

subID = get_subid(cfg);
funMat = fullfile(goldBase, 'ReorientMats', sprintf('%s_ReorientFunImgMat.mat', subID));
t1Mat = fullfile(goldBase, 'ReorientMats', sprintf('%s_ReorientT1ImgMat.mat', subID));
if ~isfield(cfg, 'reorient') || ~isstruct(cfg.reorient)
    cfg.reorient = struct();
end
if exist(funMat, 'file')
    cfg.reorient.funMatFile = funMat;
end
if exist(t1Mat, 'file')
    cfg.reorient.t1MatFile = t1Mat;
end

cfg.gold.snapshotFileResolved = snapshotFile;
cfg.gold.snapshotLoaded = true;

fprintf('[load_gold_snapshot_cfg] 已加载 Gold 快照: %s\n', snapshotFile);
end

function goldBase = get_gold_base(cfg)
goldBase = 'D:/MRI_PRO/MRILAB3/BOLDDATA';
if isfield(cfg, 'gold') && isfield(cfg.gold, 'baseDir') && ~isempty(cfg.gold.baseDir)
    goldBase = cfg.gold.baseDir;
end
end

function idx = resolve_reference_slice_idx(sliceTimingMs, refValueMs)
[~, idx] = min(abs(double(sliceTimingMs(:)) - double(refValueMs)));
idx = max(1, min(numel(sliceTimingMs), idx));
end

function cfg = refresh_normalize_compat_fields(cfg)
if ~isfield(cfg, 'normalize') || ~isfield(cfg.normalize, 'boundingBox') || ~isfield(cfg.normalize, 'voxSize')
    return;
end
bboxMin = cfg.normalize.boundingBox(1,:);
bboxMax = cfg.normalize.boundingBox(2,:);
vox = cfg.normalize.voxSize;
cfg.normalize.dims = round((bboxMax - bboxMin) ./ vox) + 1;
cfg.normalize.origin = 1 - bboxMin ./ vox;
end

function subID = get_subid(cfg)
subID = 'Sub_01';
if isfield(cfg, 'subID') && ischar(cfg.subID) && ~isempty(cfg.subID)
    subID = cfg.subID;
end
end
