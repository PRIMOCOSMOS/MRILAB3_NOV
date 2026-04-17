function cfg = config_sub01()
% config_sub01 - 被试 Sub_01 的 pipeline 全局参数配置
% 修改此文件以适配不同被试或扫描协议
% 输出: cfg - 结构体，包含所有 pipeline 参数
%
% 使用方法: cfg = config_sub01();

% ====== 基础路径 ======
cfg.baseDir  = 'D:\MRI_PRO\MRILAB3\BOLDCODE\BOLDDATA';
cfg.subID    = 'Sub_01';

% 原始数据路径
cfg.funRawDir = fullfile(cfg.baseDir, 'FunRaw', cfg.subID);
cfg.t1RawDir  = fullfile(cfg.baseDir, 'T1Raw',  cfg.subID);

% ====== 输出目录（参考 DPABI 目录命名风格）======
cfg.funImgDir       = fullfile(cfg.baseDir, 'FunImg',          cfg.subID);
cfg.funImgADir      = fullfile(cfg.baseDir, 'FunImgA',         cfg.subID);  % 去 Dummy TR
cfg.funImgARDir     = fullfile(cfg.baseDir, 'FunImgAR',        cfg.subID);  % Realigned
cfg.funImgARWDir    = fullfile(cfg.baseDir, 'FunImgARW',       cfg.subID);  % Normalized
cfg.funImgARWSDir   = fullfile(cfg.baseDir, 'FunImgARWS',      cfg.subID);  % Smoothed
cfg.t1ImgDir        = fullfile(cfg.baseDir, 'T1Img',           cfg.subID);
cfg.t1ImgCoregDir   = fullfile(cfg.baseDir, 'T1ImgCoreg',      cfg.subID);
cfg.t1SegDir        = fullfile(cfg.baseDir, 'T1ImgNewSegment', cfg.subID);
cfg.realignParamDir = fullfile(cfg.baseDir, 'RealignParameter',cfg.subID);
cfg.firstLevelDir   = fullfile(cfg.baseDir, [cfg.subID '_1stLevel']);
cfg.logDir          = fullfile(cfg.baseDir, 'Logs',            cfg.subID);

% 所有输出目录列表，用于一键创建
cfg.outDirs = {
    cfg.funImgDir,   cfg.funImgADir,  cfg.funImgARDir,
    cfg.funImgARWDir, cfg.funImgARWSDir,
    cfg.t1ImgDir,    cfg.t1ImgCoregDir, cfg.t1SegDir,
    cfg.realignParamDir, cfg.firstLevelDir, cfg.logDir
};

% ====== EPI/MOSAIC 扫描参数 ======
cfg.TR          = 2.0;      % 重复时间 (秒)
cfg.nSlices     = 36;       % 层数
cfg.sliceSize   = 64;       % 每层体素大小 (sliceSize x sliceSize)，用于MOSAIC解码
cfg.mosaicSize  = 384;      % MOSAIC图像边长（例如 384 = 6x64 每行6层）

% 层采集时间 (ms)，相对于每个 TR 开始时刻
% 此处示例为升序隔层采集，请按实际扫描协议修改
% 方法: 连续升序顺序采集
cfg.sliceTimingMs = linspace(0, cfg.TR*1000 - cfg.TR*1000/cfg.nSlices, cfg.nSlices);
% 若为隔层采集（奇数层先，偶数层后），可改为:
% oddIdx  = 1:2:cfg.nSlices;
% evenIdx = 2:2:cfg.nSlices;
% cfg.sliceTimingMs = zeros(1, cfg.nSlices);
% cfg.sliceTimingMs(oddIdx)  = linspace(0, cfg.TR*1000/2, numel(oddIdx));
% cfg.sliceTimingMs(evenIdx) = linspace(cfg.TR*1000/2, cfg.TR*1000, numel(evenIdx));

% 参考层（切片时序校正的目标时间点，通常为中间层或第1层）
cfg.refSliceIdx = 18;  % 第18层（1-based index），对应时间 cfg.sliceTimingMs(18)

% ====== 去除起始不稳定 TR ======
cfg.nDummy = 10;  % 去掉前10个TR（根据实验协议修改）

% ====== 头动校正参数 ======
cfg.realign.quality   = 0.9;   % 降采样质量因子 (0-1)
cfg.realign.sep       = 4;     % 估计时使用的体素间距 (mm)
cfg.realign.fwhm      = 5;     % 估计前预平滑 FWHM (mm)
cfg.realign.rtm       = 1;     % 是否配准到均值像 (1=是, 0=配准到第1帧)
cfg.realign.maxIter   = 24;    % Gauss-Newton 最大迭代次数
cfg.realign.tol       = 1e-8;  % 收敛阈值

% ====== 空间平滑 ======
cfg.fwhm = [6 6 6];  % 高斯核 FWHM (mm) [x y z]

% ====== 分割参数 ======
cfg.seg.nClasses  = 3;    % 分割类别数 (GM, WM, CSF)
cfg.seg.nIter     = 100;  % GMM EM 最大迭代次数
cfg.seg.mrfBeta   = 0.1;  % MRF 正则化系数 (0 = 不使用 MRF)

% ====== 非线性配准（DARTEL 替代）======
cfg.dartel.nLevels  = 3;    % 多分辨率层数
cfg.dartel.nIter    = [3 3 6];  % 各分辨率层迭代次数
cfg.dartel.reg      = 1.0;  % 正则化权重

% ====== 标准化目标空间 ======
% MNI 152 标准空间网格参数（2mm 各向同性）
cfg.mni.origin  = [91 109 91];   % MNI 模板中心体素坐标（1-based），对应世界坐标 [0 0 0] mm
                                  % 即仿射矩阵平移分量: T = -voxSize * (origin - 1)
cfg.mni.voxSize = [2 2 2];       % 体素尺寸 (mm)
cfg.mni.dims    = [91 109 91];   % 图像维度

% ====== 一阶 GLM 参数 ======
cfg.hpf   = 128;      % 高通滤波截止周期 (秒)
cfg.units = 'secs';   % onset/duration 单位

% 示例任务条件（按实际实验设计修改）
cfg.cond.names     = {'Task', 'Rest'};
cfg.cond.onsets    = {[20 60 100 140], [0 40 80 120]};   % 每个 onset (秒)
cfg.cond.durations = {[20 20 20 20],   [20 20 20 20]};   % 持续时间 (秒)

% T-contrast: Task > Rest
cfg.tcons.name   = 'Task_gt_Rest';
cfg.tcons.weight = [1 -1 zeros(1, 6)];  % 2 条件 + 6 头动参数

end
