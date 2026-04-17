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
cfg.t1ImgBetDir     = fullfile(cfg.baseDir, 'T1ImgBet',        cfg.subID);  % BET 后的 T1
cfg.t1ImgCoregDir   = fullfile(cfg.baseDir, 'T1ImgCoreg',      cfg.subID);
cfg.t1SegDir        = fullfile(cfg.baseDir, 'T1ImgNewSegment', cfg.subID);
cfg.realignParamDir = fullfile(cfg.baseDir, 'RealignParameter',cfg.subID);
cfg.firstLevelDir   = fullfile(cfg.baseDir, [cfg.subID '_1stLevel']);
cfg.logDir          = fullfile(cfg.baseDir, 'Logs',            cfg.subID);
cfg.reorientMatDir  = fullfile(cfg.baseDir, 'ReorientMats',    cfg.subID);
cfg.qcDir           = fullfile(cfg.baseDir, 'QC',              cfg.subID);
cfg.chkNormDir      = fullfile(cfg.baseDir, 'PicturesForChkNormalization', cfg.subID);
cfg.maskDir         = fullfile(cfg.baseDir, 'Masks',           cfg.subID);

% 所有输出目录列表，用于一键创建
cfg.outDirs = { ...
    cfg.funImgDir,   cfg.funImgADir,  cfg.funImgARDir, ...
    cfg.funImgARWDir, cfg.funImgARWSDir, ...
    cfg.t1ImgDir, cfg.t1ImgBetDir, cfg.t1ImgCoregDir, cfg.t1SegDir, ...
    cfg.realignParamDir, cfg.firstLevelDir, cfg.logDir, ...
    cfg.reorientMatDir, cfg.qcDir, cfg.chkNormDir, cfg.maskDir ...
};

% ====== 应用级模板与资源路径（必须配置）======
% 说明:
% 1) 以下模板不随代码仓库分发，需使用者本地准备并填写绝对路径
% 2) 所有模板会在 run_pipeline_sub01 启动时做 fail-fast 校验
% 3) 默认命名参考 SPM/DPABI 常见模板，但实现不依赖这些工具箱
cfg.templates.dartel.gmTemplateNii = 'D:\MRI_PRO\MRILAB3\Templates\EastAsian\Template_GM.nii';
cfg.templates.dartel.wmTemplateNii = 'D:\MRI_PRO\MRILAB3\Templates\EastAsian\Template_WM.nii';
cfg.templates.standard.brainMaskNii = 'D:\MRI_PRO\MRILAB3\Templates\MNI\BrainMask_2mm.nii';
cfg.templates.standard.t1TemplateNii = 'D:\MRI_PRO\MRILAB3\Templates\MNI\MNI152_T1_2mm.nii';

% Renderer（交互式3D显示）所需模板
cfg.visualization.enable = true;                     % 是否在1st-level后自动出3D交互图
cfg.visualization.tThreshold = 3.0;                 % T阈值
cfg.visualization.alphaBrain = 0.15;                % 脑壳透明度
cfg.visualization.alphaActivation = 0.85;           % 激活层透明度
cfg.visualization.outputPng = true;                 % 是否导出静态截图
cfg.visualization.brainTemplateNii = cfg.templates.standard.t1TemplateNii;

% ====== EPI/MOSAIC 扫描参数 ======
cfg.TR          = 2.0;      % 重复时间 (秒)
cfg.nSlices     = 36;       % 层数
cfg.sliceSize   = 64;       % 每层体素大小 (sliceSize x sliceSize)，用于MOSAIC解码
cfg.mosaicSize  = 384;      % MOSAIC图像边长（例如 384 = 6x64 每行6层）

% 层采集时间 (ms)，相对于每个 TR 开始时刻
% 36层、多频带因子2（MB2）隔层交错采集的实测时序
% 每18层为一组（奇数/偶数组同步采集），组内按交错顺序依次排列
% 格式：cfg.sliceTimingMs(z) 表示第 z 层相对 TR 开始的采集延迟（ms）
cfg.sliceTimingMs = [0 1430 880 330 1760 1210 660 110 1540 990 440 1870 ...
                     1320 770 220 1650 1100 550 ...
                     0 1430 880 330 1760 1210 660 110 1540 990 440 1870 ...
                     1320 770 220 1650 1100 550];

% 参考层（切片时序校正的目标时间点，通常为中间层或第1层）
cfg.refSliceIdx = 18;  % 第18层（1-based index），对应时间 cfg.sliceTimingMs(18)

% ====== 去除起始不稳定 TR ======
cfg.nDummy = 6;  % 去掉前6个TR（本实验协议）

% ====== 头动校正参数 ======
cfg.realign.quality   = 0.9;   % 降采样质量因子 (0-1)
cfg.realign.sep       = 4;     % 估计时使用的体素间距 (mm)
cfg.realign.fwhm      = 5;     % 估计前预平滑 FWHM (mm)
cfg.realign.rtm       = 1;     % 是否配准到均值像 (1=是, 0=配准到第1帧)
cfg.realign.maxIter   = 24;    % Gauss-Newton 最大迭代次数
cfg.realign.tol       = 1e-8;  % 收敛阈值

% ====== 空间平滑 ======
cfg.fwhm = [6 6 6];  % 高斯核 FWHM (mm) [x y z]

% ====== T1 脑提取（BET）参数 ======
cfg.bet.percentile = 20;   % 强度阈值百分位（越大掩模越紧）
cfg.bet.smoothSigma = 1.0; % 掩模平滑核 sigma（体素）；0 表示不平滑

% ====== 分割参数 ======
cfg.seg.nClasses  = 3;    % 分割类别数 (GM, WM, CSF)
cfg.seg.nIter     = 100;  % GMM EM 最大迭代次数
cfg.seg.mrfBeta   = 0.1;  % MRF 正则化系数 (0 = 不使用 MRF)

% ====== 非线性配准（DARTEL 替代）======
cfg.dartel.nLevels  = 4;            % 多分辨率层数
cfg.dartel.nIter    = [3 3 6 8];    % 各分辨率层迭代次数
cfg.dartel.reg      = 1.0;          % 正则化权重
cfg.dartel.svfIntegrationSteps = 6; % scaling&squaring 步数

% ====== 标准化目标空间（可配置：Bounding box + Voxel size）======
% 参照 DPABI 参数形式：
%   Bounding box: [xmin ymin zmin; xmax ymax zmax]（单位 mm）
%   Voxel size:   [vx vy vz]（单位 mm）
cfg.normalize.boundingBox = [-90 -126 -72; 90 90 108];
cfg.normalize.voxSize     = [3 3 3];

% 兼容 normalize_apply 的网格参数（由 bbox 与 voxel size 自动推导）
bboxMin = cfg.normalize.boundingBox(1,:);
bboxMax = cfg.normalize.boundingBox(2,:);
cfg.normalize.dims = round((bboxMax - bboxMin) ./ cfg.normalize.voxSize) + 1;
cfg.normalize.origin = 1 - bboxMin ./ cfg.normalize.voxSize;  % 1-based 原点

% ====== 一阶 GLM 参数 ======
cfg.hpf   = 128;      % 高通滤波截止周期 (秒)
cfg.units = 'secs';   % onset/duration 单位

% 任务设计（本实验协议）
% 总 TR=156，TR=2s；去掉前6TR 后剩余 150TR（300s）
% Task 15TR(30s) + Rest 15TR(30s)，重复 5 次
cfg.cond.names     = {'Task', 'Rest'};
cfg.cond.onsets    = {[0 60 120 180 240], [30 90 150 210 270]};  % 秒（相对于去除6个dummy TR后的时间起点）
cfg.cond.durations = {30*ones(1,5),        30*ones(1,5)};         % 秒

% T-contrast: Task > Rest
cfg.tcons.name   = 'Task_gt_Rest';
cfg.tcons.weight = [1 -1 zeros(1, 6)];  % 示例：前2列对应条件列，后续可按协变量列补零

end
