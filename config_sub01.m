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
% 3) 模板文件均来自 DPABI 或 SPM 安装目录，见下方说明
%
% 安装目录（用于自动推断模板路径）
% 用户本机常见路径：
%   DPABI: D:\DPABI_V9.0_250415
%   SPM25: D:\spm
cfg.installPaths.dpabiRoot = 'D:\DPABI_V9.0_250415';
cfg.installPaths.spmRoot   = 'D:\spm';                % 优先使用用户本机 SPM25 路径
cfg.installPaths.spmFallbackRoots = {'D:\spm25'};     % 兼容其他命名

% 结构链路后端：true=使用 SPM 默认实现（coreg/new-segment/dartel/normalize）
cfg.spm.useStructural = true;
cfg.spm.dir = cfg.installPaths.spmRoot;

% T1 结构链路：优先使用 dcm2niix 生成的 Crop_1 作为结构输入（更接近 DPABI 语义）
cfg.t1.useDcm2niixCrop = true;
cfg.t1.dcm2niixExeCandidates = { ...
    fullfile(cfg.installPaths.dpabiRoot, 'DPARSF', 'dcm2nii', 'dcm2niix.exe'), ...
    'D:\SRTP\MRIcron\Resources\dcm2niix.exe'};

% Coreg 模式：'identity' 保持几何不变（Gold-parity），'estimate' 使用 SPM 估计
cfg.coreg.mode = 'identity';
%
% ── DARTEL 模板 ──────────────────────────────────────────────────────
%   本 pipeline 为单被试模式，无法自建群组级 DARTEL 模板（DPABI 中
%   "DARTEL: Create Template" 需多被试联合）。
%   优先使用 SPM DARTEL 工具箱模板（若存在）：
%     <spm安装目录>/toolbox/DARTEL/Template_6_IXI555_MNI152.nii
%   若该文件不存在，启动时会自动回退扫描：
%     Template_6.nii 或 <spm安装目录>/tpm/TPM.nii
%   若需使用东亚人模板（提高中国受试者的分割精度），可替换为
%   由多名被试自建的 Template_6.nii，路径同样配置于此处。
%   DARTEL 模板支持两种配置方式:
%     A) 4D 单文件（推荐）: template4DNii + gmVolumeIndex + wmVolumeIndex
%        第1帧=GM，第2帧=WM（与 DARTEL 惯例一致）
%     B) 双文件: gmTemplateNii + wmTemplateNii
cfg.templates.dartel.template4DNii = fullfile(cfg.installPaths.spmRoot, 'toolbox', 'DARTEL', 'Template_6_IXI555_MNI152.nii');
cfg.templates.dartel.gmVolumeIndex = 1;   % 4D模板中 GM 所在帧（通常为第1帧）
cfg.templates.dartel.wmVolumeIndex = 2;   % 4D模板中 WM 所在帧（通常为第2帧）
%
% ── 标准脑掩模（Analyze 7.5 格式）──────────────────────────────────
%   来自 DPABI/Templates/ 目录，仅有 .hdr/.img 格式（无 .nii）
%   pipeline 的 nifti_read 已支持 Analyze 7.5 格式读取。
%   仿射矩阵自动从 .mat sidecar（如存在）或头中的 originator 字段解析。
cfg.templates.standard.brainMaskNii = fullfile(cfg.installPaths.dpabiRoot, 'Templates', 'BrainMask_05_61x73x61.hdr');
%
% ── T1 可视化模板 ───────────────────────────────────────────────────
%   来自 DPABI/Templates/ch2.nii（Colin Holmes T1 MNI 标准脑，非 ch2bet）
%   DPARSFA_run.m line 3255: Ch2Filename = fullfile(TemplatePath,'ch2.nii')
cfg.templates.standard.t1TemplateNii = fullfile(cfg.installPaths.dpabiRoot, 'Templates', 'ch2.nii');

% SPM 经典 Renderer 的 rend 模板（.mat）；用于记录参考逻辑
% 现代化 3D 渲染可不直接使用该文件，但可用于检查与兼容
cfg.visualization.spmRenderTemplateMat = fullfile(cfg.installPaths.spmRoot, 'rend', 'render_single_subj.mat');

% Renderer（交互式3D显示）所需模板
cfg.visualization.enable = true;                     % 是否在1st-level后自动出3D交互图
cfg.visualization.tThreshold = 3.0;                 % T阈值
cfg.visualization.alphaBrain = 0.15;                % 脑壳透明度
cfg.visualization.alphaActivation = 0.85;           % 激活层透明度
cfg.visualization.outputPng = true;                 % 是否导出静态截图
cfg.visualization.brainTemplateNii = cfg.templates.standard.t1TemplateNii;

% 参考流程一致性审计（对照 SPM25/DPABI 逻辑级流程）
cfg.referenceAudit.strict = false;  % true 时若关键流程缺失将 fail-fast

% ====== Pipeline 执行策略 ======
% true: 每次运行都从头重算，不因已有输出文件而跳过步骤
cfg.pipeline.forceRerun = true;

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
cfg.refSliceIdx = 1;  % 第18层（1-based index），对应时间 cfg.sliceTimingMs(18)

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
cfg.bet.method = 'dpabi_bet'; % 'dpabi_bet' / 'spm_prob' / 'percentile'
cfg.bet.percentile = 1;        % percentile 模式阈值（越大掩模越紧）
cfg.bet.smoothSigma = 1.0;   % percentile 模式掩模平滑 sigma（体素）
cfg.bet.dpabiOption = '';      % dpabi_bet 模式传给 y_Call_bet 的选项（结构像默认空）

% spm_prob 模式：利用 SPM 组织概率构建软脑图
cfg.bet.spmProbGamma = 0.2;
cfg.bet.spmProbThreshold = 0.01;
cfg.bet.spmProbScale = 200;

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
cfg.units = 'scans';  % onset/duration 单位: 'scans' 或 'secs'

% SPM风格一阶建模细节（用于提升与 SPM/DPABI 结果一致性）
cfg.glm.microtimeResolution   = 16;    % 对应 SPM.xBF.T
cfg.glm.microtimeOnsetBin     = 8;     % 对应 SPM.xBF.T0
cfg.glm.applyHighPassFilter   = true;  % 将高通滤波作用于 X 与 Y（SPM风格）
cfg.glm.useAR1Whitening       = true;  % 串行相关修正（SPM一阶默认近似）
cfg.glm.maskMethod            = 'globalFraction'; % 'percentile' | 'globalFraction'
cfg.glm.maskGlobalFraction    = 0.50;  % meanVol > 0.5 * globalMean
cfg.glm.explicitDriftRegressors = false; % false 时不在 X 中显式加入 DCT 漂移列

% 任务设计（本实验协议）
% 总 TR=156，TR=2s；去掉前6TR 后剩余 150TR（300s）
% 单任务 block 设计：Task 15TR + Rest 15TR，重复 5 次
% 与 SPM/DPABI 一致，Rest 作为隐式基线（不显式建模为第二个条件列）
cfg.cond.names     = {'Task'};
cfg.cond.onsets    = {[0 30 60 90 120]};  % 扫描点（去除6个dummy TR后）
cfg.cond.durations = {15*ones(1,5)};      % 扫描点

% T-contrast: Task > Baseline
% 仅需给出任务列权重，后续 nuisance 列会在统计函数中自动补零
cfg.tcons.name   = 'Task_gt_Baseline';
cfg.tcons.weight = [1];

end
