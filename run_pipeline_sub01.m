function run_pipeline_sub01()
% run_pipeline_sub01 - BOLD-fMRI 分析 Pipeline 主入口（一键运行）
%
% ======================================================================
% 功能说明:
%   完全 Standalone 的 fMRI 分析流程，不依赖 SPM/DPABI，
%   使用 MATLAB 内置函数实现全部算法。
%   参考 SPM12 / DPABI 的处理逻辑与目录风格。
%
% 处理步骤:
%   Step 01: DICOM → NIfTI（EPI MOSAIC + T1）
%   Step 02: 去除起始不稳定 TR（Dummy Scans）
%   Step 03: 切片时序校正（Slice Timing Correction）
%   Step 04: 头动校正（Realign: Estimate + Reslice）
%   Step 05: 重定位（Reorient: 坐标原点平移至 AC）
%   Step 06: T1 脑提取（BET）
%   Step 07: T1 配准到功能像空间（Coreg T1 to Fun）
%   Step 08: T1 组织分割（GM/WM/CSF 概率图）
%   Step 09: 非线性配准（DARTEL 替代实现: SVF 速度场）
%   Step 10: 空间标准化（Normalize: 应用形变场到 MNI）
%   Step 11: 空间平滑（Smooth: 3D Gaussian）
%   Step 12: 一阶 GLM（设计矩阵 + OLS 估计 + T-contrast）
%
% 目录风格（参考 DPABI）:
%   FunImg      - 原始 NIfTI 功能像
%   FunImgA     - 去 Dummy TR 后
%   FunImgAR    - 切片时序校正 + 头动校正后（a=after, R=realigned）
%   FunImgARW   - 空间标准化后（W=warped to MNI）
%   FunImgARWS  - 平滑后（S=smoothed）
%   T1Img            - 原始 T1
%   T1ImgBet         - BET 后的 T1（bet_*.nii）
%   T1ImgCoreg       - 配准到 Fun 空间后的 T1
%   T1ImgNewSegment  - 组织分割结果（GM/WM/CSF）
%   RealignParameter - 头动参数 rp_*.txt
%   Sub01_1stLevel   - 一阶 GLM 结果
%
% 运行前须知:
%   1. 修改 config_sub01.m 中的路径和参数
%   2. 确保原始 DICOM 数据已放置于对应目录
%   3. 首次运行前无需手动创建输出目录
%
% 使用方法:
%   run_pipeline_sub01;
%
% 版本: 1.0  日期: 2026-04
% ======================================================================

clc;
fprintf('=============================================================\n');
fprintf('  BOLD-fMRI Standalone Analysis Pipeline  v1.0\n');
fprintf('  被试: Sub_01\n');
fprintf('  时间: %s\n', datestr(now));
fprintf('=============================================================\n\n');

% -------- 添加模块路径 --------
pipelineDir = fileparts(mfilename('fullpath'));
addpath(fullfile(pipelineDir, 'utils'));
addpath(fullfile(pipelineDir, 'io'));
addpath(fullfile(pipelineDir, 'preprocess'));
% register 模块函数已并入 preprocess（例如 coreg_t1_to_fun）
addpath(fullfile(pipelineDir, 'stats'));
addpath(fullfile(pipelineDir, 'visualize'));

% -------- 载入配置 --------
cfg = config_sub01();

% -------- 按 DPABI/SPM 安装目录自动推断模板路径 --------
cfg = resolve_pipeline_template_paths(cfg);

% -------- 开箱即用配置校验（失败即退出）--------
validate_pipeline_config(cfg);

% -------- 创建输出目录 --------
fprintf('[Pipeline] 创建输出目录...\n');
for i = 1:numel(cfg.outDirs)
    ensure_dir(cfg.outDirs{i});
end

% -------- 初始化日志 --------
logFile = fullfile(cfg.logDir, sprintf('pipeline_%s.log', datestr(now,'yyyymmdd_HHMMSS')));
write_log(logFile, '=== Pipeline 开始 ===');
write_log(logFile, sprintf('被试: %s', cfg.subID));

% -------- 逻辑一致性审计（对照 SPM25/DPABI 关键流程）--------
auditFile = fullfile(cfg.logDir, sprintf('parity_audit_%s.txt', datestr(now,'yyyymmdd_HHMMSS')));
audit = audit_pipeline_parity(cfg, pipelineDir, auditFile);
write_log(logFile, sprintf('ParityAudit: %s (pass=%d)', audit.reportFile, audit.pass));

% ======================================================================
% Step 01: DICOM → NIfTI
% ======================================================================
write_log(logFile, 'Step01: DICOM → NIfTI');
fprintf('\n--- Step 01: DICOM → NIfTI ---\n');

funNii = fullfile(cfg.funImgDir, 'bold_4d.nii');
t1Nii  = fullfile(cfg.t1ImgDir,  't1.nii');

if ~exist(funNii, 'file')
    funNii = dicom2nifti_fun(cfg.funRawDir, cfg.funImgDir, cfg);
else
    write_log(logFile, '  功能像 NIfTI 已存在，跳过');
end

if ~exist(t1Nii, 'file')
    t1Nii = dicom2nifti_t1(cfg.t1RawDir, cfg.t1ImgDir);
else
    write_log(logFile, '  T1 NIfTI 已存在，跳过');
end
write_log(logFile, sprintf('  EPI: %s', funNii));
write_log(logFile, sprintf('  T1:  %s', t1Nii));

% ======================================================================
% Step 02: 去除 Dummy TR
% ======================================================================
write_log(logFile, sprintf('Step02: 去除前 %d TR', cfg.nDummy));
fprintf('\n--- Step 02: 去除 Dummy TR ---\n');

funNii_A = fullfile(cfg.funImgADir, ['a' 'bold_4d.nii']);
if ~exist(funNii_A, 'file')
    funNii_A = remove_dummy_tr(funNii, cfg.funImgADir, cfg.nDummy);
else
    write_log(logFile, '  已存在，跳过');
end
write_log(logFile, sprintf('  输出: %s', funNii_A));

% ======================================================================
% Step 03: 切片时序校正
% ======================================================================
write_log(logFile, 'Step03: Slice Timing Correction');
fprintf('\n--- Step 03: Slice Timing Correction ---\n');

[~, funBase_A] = fileparts(funNii_A);
funNii_ST = fullfile(cfg.funImgADir, ['st' funBase_A '.nii']);
if ~exist(funNii_ST, 'file')
    funNii_ST = slice_timing_corr(funNii_A, cfg.funImgADir, ...
        cfg.sliceTimingMs, cfg.refSliceIdx, cfg.TR);
else
    write_log(logFile, '  已存在，跳过');
end
write_log(logFile, sprintf('  输出: %s', funNii_ST));

% ======================================================================
% Step 04: 头动校正（Realign）
% ======================================================================
write_log(logFile, 'Step04: Realign（头动校正）');
fprintf('\n--- Step 04: Realign ---\n');

[~, funBase_ST] = fileparts(funNii_ST);
funNii_R  = fullfile(cfg.funImgARDir, ['r' funBase_ST '.nii']);
rpFile    = fullfile(cfg.realignParamDir, sprintf('rp_%s.txt', cfg.subID));

if ~exist(funNii_R, 'file')
    [funNii_R, rpFile] = realign_estimate_reslice(...
        funNii_ST, cfg.funImgARDir, cfg.realignParamDir, cfg);
else
    write_log(logFile, '  已存在，跳过');
end
write_log(logFile, sprintf('  重采样: %s', funNii_R));
write_log(logFile, sprintf('  头动参数: %s', rpFile));

% ======================================================================
% Step 05: 重定位（Reorient: AC 为原点）
% ======================================================================
write_log(logFile, 'Step05: Reorient（坐标原点 → AC）');
fprintf('\n--- Step 05: Reorient ---\n');

% 功能像重定位
[~, funBase_R] = fileparts(funNii_R);
funNii_Ro = fullfile(cfg.funImgARDir, ['reorient_' funBase_R '.nii']);
if ~exist(funNii_Ro, 'file')
    funNii_Ro = reorient_set_origin(funNii_R, cfg.funImgARDir, []);
else
    write_log(logFile, '  功能像重定位已存在，跳过');
end

% T1 重定位
[~, t1Base] = fileparts(t1Nii);
t1Nii_Ro = fullfile(cfg.t1ImgDir, ['reorient_' t1Base '.nii']);
if ~exist(t1Nii_Ro, 'file')
    t1Nii_Ro = reorient_set_origin(t1Nii, cfg.t1ImgDir, []);
else
    write_log(logFile, '  T1 重定位已存在，跳过');
end
write_log(logFile, sprintf('  Fun Reoriented: %s', funNii_Ro));
write_log(logFile, sprintf('  T1 Reoriented:  %s', t1Nii_Ro));

% ======================================================================
% Step 06: T1 脑提取（BET）
% ======================================================================
write_log(logFile, 'Step06: T1 BET');
fprintf('\n--- Step 06: BET T1 ---\n');

[~, t1Base_Ro] = fileparts(t1Nii_Ro);
t1Nii_Bet = fullfile(cfg.t1ImgBetDir, ['bet_' t1Base_Ro '.nii']);
t1Mask_Bet = fullfile(cfg.t1ImgBetDir, ['betmask_' t1Base_Ro '.nii']);
if ~exist(t1Nii_Bet, 'file')
    [t1Nii_Bet, t1Mask_Bet] = brain_extract_t1(t1Nii_Ro, cfg.t1ImgBetDir, cfg);
else
    write_log(logFile, '  已存在，跳过');
end
write_log(logFile, sprintf('  T1 BET: %s', t1Nii_Bet));
write_log(logFile, sprintf('  BET Mask: %s', t1Mask_Bet));

% ======================================================================
% Step 07: T1 配准到功能像（Coreg T1 → Fun）
% ======================================================================
write_log(logFile, 'Step07: T1 Coreg to Fun');
fprintf('\n--- Step 07: Coreg T1 → Fun ---\n');

t1Nii_Coreg = fullfile(cfg.t1ImgCoregDir, ['coreg_' t1Base_Ro '.nii']);
if ~exist(t1Nii_Coreg, 'file')
    [t1Nii_Coreg, ~] = coreg_t1_to_fun(t1Nii_Bet, funNii_Ro, cfg.t1ImgCoregDir);
else
    write_log(logFile, '  已存在，跳过');
end
write_log(logFile, sprintf('  T1 Coreg: %s', t1Nii_Coreg));

% ======================================================================
% Step 08: 组织分割（New Segment: GM/WM/CSF）
% ======================================================================
write_log(logFile, 'Step08: T1 组织分割（GMM）');
fprintf('\n--- Step 08: Segment ---\n');

gmFile  = fullfile(cfg.t1SegDir, 'c2_t1.nii');  % GM: c2
wmFile  = fullfile(cfg.t1SegDir, 'c3_t1.nii');  % WM: c3
csfFile = fullfile(cfg.t1SegDir, 'c1_t1.nii');  % CSF: c1

if ~exist(gmFile, 'file')
    [gmFile, wmFile, csfFile] = segment_tissue(t1Nii_Coreg, cfg.t1SegDir, cfg);
else
    write_log(logFile, '  已存在，跳过');
end
write_log(logFile, sprintf('  GM: %s', gmFile));
write_log(logFile, sprintf('  WM: %s', wmFile));

% ======================================================================
% Step 09: 非线性配准（DARTEL 替代：SVF 速度场）
% ======================================================================
write_log(logFile, 'Step09: DARTEL 非线性配准');
fprintf('\n--- Step 09: DARTEL-like Warp ---\n');

flowFile = fullfile(cfg.t1SegDir, 'u_rc1_Template.nii');
if ~exist(flowFile, 'file')
    [flowFile, ~] = dartel_warp(gmFile, wmFile, cfg.t1SegDir, cfg);
else
    write_log(logFile, '  流场已存在，跳过');
end
write_log(logFile, sprintf('  流场: %s', flowFile));

% ======================================================================
% Step 10: 空间标准化（Normalize: Fun → MNI）
% ======================================================================
write_log(logFile, 'Step10: Normalize（功能像 → MNI 空间）');
fprintf('\n--- Step 10: Normalize ---\n');

[~, funBase_Ro] = fileparts(funNii_Ro);
funNii_W = fullfile(cfg.funImgARWDir, ['w' funBase_Ro '.nii']);
if ~exist(funNii_W, 'file')
    funNii_W = normalize_apply(funNii_Ro, flowFile, cfg.funImgARWDir, cfg.normalize);
else
    write_log(logFile, '  已存在，跳过');
end
write_log(logFile, sprintf('  标准化输出: %s', funNii_W));

% ======================================================================
% Step 11: 空间平滑（Smooth: 3D Gaussian）
% ======================================================================
write_log(logFile, sprintf('Step11: Smooth FWHM=[%.1f %.1f %.1f]mm', cfg.fwhm));
fprintf('\n--- Step 11: Smooth ---\n');

[~, funBase_W] = fileparts(funNii_W);
funNii_S = fullfile(cfg.funImgARWSDir, ['s' funBase_W '.nii']);
if ~exist(funNii_S, 'file')
    funNii_S = smooth_3d(funNii_W, cfg.funImgARWSDir, cfg.fwhm);
else
    write_log(logFile, '  已存在，跳过');
end
write_log(logFile, sprintf('  平滑输出: %s', funNii_S));

% ======================================================================
% Step 12: 一阶 GLM（Specify + Estimate + T-contrast）
% ======================================================================
write_log(logFile, 'Step12: 一阶 GLM');
fprintf('\n--- Step 12: 1st-level GLM ---\n');

spmMat = fullfile(cfg.firstLevelDir, 'SPM.mat');
if ~exist(spmMat, 'file')
    run_firstlevel_glm(funNii_S, rpFile, cfg.firstLevelDir, cfg);
else
    write_log(logFile, '  SPM.mat 已存在，跳过');
end
write_log(logFile, sprintf('  GLM 输出: %s', cfg.firstLevelDir));

% ======================================================================
% 完成
% ======================================================================
write_log(logFile, '=== Pipeline 全部完成 ===');
fprintf('\n=============================================================\n');
fprintf('  Pipeline 完成！\n');
fprintf('  结果目录: %s\n', cfg.firstLevelDir);
fprintf('  日志文件: %s\n', logFile);
fprintf('=============================================================\n');
end
