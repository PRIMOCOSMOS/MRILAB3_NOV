# MRILAB3_NOV

完全 Standalone 的 MATLAB BOLD-fMRI 分析 Pipeline，不依赖 SPM / DPABI。

---

## 目录结构

```
MRILAB3_NOV/
├── run_pipeline_sub01.m          % 主入口（一键运行）
├── config_sub01.m                % 参数配置文件
│
├── utils/                        % 通用工具函数
│   ├── ensure_dir.m              % 目录创建
│   ├── write_log.m               % 日志记录
│   ├── nifti_read.m              % NIfTI 读取（standalone）
│   ├── nifti_write.m             % NIfTI 写出（standalone）
│   ├── nifti_default_hdr.m       % 默认 NIfTI 头构建
│   ├── trilinear_interp.m        % 3D 三线性插值
│   ├── rigid_mat.m               % 刚体变换矩阵（6参数）
│   └── validate_pipeline_config.m% 开箱即用配置/模板校验
│
├── io/                           % 格式转换
│   ├── dicom2nifti_fun.m         % EPI MOSAIC DICOM → 4D NIfTI
│   └── dicom2nifti_t1.m          % T1 DICOM → 3D NIfTI
│
├── preprocess/                   % 预处理各步骤
│   ├── remove_dummy_tr.m         % 去除起始不稳定 TR
│   ├── slice_timing_corr.m       % 切片时序校正（傅里叶时移）
│   ├── realign_estimate_reslice.m% 头动校正（Gauss-Newton）+ rp_*.txt
│   ├── reorient_set_origin.m     % 坐标原点重定位到 AC
│   ├── coreg_t1_to_fun.m         % T1 配准到功能像（互信息）
│   ├── segment_tissue.m          % GMM 组织分割（GM/WM/CSF）
│   ├── dartel_warp.m             % 非线性配准（SVF/DARTEL 替代）
│   ├── normalize_apply.m         % 应用形变场到 MNI 空间
│   └── smooth_3d.m               % 3D 高斯平滑
│
├── stats/                        % 统计分析
│   ├── hrf_canonical.m           % 双伽马 HRF 函数
│   ├── build_design_matrix.m     % GLM 设计矩阵（HRF卷积+漂移）
│   ├── glm_ols.m                 % OLS 参数估计（β̂, σ²）
│   ├── compute_tcontrast.m       % T-contrast 统计图
│   └── run_firstlevel_glm.m      % 一阶 GLM 主函数
│
├── visualize/                    % 结果可视化
│   └── render_activation_3d.m    % 交互式3D激活渲染（可旋转）
│
└── register/                     % 配准（由 preprocess 调用）
```

---

## 快速开始

### 1. 配置参数

编辑 `config_sub01.m`，修改以下关键参数（全部集中在一个文件中）：

```matlab
cfg.funRawDir = 'D:\MRI_PRO\MRILAB3\BOLDCODE\BOLDDATA\FunRaw\Sub_01';
cfg.t1RawDir  = 'D:\MRI_PRO\MRILAB3\BOLDCODE\BOLDDATA\T1Raw\Sub_01';
cfg.TR        = 2.0;      % 重复时间（秒）
cfg.nSlices   = 36;       % EPI 层数
cfg.sliceSize = 64;       % 每层体素边长
cfg.nDummy    = 10;       % 去除的起始 TR 数
cfg.fwhm      = [6 6 6];  % 平滑核 FWHM（mm）
```

还需配置模板绝对路径（必须）：

```matlab
cfg.templates.dartel.gmTemplateNii = 'D:\...\Template_GM.nii';
cfg.templates.dartel.wmTemplateNii = 'D:\...\Template_WM.nii';
cfg.templates.standard.brainMaskNii = 'D:\...\BrainMask_2mm.nii';
cfg.templates.standard.t1TemplateNii = 'D:\...\MNI152_T1_2mm.nii';
cfg.visualization.brainTemplateNii = cfg.templates.standard.t1TemplateNii;
```

> 启动时会执行 `validate_pipeline_config`，模板缺失会直接报错退出（fail-fast）。

### 2. 运行 Pipeline

在 MATLAB 命令窗口中：

```matlab
cd('path/to/MRILAB3_NOV');
run_pipeline_sub01;
```

所有步骤会自动按序执行，中间产物保存到各输出目录。若某步骤输出已存在则跳过。

---

## 处理流程（11步）

| 步骤 | 操作 | 输入 | 输出目录 |
|------|------|------|----------|
| 01 | DICOM → NIfTI | FunRaw / T1Raw | FunImg / T1Img |
| 02 | 去除 Dummy TR | FunImg | FunImgA |
| 03 | Slice Timing 校正 | FunImgA | FunImgA（st前缀）|
| 04 | 头动校正（Realign）| FunImgA | FunImgAR + RealignParameter |
| 05 | 重定位（Reorient）| FunImgAR / T1Img | 原目录（reorient前缀）|
| 06 | T1 配准到 Fun | T1Img | T1ImgCoreg |
| 07 | 组织分割 | T1ImgCoreg | T1ImgNewSegment |
| 08 | DARTEL 非线性配准 | T1ImgNewSegment | T1ImgNewSegment（流场）|
| 09 | 标准化（→ MNI）| FunImgAR | FunImgARW |
| 10 | 空间平滑 | FunImgARW | FunImgARWS |
| 11 | 一阶 GLM | FunImgARWS | Sub01_1stLevel |
| 12 | 交互式3D渲染 | spmT_*.nii + 标准脑模板 | Sub01_1stLevel |

---

## 输出文件说明

- `FunImgARWS/s*.nii` — 预处理完成的4D功能像（MNI空间，平滑后）
- `RealignParameter/rp_Sub_01.txt` — 头动参数（nT×6，单位：mm/rad）
- `T1ImgNewSegment/c1_t1.nii` — CSF 概率图
- `T1ImgNewSegment/c2_t1.nii` — GM 概率图
- `T1ImgNewSegment/c3_t1.nii` — WM 概率图
- 注意：本项目使用强度升序编号（c1=CSF,c2=GM,c3=WM），与SPM默认编号不同
- `T1ImgNewSegment/u_rc1_Template.nii` — 非线性位移场（4D, 第4维=3方向）
- `Sub01_1stLevel/SPM.mat` — GLM 结果（MATLAB 变量，SPM 命名兼容）
- `Sub01_1stLevel/spmT_*.nii` — T 统计图
- `Sub01_1stLevel/neg_log10p_*.nii` — -log10(p) 统计图
- `Sub01_1stLevel/beta_*.nii` — 各回归参数估计图
- `Sub01_1stLevel/design_matrix.png` — 设计矩阵可视化
- `Sub01_1stLevel/Renderer3D_Activation.png` — 3D激活截图（可选）

---

## 算法说明

| 模块 | 算法 |
|------|------|
| DICOM→NIfTI | MATLAB内置dicominfo/dicomread，MOSAIC自动解码 |
| Slice Timing | 傅里叶时移定理（FFT相移+IFFT） |
| Realign | Gauss-Newton优化，6-DOF刚体，三线性插值重采样 |
| Coreg | 互信息（MI）最大化，梯度上升 |
| Segment | GMM + EM算法，K-means初始化，可选MRF正则化 |
| DARTEL | 稳态速度场（SVF）+ scaling&squaring积分 + 多分辨率 |
| Normalize | 逆变换采样（push-pull），三线性插值 |
| Smooth | 可分离3D高斯卷积（3次1D卷积） |
| GLM | OLS：β̂=(X'X)⁻¹X'Y，双伽马HRF，DCT漂移基函数 |
| T-contrast | t=c'β̂/√(σ̂²c'(X'X)⁻¹c)，p值（t分布CDF） |
| 3D Renderer | isosurface+patch+camlight+rotate3d 的交互式体渲染 |

---

## 依赖

- MATLAB R2018b 或更高版本
- **无需** SPM / DPABI / FSL / FreeSurfer 等工具箱
- 仅使用 MATLAB 内置函数（dicominfo, fft, conv, gradient 等）
