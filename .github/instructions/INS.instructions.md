---
description: fMRI-BOLD 数据分析Pipeline双阶段实现：第一阶段通过代码调用SPM25和DPABI复现Gold Standard；第二阶段再移植核心算法实现Standalone，并达到高相似度。
# applyTo: 'Describe when these instructions should be loaded by the agent based on task context' # when provided, instructions will automatically be added to the request context when the pattern matches an attached file
---

<!-- Tip: Use /create-instructions in chat to generate content with agent assistance -->

fMRI-BOLD 数据分析Pipeline采用双阶段目标推进：
第一阶段，优先通过代码调用SPM25和DPABI（含DPARSF流程）复现完整pipeline，确认可以稳定复现Gold Standard结果；
第二阶段，在第一阶段复现成功后，再逐步拆解并移植核心算法与处理逻辑，最终实现Standalone版本并达到与Gold Standard的高相似度。

本地具有MATLAB，SPM25和DPABI的环境可以调用，也有由原数据通过运行这两个软件得到的示例处理结果。MATLAB环境已配置MCP; SPM25和DPABI的路径分别是：
SPM25："D:\Coding2"；DPABI："D:\DPABI_V9.0_250415"。同时，我的示例数据(Gold Standard)位于"D:\MRI_PRO\MRILAB3\BOLDDATA"下，其中文件夹F1Raw和FunRaw是原始的数据。
"D:\MRI_PRO\MRILAB3\BOLDCODE\BOLDDATA"下是本库运行得到的测试数据处理结果，也是我们优化的目标。

我们的大致思路分两步：
第一步是“调用复现”。优先编写可重复执行的MATLAB脚本，以代码方式调用SPM25和DPABI/DPARSF，完整跑通与Gold Standard一致的处理链路，固定版本、参数、路径和输出组织方式，先建立可验证的Gold基线；
第二步是“移植替换”。在Gold基线稳定后，再分析SPM25和DPABI中DPARSF的核心算法与实现细节，按模块逐步替换为本库自实现逻辑。替换过程中持续进行阶段性相似度审计，确保每次替换不破坏与Gold基线的一致性。

本库中名字带EXAMPLE的文件夹可用于参考，也可在第一阶段用于验证调用链路；但在第二阶段的Standalone目标中，这些开源实现不能作为最终核心逻辑依赖（除文件格式转换可使用dcm2nii/dcm2niix工具外）。当相似性指标偏低时，分析不能只停留在数据层面，而要追溯到算法与实现细节进行修正。

另外，我们的pipeline在逐步的数据处理结果表示形式和其他约定上可能也存在差异，我们要在这方面也寻求和SPM25以及DPABI寻求全部对齐，以免我们的一致性评估被这些差异所干扰。最后的T-contrast与gold standard的相似度必须要在0.95以上才有讨论的价值。

供参考pipeline相关内容：
# fMRI 数据预处理 Pipeline 文件生成顺序说明

以下是基于SPM/DPARSF预处理流程生成文件的严格先后顺序，划分为数据准备、结构像（T1）处理、功能像（BOLD）处理以及质控与统计四大阶段：

## 1. 第0阶段：数据准备与配置 (Data Prep)
- **DPARSFA_AutoSave_[日期].mat**：存放运行此预处理流程前自动保存的参数配置环境快照。
- **FunRaw**：存放未经任何处理的原始功能磁共振（BOLD）影像数据。
- **T1Raw**：存放未经任何处理的原始T1加权高分辨率结构磁共振影像数据。
- **ReorientMats**：存放对原始T1或功能像进行初始位置校正（如对齐前连合）时生成的原点重定向变换矩阵。

## 2. 第1阶段：结构像预处理 (T1 Pipeline 顺序)
*注：结构像处理的主要目的是为功能像提供精准的解剖定位和标准空间映射参数。*
1. **T1Img**：存放从T1Raw转换而来、格式统一且初步整理后的T1结构像。
2. **T1ImgBet**：存放经过颅骨剥离操作（剔除非脑组织）后的干净T1结构像。
3. **T1ImgCoreg**：存放与功能像（通常是头动校正后的平均帧）完成空间严格配准（Coregistration）的T1结构像。
4. **T1ImgNewSegment**：存放配准后的T1图像被分割为灰质/白质/脑脊液的结果，同时生成**供后续功能像使用**的向标准空间映射的形变场文件。

## 3. 第2阶段：功能像预处理 (BOLD Pipeline 顺序)
*注：功能像处理严格遵循累加后缀（A->R->W->S）的先后顺序生成。*
1. **FunImg**：存放从FunRaw转换而来，并剔除序列初始不稳定时间点（Dummy Scans）后的功能像。
2. **FunImgA**：存放经过时间层校正（Slice Timing），消除了同一TR内各切片时间差的功能像。
3. **RealignParameter**：存放在此阶段计算出的头动校正（Realignment）平移与旋转六参数。
4. **FunImgAR**：存放应用上述参数完成头动校正，将所有时间点对齐到同一位置的功能像。
5. **FunImgARW**：存放应用结构像（`T1ImgNewSegment`中）生成的形变场，被非线性拉伸配准到MNI标准空间（Warping/Normalization）的功能像。
6. **FunImgARWS**：存放经过高斯核空间平滑（Smoothing）处理，提高信噪比后的最终可用功能像。

## 4. 第3阶段：质控与第一阶统计 (QC & 1st-Level Stats)
- **Masks**：存放根据结构像分割结果生成的脑实质掩膜，用于限制后续统计计算的区域。
- **PicturesForChkNormalization**：存放空间标准化后生成的功能像与T1像对齐效果截图，供人工目视检查。
- **QC**：存放提取出的全脑信号、头动相对位移（FD）等用于整体数据质量评估的报表。
- **Sub01_1stLevel**：存放基于预处理终产物（`FunImgARWS`）完成的单被试个体水平一般线性模型（GLM）统计结果。
- **Sub01_1stLevel_StandaloneCheck**：存放个体水平统计模型设置或设计矩阵的独立验证与辅助检查文件。

当前主要目标：先尝试通过代码调用SPM25和DPABI功能复现整个pipeline，确认可复现Gold Standard；复现稳定后，再进入算法/逻辑移植实现阶段，最终达到与Gold Standard的高相似度。