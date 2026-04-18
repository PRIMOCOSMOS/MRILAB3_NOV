function render_activation_3d(tMapFile, brainTemplateFile, outDir, visCfg)
% render_activation_3d - 交互式3D激活渲染 (替代 SPM Renderer 的现代化展示)
% 输入:
%   tMapFile         - T统计图 NIfTI（如 spmT_*.nii）
%   brainTemplateFile- 标准脑模板 NIfTI（如 MNI152_T1_2mm）
%   outDir           - 输出目录
%   visCfg           - 可视化配置结构体（阈值、透明度、是否输出png）

if ~exist(tMapFile, 'file')
    error('[render_activation_3d] T图不存在: %s', tMapFile);
end
if ~exist(brainTemplateFile, 'file')
    error('[render_activation_3d] 脑模板不存在: %s', brainTemplateFile);
end

% SPM Renderer 逻辑解析（仅做兼容记录）:
% spm_render(dat, brt, rendfile) 的 rendfile 常为 rend/render_single_subj.mat，
% 内含多视角投影所需的渲染矩阵与深度图（不是可直接旋转的3D网格）。
% 本函数采用现代化交互3D（isosurface + rotate3d），因此不直接依赖该 .mat，
% 但如果配置了该文件会进行存在性检查并输出提示，便于追踪与开源逻辑对齐。
if isfield(visCfg, 'spmRenderTemplateMat') && ~isempty(visCfg.spmRenderTemplateMat)
    if exist(visCfg.spmRenderTemplateMat, 'file')
        fprintf('[render_activation_3d] 检测到 SPM Renderer 模板: %s\n', visCfg.spmRenderTemplateMat);
    else
        warning('[render_activation_3d] 配置了 SPM Renderer 模板但文件不存在: %s', visCfg.spmRenderTemplateMat);
    end
end

[tData, tHdr]     = nifti_read(tMapFile);
[brainData, brainHdr] = nifti_read(brainTemplateFile);
tMap   = double(tData(:,:,:,1));
brain  = double(brainData(:,:,:,1));

% 若尺寸不同，使用仿射矩阵对齐重采样，将统计图重采样到脑模板空间
if any(size(tMap) ~= size(brain))
    tMap = resample_vol_affine(tMap, tHdr.affine, brainHdr.affine, size(brain));
end

tThr = visCfg.tThreshold;
alphaBrain = visCfg.alphaBrain;
alphaAct = visCfg.alphaActivation;

% 构建脑与激活体素等值面
brainLevel = prctile(brain(brain>0), 40);
actMask = tMap >= tThr;
if ~any(actMask(:))
    warning('[render_activation_3d] 没有超过阈值 T>=%.2f 的体素，跳过3D渲染', tThr);
    return;
end

fig = figure('Name','3D Activation Renderer','Color','k');
ax = axes('Parent', fig);
hold(ax, 'on');
axis(ax, 'equal');
axis(ax, 'off');
view(ax, 3);
camlight(ax, 'headlight');
lighting(ax, 'gouraud');

% 脑壳
fvBrain = isosurface(brain, brainLevel);
if ~isempty(fvBrain.vertices)
    pBrain = patch(ax, fvBrain);
    pBrain.FaceColor = [0.75 0.75 0.78];
    pBrain.EdgeColor = 'none';
    pBrain.FaceAlpha = alphaBrain;
end

% 激活区
actVol = zeros(size(tMap));
actVol(actMask) = tMap(actMask);
fvAct = isosurface(actVol, tThr);
if ~isempty(fvAct.vertices)
    pAct = patch(ax, fvAct);
    pAct.FaceVertexCData = fvAct.vertices(:,3);
    pAct.FaceColor = 'interp';
    pAct.EdgeColor = 'none';
    pAct.FaceAlpha = alphaAct;
    colormap(ax, hot(256));
    colorbar(ax, 'Color', 'w');
end

rotate3d(ax, 'on'); % 交互旋转
title(ax, sprintf('Interactive 3D Activation (T >= %.2f)', tThr), 'Color', 'w');

if isfield(visCfg, 'outputPng') && visCfg.outputPng
    ensure_dir(outDir);
    outPng = fullfile(outDir, 'Renderer3D_Activation.png');
    exportgraphics(fig, outPng, 'Resolution', 200);
    fprintf('[render_activation_3d] 已导出截图: %s\n', outPng);
end
end
