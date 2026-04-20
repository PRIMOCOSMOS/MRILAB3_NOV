function [flowField, warpedGM] = dartel_warp(gmFile, wmFile, outDir, cfg)
% dartel_warp - DARTEL 思路的微分同胚非线性配准（standalone 替代实现）
%
% 背景:
%   DARTEL（Diffeomorphic Anatomical Registration Through Exponentiated Lie Algebra）
%   使用李代数指数映射保证变形场的微分同胚性（可逆、无折叠）。
%
% 本实现使用：
%   - 稳态速度场（Stationary Velocity Field, SVF）近似微分同胚配准
%   - 多分辨率策略（由粗到细）
%   - 代价函数：GM+WM 概率图的均方误差 + 正则化项（速度场平滑度）
%   - 梯度下降优化速度场
%   - 积分速度场：矩阵指数法（scaling & squaring）
%
% 注意：本实现面向单被试（individual-to-template），若需群体模板，
%       可对多被试迭代调用并平均形变场。
%
% 输入:
%   gmFile  - GM 概率图 NIfTI 路径（配准源）
%   wmFile  - WM 概率图 NIfTI 路径（配准源）
%   outDir  - 输出目录
%   cfg     - 配置结构体:
%               cfg.dartel.nLevels (多分辨率层数, 默认3)
%               cfg.dartel.nIter   (各层迭代次数向量)
%               cfg.dartel.reg     (正则化系数)
%
% 输出:
%   flowField - [nx ny nz 3] 的位移场 NIfTI 文件路径
%   warpedGM  - 形变后 GM 概率图 NIfTI 文件路径

fprintf('[dartel_warp] GM: %s\n', gmFile);
fprintf('[dartel_warp] WM: %s\n', wmFile);

if use_spm_structural_backend(cfg)
    [flowField, warpedGM] = dartel_with_spm(gmFile, wmFile, outDir, cfg);
    return;
end

% -------- 读取概率图 --------
[gmData, gmHdr] = nifti_read(gmFile);
[wmData, wmHdr] = nifti_read(wmFile);

gm = double(gmData(:,:,:,1));
wm = double(wmData(:,:,:,1));
[nx, ny, nz] = size(gm);

nLevels = cfg.dartel.nLevels;
nIter   = cfg.dartel.nIter;
reg     = cfg.dartel.reg;

% -------- 加载真实模板（东亚模板或用户提供模板）--------
% 返回模板原始数据及其仿射矩阵，不再强制重采样到个体尺寸
[template_gm, template_wm, tmplHdr] = load_dartel_templates(cfg);
[tnx, tny, tnz] = size(template_gm);
fprintf('[dartel_warp] 模板尺寸: [%d %d %d]\n', tnx, tny, tnz);

% -------- 将个体 GM/WM 重采样到模板坐标空间 --------
% 使用仿射矩阵正确对齐坐标系（而非仅匹配像素数量）
% 此步骤等价于 SPM New Segment 将组织概率图输出到 MNI 空间后再进行 DARTEL 优化
fprintf('[dartel_warp] 将个体 GM/WM 重采样到模板空间...\n');
gm_in_tmpl = resample_vol_affine(gm, gmHdr.affine, tmplHdr.affine, [tnx tny tnz]);
wm_in_tmpl = resample_vol_affine(wm, wmHdr.affine, tmplHdr.affine, [tnx tny tnz]);

% -------- 初始化速度场（在模板空间定义）--------
% 速度场 v: [tnx tny tnz 3]，各分量单位为模板体素位移
v = zeros(tnx, tny, tnz, 3);

% -------- 多分辨率配准 --------
ensure_dir(outDir);
for lev = 1:nLevels
    % 当前分辨率的下采样因子
    scale = 2^(nLevels - lev);  % 从粗到细

    fprintf('[dartel_warp] 分辨率层 %d/%d (缩放1/%d)\n', lev, nLevels, scale);

    % 下采样模板空间的源图和模板（均在同一坐标系下）
    gm_s  = downsample_vol(gm_in_tmpl,  scale);
    wm_s  = downsample_vol(wm_in_tmpl,  scale);
    tgm_s = downsample_vol(template_gm, scale);
    twm_s = downsample_vol(template_wm, scale);
    v_s   = downsample_flow(v, scale);

    [sx, sy, sz, ~] = size(v_s);

    % ---- 梯度下降优化速度场 ----
    lr = 0.01 / (scale^2);  % 学习率（随分辨率调整）
    for it = 1:nIter(min(lev, numel(nIter)))
        % 积分速度场得到位移场（scaling & squaring, 简化2步）
        disp_s = integrate_svf(v_s, 4);

        % 将源图形变到模板空间
        gm_warped = warp_volume(gm_s,  disp_s);
        wm_warped = warp_volume(wm_s,  disp_s);

        % 残差（MSE 梯度）
        res_gm = gm_warped - tgm_s;
        res_wm = wm_warped - twm_s;

        % 图像梯度（在形变后空间）
        [dgx_gm, dgy_gm, dgz_gm] = gradient(gm_warped);
        [dgx_wm, dgy_wm, dgz_wm] = gradient(wm_warped);

        % 速度场梯度（数据项）
        dv = zeros(sx, sy, sz, 3);
        dv(:,:,:,1) = res_gm .* dgx_gm + res_wm .* dgx_wm;
        dv(:,:,:,2) = res_gm .* dgy_gm + res_wm .* dgy_wm;
        dv(:,:,:,3) = res_gm .* dgz_gm + res_wm .* dgz_wm;

        % 正则化项（速度场的拉普拉斯平滑）
        for d = 1:3
            dv(:,:,:,d) = dv(:,:,:,d) + reg * laplacian3d(v_s(:,:,:,d));
        end

        % 更新
        v_s = v_s - lr * dv;
    end

    % 上采样速度场回模板完整分辨率
    v = upsample_flow(v_s, scale, [tnx tny tnz]);

    mse = mean((gm_s - tgm_s).^2, 'all') + mean((wm_s - twm_s).^2, 'all');
    fprintf('[dartel_warp]   层%d 完成, MSE=%.6f\n', lev, mse);
end

% -------- 最终位移场 --------
dispFinal = integrate_svf(v, 6);  % [tnx tny tnz 3]

% -------- 写出位移场（4D NIfTI，第4维=3个方向）--------
% 关键：使用模板的仿射矩阵作为位移场的空间头信息，
% 这样 normalize_apply 可以通过读取该头信息得知模板坐标系，
% 从而正确地将 MNI 世界坐标映射到模板体素坐标，再经位移场找到个体体素。
flowHdr = tmplHdr;
flowHdr.dim = int16([4, tnx, tny, tnz, 3, 1, 1, 1]);
flowHdr.nt  = 3;
flowHdr.nx  = tnx; flowHdr.ny = tny; flowHdr.nz = tnz;
flowHdr.datatype = 16;  % float32
flowHdr.descrip = 'DARTELflow_SVF';

flowFile = fullfile(outDir, 'u_rc1_Template.nii');
nifti_write(flowFile, single(dispFinal), flowHdr);
fprintf('[dartel_warp] 位移场已写出: %s  [%d %d %d 3]\n', flowFile, tnx, tny, tnz);

% -------- 输出形变后的 GM（在模板空间中验证配准质量）--------
gm_warped_final = warp_volume(gm_in_tmpl, dispFinal);
warpedHdr = tmplHdr;
warpedHdr.nt  = 1;
warpedHdr.nx  = tnx; warpedHdr.ny = tny; warpedHdr.nz = tnz;
warpedHdr.dim = int16([3, tnx, tny, tnz, 1, 1, 1, 1]);
warpedHdr.datatype = 16;
warpedHdr.descrip = 'DARTEL_warped_GM';
warpedFile = fullfile(outDir, 'warped_gm.nii');
nifti_write(warpedFile, single(gm_warped_final), warpedHdr);
fprintf('[dartel_warp] 形变 GM 已写出: %s\n', warpedFile);

flowField = flowFile;
warpedGM  = warpedFile;
end

function [flowField, warpedGM] = dartel_with_spm(gmFile, wmFile, outDir, cfg)
% 使用 SPM DARTEL 默认 6 层优化设置
spmDir = resolve_spm_dir(cfg);
if ~exist(fullfile(spmDir, 'spm.m'), 'file')
    error('[dartel_warp] 未找到 SPM: %s', spmDir);
end

addpath(spmDir);
spm('defaults', 'FMRI');
spm_jobman('initcfg');

ensure_dir(outDir);
cleanup_dartel_outputs(outDir);
rc1File = infer_rc_file(gmFile, 1);
rc2File = infer_rc_file(wmFile, 2);

matlabbatch = build_dartel_batch(rc1File, rc2File, cfg);

spm_jobman('run', matlabbatch);

uFiles = dir(fullfile(outDir, 'u_rc1*.nii'));
if isempty(uFiles)
    error('[dartel_warp] 未找到 SPM 输出流场 u_rc1*.nii: %s', outDir);
end

[~, idx] = max([uFiles.datenum]);
flowField = fullfile(uFiles(idx).folder, uFiles(idx).name);
warpedGM = '';

fprintf('[dartel_warp][SPM] 流场: %s\n', flowField);
end

function matlabbatch = build_dartel_batch(rc1File, rc2File, cfg)
% 优先使用 DPABI 官方 Jobmat，避免参数漂移；不存在则回退 SPM 默认参数
jobmatFile = resolve_dpabi_dartel_jobmat(cfg);
if ~isempty(jobmatFile) && exist(jobmatFile, 'file')
    S = load(jobmatFile);
    if isfield(S, 'matlabbatch') && iscell(S.matlabbatch) && ~isempty(S.matlabbatch)
        matlabbatch = S.matlabbatch;
        matlabbatch{1}.spm.tools.dartel.warp.images{1,1} = {rc1File};
        matlabbatch{1}.spm.tools.dartel.warp.images{1,2} = {rc2File};
        fprintf('[dartel_warp][SPM] 使用 DPABI Jobmat: %s\n', jobmatFile);
        return;
    end
    warning('[dartel_warp] DPABI Jobmat 缺少 matlabbatch，回退默认参数: %s', jobmatFile);
end

matlabbatch = {};
matlabbatch{1}.spm.tools.dartel.warp.images{1,1} = {rc1File};
matlabbatch{1}.spm.tools.dartel.warp.images{1,2} = {rc2File};
matlabbatch{1}.spm.tools.dartel.warp.settings.template = 'Template';
matlabbatch{1}.spm.tools.dartel.warp.settings.rform = 0;

matlabbatch{1}.spm.tools.dartel.warp.settings.param(1).its = 3;
matlabbatch{1}.spm.tools.dartel.warp.settings.param(1).rparam = [4 2 1e-06];
matlabbatch{1}.spm.tools.dartel.warp.settings.param(1).K = 0;
matlabbatch{1}.spm.tools.dartel.warp.settings.param(1).slam = 16;

matlabbatch{1}.spm.tools.dartel.warp.settings.param(2).its = 3;
matlabbatch{1}.spm.tools.dartel.warp.settings.param(2).rparam = [2 1 1e-06];
matlabbatch{1}.spm.tools.dartel.warp.settings.param(2).K = 0;
matlabbatch{1}.spm.tools.dartel.warp.settings.param(2).slam = 8;

matlabbatch{1}.spm.tools.dartel.warp.settings.param(3).its = 3;
matlabbatch{1}.spm.tools.dartel.warp.settings.param(3).rparam = [1 0.5 1e-06];
matlabbatch{1}.spm.tools.dartel.warp.settings.param(3).K = 1;
matlabbatch{1}.spm.tools.dartel.warp.settings.param(3).slam = 4;

matlabbatch{1}.spm.tools.dartel.warp.settings.param(4).its = 3;
matlabbatch{1}.spm.tools.dartel.warp.settings.param(4).rparam = [0.5 0.25 1e-06];
matlabbatch{1}.spm.tools.dartel.warp.settings.param(4).K = 2;
matlabbatch{1}.spm.tools.dartel.warp.settings.param(4).slam = 2;

matlabbatch{1}.spm.tools.dartel.warp.settings.param(5).its = 3;
matlabbatch{1}.spm.tools.dartel.warp.settings.param(5).rparam = [0.25 0.125 1e-06];
matlabbatch{1}.spm.tools.dartel.warp.settings.param(5).K = 4;
matlabbatch{1}.spm.tools.dartel.warp.settings.param(5).slam = 1;

matlabbatch{1}.spm.tools.dartel.warp.settings.param(6).its = 3;
matlabbatch{1}.spm.tools.dartel.warp.settings.param(6).rparam = [0.25 0.125 1e-06];
matlabbatch{1}.spm.tools.dartel.warp.settings.param(6).K = 6;
matlabbatch{1}.spm.tools.dartel.warp.settings.param(6).slam = 0.5;

matlabbatch{1}.spm.tools.dartel.warp.settings.optim.lmreg = 0.01;
matlabbatch{1}.spm.tools.dartel.warp.settings.optim.cyc = 3;
matlabbatch{1}.spm.tools.dartel.warp.settings.optim.its = 3;
fprintf('[dartel_warp][SPM] 使用内置默认 DARTEL 参数\n');
end

function jobmatFile = resolve_dpabi_dartel_jobmat(cfg)
jobmatFile = '';
dpabiRoot = '';

if isstruct(cfg) && isfield(cfg, 'installPaths') && isfield(cfg.installPaths, 'dpabiRoot') && ...
        ~isempty(cfg.installPaths.dpabiRoot)
    dpabiRoot = cfg.installPaths.dpabiRoot;
end
if isempty(dpabiRoot)
    dpabiRoot = 'D:/DPABI_V9.0_250415';
end

cand = fullfile(dpabiRoot, 'DPARSF', 'Jobmats', 'Dartel_CreateTemplate.mat');
if exist(cand, 'file')
    jobmatFile = cand;
end
end

function cleanup_dartel_outputs(outDir)
patterns = {'Template_*.nii', 'Template_*.mat', 'u_rc1*.nii'};
for i = 1:numel(patterns)
    files = dir(fullfile(outDir, patterns{i}));
    for k = 1:numel(files)
        try
            delete(fullfile(files(k).folder, files(k).name));
        catch
            % Ignore locked files; SPM will overwrite when possible.
        end
    end
end
end

function rcFile = infer_rc_file(cFile, tissueIdx)
[pth, name, ext] = fileparts(cFile);
if startsWith(name, sprintf('c%d', tissueIdx))
    rcName = regexprep(name, sprintf('^c%d', tissueIdx), sprintf('rc%d', tissueIdx));
    rcFile = fullfile(pth, [rcName ext]);
else
    rcFile = fullfile(pth, sprintf('rc%d%s', tissueIdx, name));
end

if ~exist(rcFile, 'file')
    cand = dir(fullfile(pth, sprintf('rc%d*.nii', tissueIdx)));
    if ~isempty(cand)
        rcFile = fullfile(cand(1).folder, cand(1).name);
    end
end

if ~exist(rcFile, 'file')
    error('[dartel_warp] 未找到 DARTEL 导入文件 rc%d*: %s', tissueIdx, pth);
end
end

function tf = use_spm_structural_backend(cfg)
tf = isstruct(cfg) && isfield(cfg, 'spm') && ...
     isfield(cfg.spm, 'useStructural') && logical(cfg.spm.useStructural);
end

function spmDir = resolve_spm_dir(cfg)
if isstruct(cfg) && isfield(cfg, 'spm') && isfield(cfg.spm, 'dir') && ~isempty(cfg.spm.dir)
    spmDir = cfg.spm.dir;
elseif isstruct(cfg) && isfield(cfg, 'installPaths') && isfield(cfg.installPaths, 'spmRoot')
    spmDir = cfg.installPaths.spmRoot;
else
    spmDir = 'D:/spm';
end
end

% ======================================================================
% 内部辅助函数
% ======================================================================

function disp = integrate_svf(v, nSteps)
% 将稳态速度场 v 通过 scaling & squaring 积分为位移场
% v: [nx ny nz 3]
% nSteps: 积分步数（scaling factor = 2^nSteps）
disp = v / (2^nSteps);
for s = 1:nSteps
    disp = compose_flow(disp, disp);
end
end

function phi2 = compose_flow(phi1, phi2)
% 合成两个位移场: φ2(x) = φ1(x + φ2(x)) + φ2(x)
[nx, ny, nz, ~] = size(phi1);
[Xg, Yg, Zg] = ndgrid(1:nx, 1:ny, 1:nz);

% 变形后坐标
Xd = Xg + phi2(:,:,:,1);
Yd = Yg + phi2(:,:,:,2);
Zd = Zg + phi2(:,:,:,3);

coords = [Xd(:)'; Yd(:)'; Zd(:)'];

% 在 phi1 中插值（逐方向）
phi2_out = phi2;
for d = 1:3
    phi1_d = phi1(:,:,:,d);
    phi2_out(:,:,:,d) = phi2(:,:,:,d) + ...
        reshape(trilinear_interp(phi1_d, coords), nx, ny, nz);
end
phi2 = phi2_out;
end

function Vw = warp_volume(V, disp)
% 用位移场 disp 对体数据 V 进行前向形变采样（pull interpolation）
% 对每个输出体素，其采样位置 = 当前网格坐标 + 位移
% 这等价于 "将输出网格中的点拉取（pull）到输入图像中的对应位置"
[nx, ny, nz] = size(V);
[Xg, Yg, Zg] = ndgrid(1:nx, 1:ny, 1:nz);
Xd = Xg + disp(:,:,:,1);
Yd = Yg + disp(:,:,:,2);
Zd = Zg + disp(:,:,:,3);
Vw = reshape(trilinear_interp(V, [Xd(:)'; Yd(:)'; Zd(:)']), nx, ny, nz);
end

function V_down = downsample_vol(V, scale)
% 简单平均池化下采样
if scale == 1, V_down = V; return; end
[nx, ny, nz] = size(V);
nx2 = floor(nx/scale);
ny2 = floor(ny/scale);
nz2 = floor(nz/scale);
V_down = zeros(nx2, ny2, nz2);
for x = 1:nx2
    for y = 1:ny2
        for z = 1:nz2
            xs = (x-1)*scale+1 : min(x*scale, nx);
            ys = (y-1)*scale+1 : min(y*scale, ny);
            zs = (z-1)*scale+1 : min(z*scale, nz);
            V_down(x,y,z) = mean(V(xs,ys,zs),'all');
        end
    end
end
end

function v_down = downsample_flow(v, scale)
if scale == 1, v_down = v; return; end
v_down = zeros(floor(size(v,1)/scale), floor(size(v,2)/scale), ...
               floor(size(v,3)/scale), 3);
for d = 1:3
    v_down(:,:,:,d) = downsample_vol(v(:,:,:,d), scale) / scale;
end
end

function v_up = upsample_flow(v_s, scale, targetSize)
if scale == 1, v_up = v_s; return; end
[sx, sy, sz, ~] = size(v_s);
nx = targetSize(1); ny = targetSize(2); nz = targetSize(3);
v_up = zeros(nx, ny, nz, 3);
[Xg, Yg, Zg] = ndgrid(1:nx, 1:ny, 1:nz);
coordsSrc = [Xg(:)'/scale; Yg(:)'/scale; Zg(:)'/scale];
for d = 1:3
    v_up(:,:,:,d) = scale * reshape(trilinear_interp(v_s(:,:,:,d), coordsSrc), nx, ny, nz);
end
end

function L = laplacian3d(V)
% 3D 拉普拉斯算子（6-邻域有限差分）
L = -6*V;
L(2:end,:,:)   = L(2:end,:,:)   + V(1:end-1,:,:);
L(1:end-1,:,:) = L(1:end-1,:,:) + V(2:end,:,:);
L(:,2:end,:)   = L(:,2:end,:)   + V(:,1:end-1,:);
L(:,1:end-1,:) = L(:,1:end-1,:) + V(:,2:end,:);
L(:,:,2:end)   = L(:,:,2:end)   + V(:,:,1:end-1);
L(:,:,1:end-1) = L(:,:,1:end-1) + V(:,:,2:end);
end


function [template_gm, template_wm, tmplHdr] = load_dartel_templates(cfg)
% load_dartel_templates - 加载 DARTEL 模板并返回原始数据及其仿射矩阵
% 支持两种配置方式：4D 单文件模板 或 GM/WM 双文件模板
% 不再对模板进行尺寸重采样；由调用方负责将个体图像重采样到模板空间
if ~isfield(cfg, 'templates') || ~isfield(cfg.templates, 'dartel')
    error('[dartel_warp] 缺少 cfg.templates.dartel 配置');
end

has4DTemplate = isfield(cfg.templates.dartel, 'template4DNii') && ...
                ~isempty(cfg.templates.dartel.template4DNii);
hasDualFiles = isfield(cfg.templates.dartel, 'gmTemplateNii') && ...
               isfield(cfg.templates.dartel, 'wmTemplateNii') && ...
               ~isempty(cfg.templates.dartel.gmTemplateNii) && ...
               ~isempty(cfg.templates.dartel.wmTemplateNii);

if has4DTemplate
    template4DFile = cfg.templates.dartel.template4DNii;
    if ~exist(template4DFile, 'file')
        error('[dartel_warp] DARTEL 4D模板不存在: %s', template4DFile);
    end
    if ~isfield(cfg.templates.dartel, 'gmVolumeIndex') || ~isfield(cfg.templates.dartel, 'wmVolumeIndex')
        error('[dartel_warp] 4D模板模式下必须提供 gmVolumeIndex / wmVolumeIndex');
    end
    gmIdx = cfg.templates.dartel.gmVolumeIndex;
    wmIdx = cfg.templates.dartel.wmVolumeIndex;
    if any([gmIdx, wmIdx] < 1) || any(mod([gmIdx, wmIdx], 1) ~= 0)
        error('[dartel_warp] gmVolumeIndex / wmVolumeIndex 必须为正整数');
    end
    if gmIdx == wmIdx
        error('[dartel_warp] gmVolumeIndex 与 wmVolumeIndex 不能相同');
    end

    [template4D, tmplHdr] = nifti_read(template4DFile);
    nVols = size(template4D, 4);
    if gmIdx > nVols || wmIdx > nVols
        error('[dartel_warp] 4D模板帧数不足: 需要 GM=%d WM=%d, 实际=%d', gmIdx, wmIdx, nVols);
    end

    fprintf('[dartel_warp] 使用4D模板: %s (GM=%d, WM=%d)\n', template4DFile, gmIdx, wmIdx);
    template_gm = double(template4D(:,:,:,gmIdx));
    template_wm = double(template4D(:,:,:,wmIdx));
    return;
end

if ~hasDualFiles
    error('[dartel_warp] 未提供有效DARTEL模板：请配置4D模板或GM/WM双文件模板');
end

gmTemplateFile = cfg.templates.dartel.gmTemplateNii;
wmTemplateFile = cfg.templates.dartel.wmTemplateNii;
if ~exist(gmTemplateFile, 'file') || ~exist(wmTemplateFile, 'file')
    error('[dartel_warp] 模板文件不存在，请检查配置: GM=%s WM=%s', gmTemplateFile, wmTemplateFile);
end
fprintf('[dartel_warp] 使用双文件模板: GM=%s, WM=%s\n', gmTemplateFile, wmTemplateFile);
[template_gm_raw, tmplHdr] = nifti_read(gmTemplateFile);
[template_wm_raw, ~]       = nifti_read(wmTemplateFile);
template_gm = double(template_gm_raw(:,:,:,1));
template_wm = double(template_wm_raw(:,:,:,1));
end
