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

% -------- 读取概率图 --------
[gmData, gmHdr] = nifti_read(gmFile);
[wmData, ~]     = nifti_read(wmFile);

gm = double(gmData(:,:,:,1));
wm = double(wmData(:,:,:,1));
[nx, ny, nz] = size(gm);

nLevels = cfg.dartel.nLevels;
nIter   = cfg.dartel.nIter;
reg     = cfg.dartel.reg;

% -------- 加载真实模板（东亚模板或用户提供模板）--------
[template_gm, template_wm] = load_dartel_templates(cfg, [nx ny nz]);

% -------- 初始化速度场 --------
% 速度场 v: [nx ny nz 3]，单位为体素位移
v = zeros(nx, ny, nz, 3);

% -------- 多分辨率配准 --------
ensure_dir(outDir);
for lev = 1:nLevels
    % 当前分辨率的下采样因子
    scale = 2^(nLevels - lev);  % 从粗到细

    fprintf('[dartel_warp] 分辨率层 %d/%d (缩放1/%d)\n', lev, nLevels, scale);

    % 下采样源图和模板
    gm_s  = downsample_vol(gm,          scale);
    wm_s  = downsample_vol(wm,          scale);
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

    % 上采样速度场回原分辨率
    if lev < nLevels
        v = upsample_flow(v_s, scale, [nx ny nz]);
    else
        v = upsample_flow(v_s, scale, [nx ny nz]);
    end

    mse = mean((gm_s - tgm_s).^2) + mean((wm_s - twm_s).^2);
    fprintf('[dartel_warp]   层%d 完成, MSE=%.6f\n', lev, mse);
end

% -------- 最终位移场 --------
dispFinal = integrate_svf(v, 6);  % [nx ny nz 3]

% -------- 写出位移场（4D NIfTI，第4维=3个方向）--------
flowHdr = gmHdr;
flowHdr.dim = int16([4, nx, ny, nz, 3, 1, 1, 1]);
flowHdr.nt  = 3;
flowHdr.descrip = 'DARTELflow_SVF';

flowFile = fullfile(outDir, 'u_rc1_Template.nii');
nifti_write(flowFile, single(dispFinal), flowHdr);
fprintf('[dartel_warp] 位移场已写出: %s\n', flowFile);

% -------- 输出形变后的 GM --------
gm_warped_final = warp_volume(gm, dispFinal);
warpedHdr = gmHdr;
warpedHdr.descrip = 'DARTEL_warped_GM';
warpedFile = fullfile(outDir, 'warped_gm.nii');
nifti_write(warpedFile, single(gm_warped_final), warpedHdr);
fprintf('[dartel_warp] 形变 GM 已写出: %s\n', warpedFile);

flowField = flowFile;
warpedGM  = warpedFile;
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

function V2 = resample_vol_to_size(V, targetSize)
% 三线性重采样到目标尺寸（仅尺寸匹配，不改变强度统计趋势）
[nx, ny, nz] = size(V);
tx = targetSize(1); ty = targetSize(2); tz = targetSize(3);
if nx == tx && ny == ty && nz == tz
    V2 = V;
    return;
end
[Xt, Yt, Zt] = ndgrid(1:tx, 1:ty, 1:tz);
Xs = (Xt - 1) * (nx - 1) / max(tx - 1, 1) + 1;
Ys = (Yt - 1) * (ny - 1) / max(ty - 1, 1) + 1;
Zs = (Zt - 1) * (nz - 1) / max(tz - 1, 1) + 1;
V2 = reshape(trilinear_interp(V, [Xs(:)'; Ys(:)'; Zs(:)']), tx, ty, tz);
end

function [template_gm, template_wm] = load_dartel_templates(cfg, targetSize)
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

    [template4D, ~] = nifti_read(template4DFile);
    nVols = size(template4D, 4);
    if gmIdx > nVols || wmIdx > nVols
        error('[dartel_warp] 4D模板帧数不足: 需要 GM=%d WM=%d, 实际=%d', gmIdx, wmIdx, nVols);
    end

    fprintf('[dartel_warp] 使用4D模板: %s (GM=%d, WM=%d)\n', template4DFile, gmIdx, wmIdx);
    template_gm = resample_vol_to_size(double(template4D(:,:,:,gmIdx)), targetSize);
    template_wm = resample_vol_to_size(double(template4D(:,:,:,wmIdx)), targetSize);
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
[template_gm_raw, ~] = nifti_read(gmTemplateFile);
[template_wm_raw, ~] = nifti_read(wmTemplateFile);
template_gm = resample_vol_to_size(double(template_gm_raw(:,:,:,1)), targetSize);
template_wm = resample_vol_to_size(double(template_wm_raw(:,:,:,1)), targetSize);
end
