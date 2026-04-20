function outFile = normalize_apply(inFile, flowFile, outDir, mniCfg)
% normalize_apply - 将 4D 功能像应用形变场，标准化到 MNI 空间
%
% 背景:
%   空间标准化（Normalize）将个体脑图像"扭曲"变形到公共参考空间（MNI），
%   使不同被试的脑区在空间上对齐，从而支持群体级统计分析。
%
% 方法:
%   对 MNI 目标网格中的每个体素，通过形变场（位移场）
%   找到其在个体空间中的对应位置，然后三线性插值采样。
%   （逆变换插值，避免采样空洞）
%
% 输入:
%   inFile   - 输入 4D NIfTI 文件（个体空间）
%   flowFile - 位移场 NIfTI 文件（由 dartel_warp 生成, [nx ny nz 3]）
%   outDir   - 输出目录
%   mniCfg   - 目标网格配置（两种格式均支持）:
%              A) DPABI风格:
%                 mniCfg.boundingBox = [xmin ymin zmin; xmax ymax zmax]（mm）
%                 mniCfg.voxSize     = [vx vy vz]（mm）
%              B) 传统网格:
%                 mniCfg.dims    = [nx ny nz]
%                 mniCfg.voxSize = [vx vy vz]
%                 mniCfg.origin  = [ox oy oz]（1-based）
%
% 输出:
%   outFile - 标准化后的 4D NIfTI 文件路径（前缀 'w'）

fprintf('[normalize_apply] 读取功能像: %s\n', inFile);
[data, hdr] = nifti_read(inFile);
[nx_in, ny_in, nz_in, nt] = size(data);

fprintf('[normalize_apply] 读取位移场: %s\n', flowFile);
[flowData, flowHdr] = nifti_read(flowFile);
% flowData: [tnx tny tnz 3]，定义在模板坐标空间，各分量为模板体素位移
% flowHdr.affine 即模板的仿射矩阵（由 dartel_warp 写入）

% 位移场尺寸（模板空间）
[nx_f, ny_f, nz_f, nd] = size(flowData);
if nd ~= 3
    error('[normalize_apply] 位移场第4维应为3（x/y/z方向）');
end
if ~isfield(flowHdr, 'affine') || any(~isfinite(flowHdr.affine(:)))
    error(['[normalize_apply] 位移场缺少有效 affine，无法建立模板坐标系。', ...
           '请重新运行 dartel_warp.m 生成位移场文件: %s。', ...
           '该问题常见于旧版位移场文件或 DARTEL 过程被中断。'], flowFile);
end
if isfield(flowHdr, 'nx') && isfield(flowHdr, 'ny') && isfield(flowHdr, 'nz')
    if any([flowHdr.nx, flowHdr.ny, flowHdr.nz] ~= [nx_f, ny_f, nz_f])
        warning('[normalize_apply] flow header 维度与数据维度不一致，使用数据维度继续: hdr=[%d %d %d], data=[%d %d %d]', ...
            flowHdr.nx, flowHdr.ny, flowHdr.nz, nx_f, ny_f, nz_f);
    end
end

% -------- MNI 目标网格参数 --------
[mni_dims, mni_vox, mni_origin] = parse_target_grid(mniCfg);

mx = mni_dims(1); my = mni_dims(2); mz = mni_dims(3);

% MNI 体素 → 世界坐标（mm）
% MNI 仿射：以 origin 为零点
mni_affine = [-mni_vox(1) 0 0  mni_vox(1)*(mni_origin(1)-1);
              0  mni_vox(2) 0 -mni_vox(2)*(mni_origin(2)-1);
              0 0  mni_vox(3) -mni_vox(3)*(mni_origin(3)-1);
              0 0 0           1                             ];

% 输入图像仿射（个体空间世界坐标 → 体素坐标，用于最终采样步骤）
src_affine_inv = inv(hdr.affine);

% 从位移场头信息读取模板仿射矩阵
% dartel_warp 已将模板仿射写入位移场头，确保坐标系一致
tmpl_affine     = flowHdr.affine;
tmpl_affine_inv = inv(tmpl_affine);

% -------- 构建 MNI 网格 --------
[Xm, Ym, Zm] = ndgrid(1:mx, 1:my, 1:mz);
nMNI = mx * my * mz;

% MNI 体素（1-based）→ MNI 世界坐标（0-based 体素索引送入仿射）
mni_vox_0 = [Xm(:)'-1; Ym(:)'-1; Zm(:)'-1; ones(1, nMNI)];
mni_world  = mni_affine * mni_vox_0;  % [4 × nMNI]

% -------- 坐标映射链：MNI → 模板 → 个体 --------
% 步骤1：MNI 世界坐标 → 模板体素坐标（1-based）
%   DARTEL 模板在 MNI 空间对齐，故世界坐标可直接经模板仿射逆变换得到模板体素坐标
tmpl_vox_0 = tmpl_affine_inv * mni_world;  % [4 × nMNI]，0-based
tmpl_vox_1 = tmpl_vox_0(1:3,:) + 1;        % 0-based → 1-based

% 步骤2：在模板体素坐标处采样位移场 D
%   D(t) 含义：对模板位置 t，应在模板空间中的 t+D(t) 处采样个体图像
%   即 gm_in_tmpl(t + D(t)) ≈ template_gm(t)（DARTEL 优化目标）
dx = reshape(trilinear_interp(flowData(:,:,:,1), tmpl_vox_1), mx, my, mz);
dy = reshape(trilinear_interp(flowData(:,:,:,2), tmpl_vox_1), mx, my, mz);
dz = reshape(trilinear_interp(flowData(:,:,:,3), tmpl_vox_1), mx, my, mz);

% 步骤3：位移后的模板体素坐标（在 gm_in_tmpl 中实际采样位置）
Xt_disp = reshape(tmpl_vox_1(1,:), mx, my, mz) + dx;  % 1-based
Yt_disp = reshape(tmpl_vox_1(2,:), mx, my, mz) + dy;
Zt_disp = reshape(tmpl_vox_1(3,:), mx, my, mz) + dz;

% 步骤4：模板体素坐标（0-based）→ 世界坐标 → 个体功能像体素坐标（1-based）
%   gm_in_tmpl 是个体 GM 经仿射 individual_affine→template_affine 重采样的结果，
%   所以将模板体素坐标经模板仿射转回世界坐标，再经个体仿射逆变换可得个体体素坐标
tgt_tmpl_0 = [(Xt_disp(:)-1)'; (Yt_disp(:)-1)'; (Zt_disp(:)-1)'; ones(1, nMNI)];
tgt_world   = tmpl_affine * tgt_tmpl_0;   % 世界坐标
tgt_src_0   = src_affine_inv * tgt_world; % 个体体素坐标（0-based）
sampCoords  = tgt_src_0(1:3,:) + 1;       % 0-based → 1-based，送入 trilinear_interp

% -------- 逐时间点插值 --------
fprintf('[normalize_apply] 开始对 %d 个时间点进行空间标准化...\n', nt);
dataOut = zeros(mx, my, mz, nt, 'single');

for t = 1:nt
    vol_t = double(data(:,:,:,t));
    vals  = trilinear_interp(vol_t, sampCoords);
    dataOut(:,:,:,t) = single(reshape(vals, mx, my, mz));

    if mod(t,50)==0
        fprintf('[normalize_apply]  已处理 %d/%d 时间点\n', t, nt);
    end
end

% -------- 构建输出 NIfTI 头 --------
hdrOut = nifti_default_hdr([mx my mz nt], [mni_vox(:)' hdr.pixdim(5)]);
hdrOut.affine = mni_affine;
hdrOut.srow_x = mni_affine(1,:);
hdrOut.srow_y = mni_affine(2,:);
hdrOut.srow_z = mni_affine(3,:);
hdrOut.descrip = 'NormalizedToMNI';

% -------- 写出 --------
ensure_dir(outDir);
[~, fname, ext] = fileparts(inFile);
outFile = fullfile(outDir, ['w' fname ext]);
nifti_write(outFile, dataOut, hdrOut);
fprintf('[normalize_apply] 已写出: %s  [%d %d %d %d]\n', outFile, mx, my, mz, nt);
end

function [dims, vox, origin] = parse_target_grid(cfgIn)
% 兼容两类标准化目标参数格式
if isfield(cfgIn, 'boundingBox')
    bbox = cfgIn.boundingBox;
    if ~isequal(size(bbox), [2 3])
        error('[normalize_apply] boundingBox 必须为 2x3 矩阵');
    end
    if ~isfield(cfgIn, 'voxSize')
        error('[normalize_apply] 使用 boundingBox 时必须提供 voxSize');
    end
    vox = cfgIn.voxSize(:)';
    if numel(vox) ~= 3 || any(vox <= 0)
        error('[normalize_apply] voxSize 必须为 1x3 正数向量');
    end
    bboxMin = bbox(1,:);
    bboxMax = bbox(2,:);
    dims = round((bboxMax - bboxMin) ./ vox) + 1;
    % 将 mm 坐标的 bbox 最小角转换为 1-based 体素原点：
    % 令体素坐标 i=1 时对应 bboxMin，可得
    % bboxMin_i = vox_i * (1 - origin_i)  =>  origin_i = 1 - bboxMin_i/vox_i
    origin = 1 - bboxMin ./ vox;
else
    if ~isfield(cfgIn, 'dims') || ~isfield(cfgIn, 'voxSize') || ~isfield(cfgIn, 'origin')
        error('[normalize_apply] 目标网格配置缺少字段，需提供 boundingBox+voxSize 或 dims+voxSize+origin');
    end
    dims = cfgIn.dims(:)';
    vox = cfgIn.voxSize(:)';
    origin = cfgIn.origin(:)';
end
end
