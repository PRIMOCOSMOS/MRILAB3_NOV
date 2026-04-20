function outFile = reorient_set_origin(inFile, outDir, acVoxCoord)
% reorient_set_origin - 通过修改 NIfTI 仿射矩阵将坐标原点平移到前连合（AC）
%
% 物理背景:
%   MNI 标准空间的坐标原点位于前连合（Anterior Commissure, AC）。
%   重定位操作仅修改 NIfTI 的仿射矩阵（srow_x/y/z），
%   不改变图像体素值，也不改变图像内容。
%
% 操作原理:
%   新仿射 = 旧仿射 * T_shift
%   其中 T_shift 将 AC 的体素坐标平移到世界坐标 [0 0 0]
%
% 输入:
%   inFile     - 输入 NIfTI 文件路径
%   outDir     - 输出目录
%   acVoxCoord - AC 在图像中的体素坐标（1-based），[x y z] 向量
%                例如 [32 53 28]（64×64×36 图像中心附近）
%                若为空 []，则弹出三正交视图供用户点选（无 GUI 时回退图像中心）
%
% 输出:
%   outFile - 重定位后的 NIfTI 文件路径（前缀 'reorient_'）

fprintf('[reorient_set_origin] 读取: %s\n', inFile);
[data, hdr] = nifti_read(inFile);
[nx, ny, nz, ~] = size(data);

% -------- 确定 AC 体素坐标 --------
if nargin < 3 || isempty(acVoxCoord)
    % 参照 SPM/DPABI 交互逻辑：三正交视图点击选择 AC（而非纯文本输入）
    defaultAC = round([nx/2, ny/2, nz/2]);
    acVoxCoord = select_ac_orthviews(data(:,:,:,1), defaultAC, inFile, hdr.affine);
    fprintf('[reorient_set_origin] 交互选择 AC: [%d %d %d]\n', ...
        acVoxCoord(1), acVoxCoord(2), acVoxCoord(3));
else
    % 外部传入 AC 坐标时也做边界校验
    acVoxCoord = round(acVoxCoord(:))';
    if numel(acVoxCoord) ~= 3 || ...
       acVoxCoord(1)<1 || acVoxCoord(1)>nx || ...
       acVoxCoord(2)<1 || acVoxCoord(2)>ny || ...
       acVoxCoord(3)<1 || acVoxCoord(3)>nz
        defaultAC = round([nx/2, ny/2, nz/2]);
        warning('[reorient_set_origin] 传入 AC 坐标非法，回退图像中心 [%d %d %d]', ...
            defaultAC(1), defaultAC(2), defaultAC(3));
        acVoxCoord = defaultAC;
    end
end
ac = acVoxCoord(:);

% -------- 计算 AC 在世界坐标中的位置 --------
% NIfTI 1-based 体素坐标转为 0-based 再乘仿射
ac0 = [ac - 1; 1];  % 0-based 齐次坐标
acWorld = hdr.affine * ac0;  % 世界坐标 [x y z 1]

fprintf('[reorient_set_origin] AC 世界坐标: [%.1f %.1f %.1f] mm\n', ...
    acWorld(1), acWorld(2), acWorld(3));

% -------- 修改仿射矩阵的平移分量（最后一列）--------
% 使 AC 体素对应世界坐标 [0 0 0]
% new_affine * [ac_vox_0based; 1] = [0 0 0 1]
% => new_T = -R * ac_0based  (其中 R 为仿射的旋转/缩放部分)
R = hdr.affine(1:3, 1:3);
newT = -R * (ac - 1);  % 注意 NIfTI 仿射是 0-based 的

newAffine = hdr.affine;
newAffine(1:3, 4) = newT;

fprintf('[reorient_set_origin] 仿射平移修正: [%.2f %.2f %.2f] → [%.2f %.2f %.2f]\n', ...
    hdr.affine(1,4), hdr.affine(2,4), hdr.affine(3,4), ...
    newT(1), newT(2), newT(3));

% -------- 更新头信息（不改变体素值）--------
hdr.affine = newAffine;
hdr.srow_x = newAffine(1,:);
hdr.srow_y = newAffine(2,:);
hdr.srow_z = newAffine(3,:);
hdr.qform_code = 0;
hdr.sform_code = 1;
hdr.descrip = sprintf('Reoriented AC=[%d %d %d]', ac(1), ac(2), ac(3));

% -------- 写出（数据不变）--------
ensure_dir(outDir);
[~, fname, ext] = fileparts(inFile);
outFile = fullfile(outDir, ['reorient_' fname ext]);
nifti_write(outFile, single(data), hdr);
fprintf('[reorient_set_origin] 已写出: %s\n', outFile);
end

function acVoxCoord = select_ac_orthviews(vol, defaultAC, inFile, affine)
% 三正交视图点击选点（sagittal/coronal/axial），用于 AC 选择
[nx, ny, nz] = size(vol);
acVoxCoord = defaultAC;
orient = compute_view_orientation(affine);

% 无图形环境时回退
if ~(usejava('jvm') && feature('ShowFigureWindows'))
    warning('[reorient_set_origin] 无图形界面可用，回退图像中心 AC=[%d %d %d]', ...
        defaultAC(1), defaultAC(2), defaultAC(3));
    return;
end

[~, fname, ~] = fileparts(inFile);
fig = figure('Name', sprintf('选择 AC 点（%s）', fname), ...
    'NumberTitle', 'off', 'Color', 'w', ...
    'MenuBar', 'none', 'ToolBar', 'none', ...
    'WindowScrollWheelFcn', @onScroll, ...
    'Units', 'normalized', 'Position', [0.1 0.15 0.8 0.7]);

axSag = subplot(1,3,1,'Parent',fig);
axCor = subplot(1,3,2,'Parent',fig);
axAxi = subplot(1,3,3,'Parent',fig);

infoText = uicontrol(fig, 'Style', 'text', 'Units', 'normalized', ...
    'Position', [0.02 0.01 0.62 0.07], 'HorizontalAlignment', 'left', ...
    'BackgroundColor', 'w', 'FontSize', 10, ...
    'String', '点击任一视图更新 AC；鼠标滚轮可切换当前视图切片；完成后点击“确认AC”');

uicontrol(fig, 'Style', 'pushbutton', 'String', '确认AC', ...
    'Units', 'normalized', 'Position', [0.68 0.02 0.10 0.06], ...
    'Callback', @onConfirm);
uicontrol(fig, 'Style', 'pushbutton', 'String', '使用中心', ...
    'Units', 'normalized', 'Position', [0.80 0.02 0.08 0.06], ...
    'Callback', @onCenter);
uicontrol(fig, 'Style', 'pushbutton', 'String', '取消', ...
    'Units', 'normalized', 'Position', [0.90 0.02 0.08 0.06], ...
    'Callback', @onCancel);

% 显示与交互
refresh_views();
uiwait(fig);
if isgraphics(fig), delete(fig); end

    function refresh_views()
        % Sagittal (X 固定，显示 Y-Z)
        cla(axSag);
        sag = extract_slice_2d(vol, acVoxCoord, 'sag');
        sag = apply_view_flip(sag, orient.sag.flipLR, orient.sag.flipUD);
        hSag = imagesc(axSag, sag);
        axis(axSag, 'image'); set(axSag, 'YDir', 'normal'); axis(axSag, 'off');
        colormap(axSag, gray);
        hold(axSag, 'on');
        yDisp = voxel_to_display(acVoxCoord(2), ny, orient.sag.flipLR);
        zDisp = voxel_to_display(acVoxCoord(3), nz, orient.sag.flipUD);
        plot(axSag, [yDisp yDisp], [1 nz], 'r-', 'LineWidth', 1);
        plot(axSag, [1 ny], [zDisp zDisp], 'r-', 'LineWidth', 1);
        title(axSag, sprintf('Sagittal (X=%d) 点击更新Y/Z, 滚轮切X', acVoxCoord(1)));
        set(hSag, 'ButtonDownFcn', @onClickSag);
        set(axSag, 'ButtonDownFcn', @onClickSag);

        % Coronal (Y 固定，显示 X-Z)
        cla(axCor);
        cor = extract_slice_2d(vol, acVoxCoord, 'cor');
        cor = apply_view_flip(cor, orient.cor.flipLR, orient.cor.flipUD);
        hCor = imagesc(axCor, cor);
        axis(axCor, 'image'); set(axCor, 'YDir', 'normal'); axis(axCor, 'off');
        colormap(axCor, gray);
        hold(axCor, 'on');
        xDisp = voxel_to_display(acVoxCoord(1), nx, orient.cor.flipLR);
        zDisp = voxel_to_display(acVoxCoord(3), nz, orient.cor.flipUD);
        plot(axCor, [xDisp xDisp], [1 nz], 'r-', 'LineWidth', 1);
        plot(axCor, [1 nx], [zDisp zDisp], 'r-', 'LineWidth', 1);
        title(axCor, sprintf('Coronal (Y=%d) 点击更新X/Z, 滚轮切Y', acVoxCoord(2)));
        set(hCor, 'ButtonDownFcn', @onClickCor);
        set(axCor, 'ButtonDownFcn', @onClickCor);

        % Axial (Z 固定，显示 X-Y)
        cla(axAxi);
        axi = extract_slice_2d(vol, acVoxCoord, 'axi');
        axi = apply_view_flip(axi, orient.axi.flipLR, orient.axi.flipUD);
        hAxi = imagesc(axAxi, axi);
        axis(axAxi, 'image'); set(axAxi, 'YDir', 'normal'); axis(axAxi, 'off');
        colormap(axAxi, gray);
        hold(axAxi, 'on');
        xDisp = voxel_to_display(acVoxCoord(1), nx, orient.axi.flipLR);
        yDisp = voxel_to_display(acVoxCoord(2), ny, orient.axi.flipUD);
        plot(axAxi, [xDisp xDisp], [1 ny], 'r-', 'LineWidth', 1);
        plot(axAxi, [1 nx], [yDisp yDisp], 'r-', 'LineWidth', 1);
        title(axAxi, sprintf('Axial (Z=%d) 点击更新X/Y, 滚轮切Z', acVoxCoord(3)));
        set(hAxi, 'ButtonDownFcn', @onClickAxi);
        set(axAxi, 'ButtonDownFcn', @onClickAxi);

        set(infoText, 'String', sprintf('当前 AC (1-based voxel): [%d %d %d]', ...
            acVoxCoord(1), acVoxCoord(2), acVoxCoord(3)));
        drawnow;
    end

    function onClickSag(~, ~)
        cp = get(axSag, 'CurrentPoint');
        yDisp = clamp(round(cp(1,1)), 1, ny);
        zDisp = clamp(round(cp(1,2)), 1, nz);
        acVoxCoord(2) = display_to_voxel(yDisp, ny, orient.sag.flipLR);
        acVoxCoord(3) = display_to_voxel(zDisp, nz, orient.sag.flipUD);
        refresh_views();
    end

    function onClickCor(~, ~)
        cp = get(axCor, 'CurrentPoint');
        xDisp = clamp(round(cp(1,1)), 1, nx);
        zDisp = clamp(round(cp(1,2)), 1, nz);
        acVoxCoord(1) = display_to_voxel(xDisp, nx, orient.cor.flipLR);
        acVoxCoord(3) = display_to_voxel(zDisp, nz, orient.cor.flipUD);
        refresh_views();
    end

    function onClickAxi(~, ~)
        cp = get(axAxi, 'CurrentPoint');
        xDisp = clamp(round(cp(1,1)), 1, nx);
        yDisp = clamp(round(cp(1,2)), 1, ny);
        acVoxCoord(1) = display_to_voxel(xDisp, nx, orient.axi.flipLR);
        acVoxCoord(2) = display_to_voxel(yDisp, ny, orient.axi.flipUD);
        refresh_views();
    end

    function onScroll(~, evt)
        % 鼠标滚轮：在当前悬停视图上切换该视图对应方向的切片
        obj = hittest(fig);
        ax = ancestor(obj, 'axes');
        if isempty(ax) || ~isgraphics(ax)
            return;
        end
        step = sign(evt.VerticalScrollCount);
        if step == 0
            return;
        end

        if isequal(ax, axSag)
            acVoxCoord(1) = clamp(acVoxCoord(1) - step, 1, nx);
        elseif isequal(ax, axCor)
            acVoxCoord(2) = clamp(acVoxCoord(2) - step, 1, ny);
        elseif isequal(ax, axAxi)
            acVoxCoord(3) = clamp(acVoxCoord(3) - step, 1, nz);
        end
        refresh_views();
    end

    function onConfirm(~, ~)
        uiresume(fig);
    end

    function onCenter(~, ~)
        acVoxCoord = defaultAC;
        uiresume(fig);
    end

    function onCancel(~, ~)
        acVoxCoord = defaultAC;
        uiresume(fig);
    end
end

function out = apply_view_flip(in, flipLR, flipUD)
out = in;
if flipLR
    out = fliplr(out);
end
if flipUD
    out = flipud(out);
end
end

function idxDisp = voxel_to_display(idxVox, n, isFlipped)
if isFlipped
    idxDisp = n - idxVox + 1;
else
    idxDisp = idxVox;
end
end

function idxVox = display_to_voxel(idxDisp, n, isFlipped)
if isFlipped
    idxVox = n - idxDisp + 1;
else
    idxVox = idxDisp;
end
end

function orient = compute_view_orientation(affine)
% 根据仿射矩阵推断各视图是否需要左右/上下翻转，减少“图像颠倒”
R = affine(1:3,1:3);
ex = [1; 0; 0];
ey = [0; 1; 0];
ez = [0; 0; 1];

vx = R(:,1);  % voxel-x 在世界坐标中的方向
vy = R(:,2);  % voxel-y 在世界坐标中的方向
vz = R(:,3);  % voxel-z 在世界坐标中的方向

orient.sag.flipLR = dot(vy, ey) < 0;  % Sag 水平轴对应 voxel-y
orient.sag.flipUD = dot(vz, ez) < 0;  % Sag 垂直轴对应 voxel-z

orient.cor.flipLR = dot(vx, ex) < 0;  % Cor 水平轴对应 voxel-x
orient.cor.flipUD = dot(vz, ez) < 0;  % Cor 垂直轴对应 voxel-z

orient.axi.flipLR = dot(vx, ex) < 0;  % Axi 水平轴对应 voxel-x
orient.axi.flipUD = dot(vy, ey) < 0;  % Axi 垂直轴对应 voxel-y
end

function y = clamp(x, lo, hi)
y = min(max(x, lo), hi);
end

function I = extract_slice_2d(vol, acVoxCoord, mode)
switch lower(mode)
    case 'sag' % X 固定，显示 Y-Z
        I = squeeze(vol(acVoxCoord(1),:,:))';
    case 'cor' % Y 固定，显示 X-Z
        I = squeeze(vol(:,acVoxCoord(2),:))';
    case 'axi' % Z 固定，显示 X-Y
        I = squeeze(vol(:,:,acVoxCoord(3)))';
    otherwise
        error('extract_slice_2d: 未知 mode=%s', mode);
end
end
