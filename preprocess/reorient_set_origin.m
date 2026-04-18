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
%                若为空 []，则弹出三正交视图供用户点选（无GUI时回退图像中心）
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
    acVoxCoord = select_ac_orthviews(data(:,:,:,1), defaultAC, inFile);
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

function acVoxCoord = select_ac_orthviews(vol, defaultAC, inFile)
% 三正交视图点击选点（sagittal/coronal/axial），用于 AC 选择
[nx, ny, nz] = size(vol);
acVoxCoord = defaultAC;

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
    'Units', 'normalized', 'Position', [0.1 0.15 0.8 0.7]);

axSag = subplot(1,3,1,'Parent',fig);
axCor = subplot(1,3,2,'Parent',fig);
axAxi = subplot(1,3,3,'Parent',fig);

infoText = uicontrol(fig, 'Style', 'text', 'Units', 'normalized', ...
    'Position', [0.02 0.01 0.62 0.07], 'HorizontalAlignment', 'left', ...
    'BackgroundColor', 'w', 'FontSize', 10, ...
    'String', '点击任一视图更新 AC；完成后点击“确认AC”');

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
        sag = squeeze(vol(acVoxCoord(1),:,:))';
        hSag = imagesc(axSag, sag);
        axis(axSag, 'image'); axis(axSag, 'ij'); axis(axSag, 'off');
        colormap(axSag, gray);
        hold(axSag, 'on');
        plot(axSag, [acVoxCoord(2) acVoxCoord(2)], [1 nz], 'r-', 'LineWidth', 1);
        plot(axSag, [1 ny], [acVoxCoord(3) acVoxCoord(3)], 'r-', 'LineWidth', 1);
        title(axSag, sprintf('Sagittal (X=%d) 点击更新Y/Z', acVoxCoord(1)));
        set(hSag, 'ButtonDownFcn', @onClickSag);
        set(axSag, 'ButtonDownFcn', @onClickSag);

        % Coronal (Y 固定，显示 X-Z)
        cla(axCor);
        cor = squeeze(vol(:,acVoxCoord(2),:))';
        hCor = imagesc(axCor, cor);
        axis(axCor, 'image'); axis(axCor, 'ij'); axis(axCor, 'off');
        colormap(axCor, gray);
        hold(axCor, 'on');
        plot(axCor, [acVoxCoord(1) acVoxCoord(1)], [1 nz], 'r-', 'LineWidth', 1);
        plot(axCor, [1 nx], [acVoxCoord(3) acVoxCoord(3)], 'r-', 'LineWidth', 1);
        title(axCor, sprintf('Coronal (Y=%d) 点击更新X/Z', acVoxCoord(2)));
        set(hCor, 'ButtonDownFcn', @onClickCor);
        set(axCor, 'ButtonDownFcn', @onClickCor);

        % Axial (Z 固定，显示 X-Y)
        cla(axAxi);
        axi = squeeze(vol(:,:,acVoxCoord(3)))';
        hAxi = imagesc(axAxi, axi);
        axis(axAxi, 'image'); axis(axAxi, 'ij'); axis(axAxi, 'off');
        colormap(axAxi, gray);
        hold(axAxi, 'on');
        plot(axAxi, [acVoxCoord(1) acVoxCoord(1)], [1 ny], 'r-', 'LineWidth', 1);
        plot(axAxi, [1 nx], [acVoxCoord(2) acVoxCoord(2)], 'r-', 'LineWidth', 1);
        title(axAxi, sprintf('Axial (Z=%d) 点击更新X/Y', acVoxCoord(3)));
        set(hAxi, 'ButtonDownFcn', @onClickAxi);
        set(axAxi, 'ButtonDownFcn', @onClickAxi);

        set(infoText, 'String', sprintf('当前 AC (1-based voxel): [%d %d %d]', ...
            acVoxCoord(1), acVoxCoord(2), acVoxCoord(3)));
        drawnow;
    end

    function onClickSag(~, ~)
        cp = get(axSag, 'CurrentPoint');
        acVoxCoord(2) = clamp(round(cp(1,1)), 1, ny);
        acVoxCoord(3) = clamp(round(cp(1,2)), 1, nz);
        refresh_views();
    end

    function onClickCor(~, ~)
        cp = get(axCor, 'CurrentPoint');
        acVoxCoord(1) = clamp(round(cp(1,1)), 1, nx);
        acVoxCoord(3) = clamp(round(cp(1,2)), 1, nz);
        refresh_views();
    end

    function onClickAxi(~, ~)
        cp = get(axAxi, 'CurrentPoint');
        acVoxCoord(1) = clamp(round(cp(1,1)), 1, nx);
        acVoxCoord(2) = clamp(round(cp(1,2)), 1, ny);
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

function y = clamp(x, lo, hi)
y = min(max(x, lo), hi);
end
