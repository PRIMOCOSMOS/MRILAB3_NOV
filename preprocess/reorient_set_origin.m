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
%                若为空 []，则自动估计为图像中心
%
% 输出:
%   outFile - 重定位后的 NIfTI 文件路径（前缀 'reorient_'）

fprintf('[reorient_set_origin] 读取: %s\n', inFile);
[data, hdr] = nifti_read(inFile);
[nx, ny, nz, ~] = size(data);

% -------- 确定 AC 体素坐标 --------
if isempty(acVoxCoord) || nargin < 3
    % 弹出对话框，提示用户手动输入前连合 AC 的体素坐标（参照 DPABI 逻辑）
    defaultAC = round([nx/2, ny/2, nz/2]);
    [~, fname_disp, ~] = fileparts(inFile);
    prompt = { ...
        sprintf('图像文件: %s\n图像维度: %d × %d × %d\n\n请在图像中定位前连合（AC）并输入体素坐标（1-based）。\n若不确定请保留估计值。\n\nAC 体素坐标 X（左右，建议 %d）：', ...
            fname_disp, nx, ny, nz, defaultAC(1)), ...
        sprintf('AC 体素坐标 Y（前后，建议 %d）：', defaultAC(2)), ...
        sprintf('AC 体素坐标 Z（上下，建议 %d）：', defaultAC(3)) };
    dlgtitle = sprintf('重定位 AC 坐标输入 — %s', fname_disp);
    dims = [1 70];
    defInput = { num2str(defaultAC(1)), num2str(defaultAC(2)), num2str(defaultAC(3)) };
    answer = inputdlg(prompt, dlgtitle, dims, defInput);

    if isempty(answer)
        % 用户取消对话框：回退到图像中心
        acVoxCoord = defaultAC;
        fprintf('[reorient_set_origin] 用户取消输入，使用图像中心: [%d %d %d]\n', ...
            acVoxCoord(1), acVoxCoord(2), acVoxCoord(3));
    else
        acVoxCoord = [str2double(answer{1}), str2double(answer{2}), str2double(answer{3})];
        % 校验输入是否为数值
        if any(isnan(acVoxCoord))
            warning('[reorient_set_origin] 输入包含非数值字符，回退到图像中心 [%d %d %d]', ...
                defaultAC(1), defaultAC(2), defaultAC(3));
            acVoxCoord = defaultAC;
        elseif acVoxCoord(1)<1 || acVoxCoord(1)>nx || ...
                acVoxCoord(2)<1 || acVoxCoord(2)>ny || acVoxCoord(3)<1 || acVoxCoord(3)>nz
            warning('[reorient_set_origin] 输入坐标 [%d %d %d] 超出图像范围 [1..%d, 1..%d, 1..%d]，回退到图像中心', ...
                acVoxCoord(1), acVoxCoord(2), acVoxCoord(3), nx, ny, nz);
            acVoxCoord = defaultAC;
        end
        fprintf('[reorient_set_origin] 用户指定 AC: [%d %d %d]\n', ...
            acVoxCoord(1), acVoxCoord(2), acVoxCoord(3));
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
