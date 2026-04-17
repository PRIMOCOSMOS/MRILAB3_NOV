function outFile = remove_dummy_tr(inFile, outDir, nDummy)
% remove_dummy_tr - 去除 4D fMRI 数据中起始不稳定的 TR（Dummy Scans）
%
% 物理背景:
%   MRI 序列开始时，纵向磁化矢量 Mz 尚未达到稳态（Steady State），
%   前几个 TR 的信号偏高且不稳定，需在分析前去除。
%
% 实现原理:
%   仅在时间维进行矩阵截取: data_out = data(:,:,:, nDummy+1:end)
%
% 输入:
%   inFile  - 输入 4D NIfTI 文件路径（.nii）
%   outDir  - 输出目录
%   nDummy  - 去除的 TR 数量（整数，如 10）
%
% 输出:
%   outFile - 输出 NIfTI 文件路径（前缀 'a' 表示 after dummy removal）

fprintf('[remove_dummy_tr] 读取: %s\n', inFile);
[data, hdr] = nifti_read(inFile);

[nx, ny, nz, nt] = size(data);
fprintf('[remove_dummy_tr] 原始维度: [%d %d %d %d]，去除前 %d TR\n', ...
    nx, ny, nz, nt, nDummy);

if nt <= nDummy
    error('[remove_dummy_tr] 时间点数 %d ≤ 去除量 %d，请检查参数', nt, nDummy);
end

% -------- 时间维截取 --------
dataOut = data(:,:,:, nDummy+1:end);
ntNew   = size(dataOut, 4);
fprintf('[remove_dummy_tr] 截取后时间点数: %d\n', ntNew);

% -------- 更新头信息 --------
hdr.nt      = ntNew;
hdr.dim     = int16([4, nx, ny, nz, ntNew, 1, 1, 1]);
hdr.descrip = sprintf('DummyTR_removed=%d', nDummy);

% -------- 写出 --------
ensure_dir(outDir);
[~, fname, ext] = fileparts(inFile);
outFile = fullfile(outDir, ['a' fname ext]);
nifti_write(outFile, single(dataOut), hdr);
fprintf('[remove_dummy_tr] 已写出: %s\n', outFile);
end
