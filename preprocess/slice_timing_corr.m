function outFile = slice_timing_corr(inFile, outDir, sliceTimingMs, refSliceIdx, TR)
% slice_timing_corr - 基于傅里叶时移性质的切片时序校正（Slice Timing Correction）
%
% 物理背景:
%   EPI 序列中，同一 TR 内不同层的采集时刻不同，
%   导致一个脑体积内各层的时间参考点不统一。
%   时序校正将所有层的时间基准统一到参考层（refSlice）的采集时刻。
%
% 数学原理（傅里叶时移定理）:
%   若 f(t) 的傅里叶变换为 F(ω)，则 f(t - Δt) 的傅里叶变换为 F(ω)·e^{-iωΔt}。
%   因此，在频域对每层时间序列乘以相移因子，然后 IFFT 即可实现时间插值。
%
% 输入:
%   inFile       - 输入 4D NIfTI 文件路径
%   outDir       - 输出目录
%   sliceTimingMs- [1×nSlices] 每层相对于 TR 开始的采集时间（ms）
%   refSliceIdx  - 参考层索引（1-based），所有层校正到该层的采集时刻
%   TR           - 重复时间（秒）
%
% 输出:
%   outFile - 输出 NIfTI 文件路径（前缀 'st' 表示 slice-timing corrected）

fprintf('[slice_timing_corr] 读取: %s\n', inFile);
[data, hdr] = nifti_read(inFile);
[nx, ny, nz, nt] = size(data);

if numel(sliceTimingMs) ~= nz
    error('[slice_timing_corr] sliceTimingMs 长度 %d 与层数 %d 不匹配', ...
        numel(sliceTimingMs), nz);
end

fprintf('[slice_timing_corr] 维度=[%d %d %d %d], TR=%.2fs, 参考层=%d\n', ...
    nx, ny, nz, nt, TR, refSliceIdx);

% -------- 参考时间（ms）--------
refTime = sliceTimingMs(refSliceIdx);

% -------- FFT 频率向量 --------
% 频率 (Hz): k/nt/TR，其中 k = [0,1,...,nt/2,-nt/2+1,...,-1]
TR_ms = TR * 1000;  % 转换为 ms

% 构造频率轴（用于相移计算）
freq_hz = [0:floor((nt-1)/2), -floor(nt/2):-1] / (nt * TR);  % Hz

dataOut = single(data);

% -------- 逐层校正 --------
for z = 1:nz
    % 该层相对于参考层的时间差（秒），需要提前该量
    dt_sec = (sliceTimingMs(z) - refTime) / 1000;  % 转换为秒

    % 对该层所有体素的时间序列进行相移
    % 取出该层: [nx*ny × nt]
    slice2d = reshape(data(:,:,z,:), nx*ny, nt);  % [nVox × nt]

    % 1D FFT 沿时间维
    F = fft(slice2d, [], 2);  % [nVox × nt]

    % 相移因子: e^{-i * 2π * freq * dt}（负号表示时间超前）
    phaseShift = exp(-1i * 2 * pi * freq_hz * dt_sec);  % [1 × nt]
    phaseShift = repmat(phaseShift, nx*ny, 1);

    % 乘以相移并逆变换
    corrected = real(ifft(F .* phaseShift, [], 2));

    dataOut(:,:,z,:) = reshape(single(corrected), nx, ny, 1, nt);

    if mod(z,6)==0
        fprintf('[slice_timing_corr]  已处理 %d/%d 层 (dt=%.1fms)\n', z, nz, dt_sec*1000);
    end
end

% -------- 写出 --------
hdr.descrip = sprintf('SliceTimingCorr refSlice=%d TR=%.2f', refSliceIdx, TR);
ensure_dir(outDir);
[~, fname, ext] = fileparts(inFile);
outFile = fullfile(outDir, ['st' fname ext]);
nifti_write(outFile, dataOut, hdr);
fprintf('[slice_timing_corr] 已写出: %s\n', outFile);
end
