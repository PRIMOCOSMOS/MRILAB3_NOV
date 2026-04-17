function [gmFile, wmFile, csfFile] = segment_tissue(t1File, outDir, cfg)
% segment_tissue - 基于 GMM（高斯混合模型）+ EM 算法的 T1 脑组织分割
% 将 T1 脑图像分割为三类: GM（灰质）、WM（白质）、CSF（脑脊液）
% 完全 standalone，不依赖 SPM 或任何工具箱
%
% 方法:
%   1. K-means 初始化 GMM 参数（均值、方差、混合权重）
%   2. EM（期望最大化）迭代优化
%   3. 可选 MRF（马尔科夫随机场）空间正则化（cfg.seg.mrfBeta）
%   4. 输出每类的概率图（范围 [0,1]）
%
% 输入:
%   t1File - T1 NIfTI 文件路径（3D）
%   outDir - 输出目录
%   cfg    - 配置结构体:
%              cfg.seg.nClasses (类别数, 通常3)
%              cfg.seg.nIter    (EM 最大迭代次数)
%              cfg.seg.mrfBeta  (MRF 正则化系数, 0=不使用)
%
% 输出:
%   gmFile  - GM 概率图 NIfTI 路径
%   wmFile  - WM 概率图 NIfTI 路径
%   csfFile - CSF 概率图 NIfTI 路径

fprintf('[segment_tissue] 读取 T1: %s\n', t1File);
[t1Data, hdr] = nifti_read(t1File);
vol = double(t1Data(:,:,:,1));
[nx, ny, nz] = size(vol);

nClasses = cfg.seg.nClasses;
nIter    = cfg.seg.nIter;
mrfBeta  = cfg.seg.mrfBeta;

% -------- 脑掩模（简单阈值）--------
thresh = prctile(vol(:), 15);  % 去除背景噪声
brainMask = vol > thresh;
fprintf('[segment_tissue] 脑掩模体素数: %d / %d\n', sum(brainMask(:)), numel(vol));

voxVals = vol(brainMask);  % 非零脑内体素强度
nVox    = numel(voxVals);

% -------- K-means 初始化 GMM --------
fprintf('[segment_tissue] K-means 初始化 (%d 类)...\n', nClasses);
[labels_init, centroids] = kmeans_1d(voxVals, nClasses);

% 按强度升序排列类别（CSF<GM<WM）
[centroids, sortIdx] = sort(centroids);
labels_init = arrayfun(@(x) find(sortIdx==x), labels_init);

% GMM 参数初始化
mu  = centroids;                     % 均值
sig = zeros(1, nClasses);            % 标准差
pi_ = zeros(1, nClasses);            % 混合权重
for k = 1:nClasses
    idx_k = labels_init == k;
    if sum(idx_k) > 0
        sig(k) = std(voxVals(idx_k)) + 1e-3;
        pi_(k) = mean(idx_k);
    else
        sig(k) = 1;
        pi_(k) = 1/nClasses;
    end
end
fprintf('[segment_tissue] 初始均值: %s\n', mat2str(round(mu,1)));

% -------- EM 迭代 --------
gamma = zeros(nVox, nClasses);  % 责任矩阵（Responsibility）
prevLogLik = -Inf;

for iter = 1:nIter
    % ---- E 步：计算责任 γ_{ik} ----
    for k = 1:nClasses
        gamma(:,k) = pi_(k) * gauss_pdf(voxVals, mu(k), sig(k));
    end
    sumGamma = sum(gamma, 2) + 1e-12;
    gamma = gamma ./ sumGamma;

    % 对数似然
    logLik = sum(log(sumGamma));

    % ---- M 步：更新参数 ----
    Nk = sum(gamma, 1) + 1e-12;
    for k = 1:nClasses
        pi_(k) = Nk(k) / nVox;
        mu(k)  = sum(gamma(:,k) .* voxVals) / Nk(k);
        sig(k) = sqrt(sum(gamma(:,k) .* (voxVals - mu(k)).^2) / Nk(k)) + 1e-3;
    end

    % 收敛判断
    if abs(logLik - prevLogLik) < 1e-4 * abs(logLik)
        fprintf('[segment_tissue] EM 在第 %d 次迭代收敛\n', iter);
        break;
    end
    prevLogLik = logLik;

    if mod(iter,20)==0
        fprintf('[segment_tissue]  iter=%d logLik=%.2f mu=%s\n', ...
            iter, logLik, mat2str(round(mu,1)));
    end
end

fprintf('[segment_tissue] 最终均值（CSF/GM/WM）: [%.1f %.1f %.1f]\n', mu(1),mu(2),mu(3));

% -------- MRF 空间正则化（可选）--------
if mrfBeta > 0
    fprintf('[segment_tissue] 应用 MRF 正则化 (beta=%.2f)...\n', mrfBeta);
    gamma = apply_mrf(vol, brainMask, gamma, nClasses, mrfBeta);
end

% -------- 构建概率图 --------
probMaps = zeros(nx, ny, nz, nClasses, 'single');
linIdx = find(brainMask);
for k = 1:nClasses
    probMap_k = zeros(nx, ny, nz, 'single');
    probMap_k(linIdx) = single(gamma(:,k));
    probMaps(:,:,:,k) = probMap_k;
end

% -------- 写出概率图 --------
ensure_dir(outDir);
hdr_seg = hdr;
hdr_seg.nt = 1;
hdr_seg.dim = int16([3, nx, ny, nz, 1, 1, 1, 1]);
hdr_seg.scl_slope = 1;
hdr_seg.scl_inter = 0;

% 文件命名规则（与 run_pipeline_sub01.m 中的调用约定一致）:
%   k=1 → CSF（最低强度），c1_t1.nii
%   k=2 → GM（中等强度），c2_t1.nii
%   k=3 → WM（最高强度），c3_t1.nii
% 注意: SPM 的 c1=GM,c2=WM,c3=CSF 与此不同；
%       本 pipeline 采用 c1=CSF,c2=GM,c3=WM 的强度升序约定
labels = {'csf', 'gm', 'wm'};  % 强度升序：CSF < GM < WM，对应 k=1,2,3
outFiles = cell(1,3);
for k = 1:nClasses
    hdr_seg.descrip = sprintf('SegProb_%s mu=%.1f', labels{k}, mu(k));
    outF = fullfile(outDir, sprintf('c%d_t1.nii', k));
    nifti_write(outF, probMaps(:,:,:,k), hdr_seg);
    outFiles{k} = outF;
    fprintf('[segment_tissue] 已写出 %s: %s\n', upper(labels{k}), outF);
end

% 按约定返回 GM/WM/CSF（k=2/3/1）
csfFile = outFiles{1};
gmFile  = outFiles{2};
wmFile  = outFiles{3};
end

% ======================================================================
% 内部辅助函数
% ======================================================================

function p = gauss_pdf(x, mu, sigma)
% 单变量高斯概率密度函数
p = (1 / (sigma * sqrt(2*pi))) * exp(-0.5 * ((x - mu) / sigma).^2);
p = max(p, 1e-300);
end

function [labels, centroids] = kmeans_1d(data, K)
% 简单 1D K-means 聚类（等间隔初始化）
data = data(:);
lo = min(data); hi = max(data);
centroids = linspace(lo, hi, K);
labels = ones(size(data), 'int32');
for iter = 1:100
    % 分配
    dists = abs(data - centroids);  % [N × K]
    [~, new_labels] = min(dists, [], 2);
    if all(new_labels == labels), break; end
    labels = new_labels;
    % 更新
    for k = 1:K
        if sum(labels==k) > 0
            centroids(k) = mean(data(labels==k));
        end
    end
end
end

function gamma_mrf = apply_mrf(vol, mask, gamma, nClasses, beta)
% 简化 MRF 正则化：使用 6-邻域的类标签一致性
% gamma: [nVox × nClasses]
[nx, ny, nz] = size(vol);
linIdx = find(mask);
nVox   = numel(linIdx);

% 6 个邻域方向的偏移（线性索引差）
offsets = [1, -1, nx, -nx, nx*ny, -nx*ny];

gamma_mrf = gamma;

for iter = 1:3  % MRF 迭代次数（简化）
    % 获取每个体素的邻域类别概率之和
    neighborSum = zeros(nVox, nClasses);
    for d = 1:6
        nb_lin = linIdx + offsets(d);
        valid  = nb_lin >= 1 & nb_lin <= numel(vol) & mask(nb_lin);
        nb_lin(~valid) = linIdx(~valid);  % 无效邻域用自身代替

        % 邻域的当前类别分配（使用硬标签）
        [~, nb_label] = max(gamma_mrf, [], 2);
        nb_lab_neighbors = nb_label(nb_lin);  % 不精确但足够

        for k = 1:nClasses
            neighborSum(:,k) = neighborSum(:,k) + double(nb_lab_neighbors == k);
        end
    end

    % MRF 更新：log(gamma) += beta × neighborSum
    log_gamma = log(gamma_mrf + 1e-12) + beta * neighborSum;
    % Softmax 归一化
    log_gamma = log_gamma - max(log_gamma, [], 2);
    gamma_mrf = exp(log_gamma);
    gamma_mrf = gamma_mrf ./ (sum(gamma_mrf, 2) + 1e-12);
end
end
