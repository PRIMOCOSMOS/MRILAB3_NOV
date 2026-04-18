function Vo = trilinear_interp(V, coords)
% trilinear_interp - 对3D体数据进行三线性插值
% 完全 standalone，不依赖任何工具箱
%
% 输入:
%   V      - 3D 数组 [nx ny nz]，参考体数据
%   coords - [3 × N] 矩阵，每列为一个插值点的 (x,y,z) 坐标（1-based 体素坐标）
%
% 输出:
%   Vo     - [1 × N] 插值结果，超出边界的点赋值 0（NaN 可选）
%
% 使用示例:
%   V = rand(64,64,36);
%   pts = [32.3; 32.1; 18.7];   % 单点
%   val = trilinear_interp(V, pts);

[nx, ny, nz] = size(V);
if size(coords,1) ~= 3
    error('trilinear_interp: coords 必须为 [3 x N]');
end

N = size(coords, 2);
if N == 0
    Vo = zeros(1,0);
    return;
end

% 统一使用列向量，避免 MATLAB 隐式扩展导致 N×N 内存爆炸
Vo = zeros(N, 1);
V_flat = double(V(:));
% 分块大小经验值：在常见 fMRI 体数据下显著降低峰值内存且保持速度
blockSize = 2e5;

maxX = max(nx - 1, 1);
maxY = max(ny - 1, 1);
maxZ = max(nz - 1, 1);

for i0 = 1:blockSize:N
    i1 = min(i0 + blockSize - 1, N);
    r = i0:i1;

    x = double(coords(1, r)).';
    y = double(coords(2, r)).';
    z = double(coords(3, r)).';

    % 判断边界内的点
    inBounds = (x >= 1) & (x <= nx) & ...
               (y >= 1) & (y <= ny) & ...
               (z >= 1) & (z <= nz);

    % 裁剪坐标到有效范围（边界外点先设为1，最后清零）
    x = max(min(x, maxX), 1);
    y = max(min(y, maxY), 1);
    z = max(min(z, maxZ), 1);

    % 取整和小数部分
    x0 = floor(x); x1 = min(x0 + 1, nx);
    y0 = floor(y); y1 = min(y0 + 1, ny);
    z0 = floor(z); z1 = min(z0 + 1, nz);

    dx = x - x0;
    dy = y - y0;
    dz = z - z0;

    % 将 (ix, iy, iz) 转换为线性索引
    idx000 = sub2ind([nx ny nz], x0, y0, z0);
    idx100 = sub2ind([nx ny nz], x1, y0, z0);
    idx010 = sub2ind([nx ny nz], x0, y1, z0);
    idx110 = sub2ind([nx ny nz], x1, y1, z0);
    idx001 = sub2ind([nx ny nz], x0, y0, z1);
    idx101 = sub2ind([nx ny nz], x1, y0, z1);
    idx011 = sub2ind([nx ny nz], x0, y1, z1);
    idx111 = sub2ind([nx ny nz], x1, y1, z1);

    % 8个角点的值
    c000 = V_flat(idx000);
    c100 = V_flat(idx100);
    c010 = V_flat(idx010);
    c110 = V_flat(idx110);
    c001 = V_flat(idx001);
    c101 = V_flat(idx101);
    c011 = V_flat(idx011);
    c111 = V_flat(idx111);

    % 三线性插值公式
    vals = c000 .* (1-dx) .* (1-dy) .* (1-dz) + ...
           c100 .*    dx  .* (1-dy) .* (1-dz) + ...
           c010 .* (1-dx) .*    dy  .* (1-dz) + ...
           c110 .*    dx  .*    dy  .* (1-dz) + ...
           c001 .* (1-dx) .* (1-dy) .*    dz  + ...
           c101 .*    dx  .* (1-dy) .*    dz  + ...
           c011 .* (1-dx) .*    dy  .*    dz  + ...
           c111 .*    dx  .*    dy  .*    dz;

    % 边界外置零
    vals(~inBounds) = 0;
    Vo(r) = vals;
end

Vo = Vo.';
end
