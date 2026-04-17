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
N = size(coords, 2);

x = coords(1,:);
y = coords(2,:);
z = coords(3,:);

% 判断边界内的点
inBounds = (x >= 1) & (x <= nx) & ...
           (y >= 1) & (y <= ny) & ...
           (z >= 1) & (z <= nz);

% 裁剪坐标到有效范围（边界外点先设为1，最后清零）
x = max(min(x, nx-1), 1);
y = max(min(y, ny-1), 1);
z = max(min(z, nz-1), 1);

% 取整和小数部分
x0 = floor(x); x1 = x0 + 1;
y0 = floor(y); y1 = y0 + 1;
z0 = floor(z); z1 = z0 + 1;

dx = x - x0;
dy = y - y0;
dz = z - z0;

% 防止越界（x1, y1, z1 不超过最大值）
x1 = min(x1, nx);
y1 = min(y1, ny);
z1 = min(z1, nz);

% 展平V以便快速索引
V_flat = V(:);

% 将 (ix, iy, iz) 转换为线性索引
idx = @(ix, iy, iz) sub2ind([nx ny nz], ix, iy, iz);

% 8个角点的值
c000 = V_flat(idx(x0, y0, z0));
c100 = V_flat(idx(x1, y0, z0));
c010 = V_flat(idx(x0, y1, z0));
c110 = V_flat(idx(x1, y1, z0));
c001 = V_flat(idx(x0, y0, z1));
c101 = V_flat(idx(x1, y0, z1));
c011 = V_flat(idx(x0, y1, z1));
c111 = V_flat(idx(x1, y1, z1));

% 三线性插值公式
Vo = c000 .* (1-dx) .* (1-dy) .* (1-dz) + ...
     c100 .*    dx  .* (1-dy) .* (1-dz) + ...
     c010 .* (1-dx) .*    dy  .* (1-dz) + ...
     c110 .*    dx  .*    dy  .* (1-dz) + ...
     c001 .* (1-dx) .* (1-dy) .*    dz  + ...
     c101 .*    dx  .* (1-dy) .*    dz  + ...
     c011 .* (1-dx) .*    dy  .*    dz  + ...
     c111 .*    dx  .*    dy  .*    dz;

% 边界外置零
Vo(~inBounds) = 0;
end
