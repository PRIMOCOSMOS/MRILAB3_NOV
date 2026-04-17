function M = rigid_mat(p)
% rigid_mat - 由6个刚体参数构造 4×4 齐次仿射变换矩阵
% 参数约定与 SPM 一致: [tx ty tz rx ry rz]
%
% 输入:
%   p - [6×1] 或 [1×6] 向量
%       p(1:3) = 平移量 [tx ty tz]（mm）
%       p(4:6) = 旋转角 [rx ry rz]（弧度），绕 x/y/z 轴
%
% 输出:
%   M - [4×4] 刚体变换矩阵（右乘规范）
%       将参考坐标变换到当前坐标: x_ref = M * x_curr
%
% 使用示例:
%   M = rigid_mat([1 2 3, 0.01 0 -0.01]);  % 平移1,2,3mm，绕X轴旋转0.01rad

p = p(:)';  % 确保行向量

tx = p(1); ty = p(2); tz = p(3);
rx = p(4); ry = p(5); rz = p(6);

% 绕X轴旋转矩阵
Rx = [1,      0,       0;
      0,  cos(rx), -sin(rx);
      0,  sin(rx),  cos(rx)];

% 绕Y轴旋转矩阵
Ry = [ cos(ry), 0, sin(ry);
            0,  1,      0;
      -sin(ry), 0, cos(ry)];

% 绕Z轴旋转矩阵
Rz = [cos(rz), -sin(rz), 0;
      sin(rz),  cos(rz), 0;
           0,        0,  1];

% 总旋转：R = Rz * Ry * Rx（与SPM约定一致）
R = Rz * Ry * Rx;

% 构建 4×4 矩阵
M = eye(4);
M(1:3, 1:3) = R;
M(1:3, 4)   = [tx; ty; tz];
end
