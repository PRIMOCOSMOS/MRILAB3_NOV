function hrf = hrf_canonical(t, params)
% hrf_canonical - 计算 SPM 双伽马函数（Double Gamma）标准 HRF
%
% SPM 标准 HRF 是两个伽马函数的差:
%   h(t) = (t/d1)^a1 * exp(-(t-d1)/b1) / (a1! * b1)
%          - c * (t/d2)^a2 * exp(-(t-d2)/b2) / (a2! * b2)
%   其中 a1=6, a2=16, b1=b2=1, c=1/6（SPM 默认值）
%
% 输入:
%   t      - 时间向量（秒），例如 0:0.1:30
%   params - （可选）HRF 参数结构体，含:
%              params.peak_delay    (峰值延迟, 默认 6s)
%              params.under_delay   (下冲延迟, 默认 16s)
%              params.peak_disp     (峰值弥散, 默认 1)
%              params.under_disp    (下冲弥散, 默认 1)
%              params.ratio         (峰/下冲比, 默认 6)
%              params.onset         (HRF 开始延迟, 默认 0)
%
% 输出:
%   hrf - HRF 时间序列（与 t 等长，峰值归一化为1）
%
% 使用示例:
%   t = 0:0.1:32;
%   h = hrf_canonical(t);
%   plot(t, h);

% -------- 默认参数（与 SPM 一致）--------
p.peak_delay  = 6;    % 峰值延迟（秒）
p.under_delay = 16;   % 下冲延迟（秒）
p.peak_disp   = 1;    % 峰值色散（b1）
p.under_disp  = 1;    % 下冲色散（b2）
p.ratio       = 6;    % 峰/下冲幅度比
p.onset       = 0;    % HRF 开始偏移

if nargin >= 2 && isstruct(params)
    fields = fieldnames(params);
    for i = 1:numel(fields)
        p.(fields{i}) = params.(fields{i});
    end
end

t = t(:)' - p.onset;

% -------- 伽马 PDF（MATLAB 内置 gampdf 可用，但此处自行实现保持 standalone）--------
% Γ(a,b): 均值=a*b，峰值在 t=(a-1)*b
a1 = p.peak_delay  / p.peak_disp;    % 形状参数
b1 = p.peak_disp;                     % 尺度参数
a2 = p.under_delay / p.under_disp;
b2 = p.under_disp;

% 伽马 PDF: f(t; a, b) = t^(a-1) * exp(-t/b) / (b^a * Gamma(a))
hrf1 = gamma_pdf(t, a1, b1);
hrf2 = gamma_pdf(t, a2, b2);

hrf = hrf1 - hrf2 / p.ratio;

% 负时间置零
hrf(t < 0) = 0;

% 峰值归一化
peak_val = max(hrf);
if peak_val > 0
    hrf = hrf / peak_val;
end
end

% ======================================================================
% 内部辅助函数
% ======================================================================

function p = gamma_pdf(x, a, b)
% 伽马分布概率密度函数（standalone，不依赖统计工具箱）
% f(x; a, b) = x^(a-1) * exp(-x/b) / (b^a * Gamma(a))
p = zeros(size(x));
idx = x > 0;
p(idx) = x(idx).^(a-1) .* exp(-x(idx)/b) / (b^a * gamma_func(a));
end

function g = gamma_func(n)
% 伽马函数 Γ(n) = (n-1)! 对正实数的推广
% 使用 Stirling 近似或精确递推
if n == round(n) && n > 0
    % 整数情况: Γ(n) = (n-1)!，使用 gammaln 提升数值稳定性
    g = exp(gammaln(n));
else
    % 使用 MATLAB 内置（仅 gamma 本身不需要工具箱）
    g = gamma(n);
end
end
