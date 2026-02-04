function nlfunc = autoFitNonlinearity(y)
% autoFitNonlinearityTruePDF
% 通过真实信号y推算出 p(y)，并计算得出新的非线性函数 f(y)
%
% 输入:
%   y - (nComponents x nSamples) 分离后的信号矩阵
%
% 输出:
%   nlfunc - 插值得到的新 score function，可以直接作为 dynamicOrica 的 nlfunc 输入

fprintf('开始基于数据推算真正的非线性函数...\n');

if nargin < 1
    error('必须输入分离信号矩阵 y！示例：autoFitNonlinearityTruePDF(res)');
end

% ==== 选择处理的通道 ====
channelIdx = 1; % 默认只处理第一个分离出来的分量
yy = y(channelIdx,:);

% ==== 核密度估计 p(y) ====
[f, xi] = ksdensity(yy, 'NumPoints', 200); % 用200个点来拟合更平滑

% 避免log(0)的问题
f(f <= 0) = min(f(f > 0));

% ==== 计算 log(p(y)) 和数值导数 ====
logf = log(f);
dlogf = gradient(logf, xi); % 对xi求梯度
f_y = -dlogf; % 取负号，得到 score function

% ==== 绘制估计结果 ====
figure('Name', '真实 p(y) 和推导的 f(y)');
subplot(2,1,1);
plot(xi, f, 'b-', 'LineWidth', 1.5);
xlabel('y'); ylabel('p(y)');
title('核密度估计 p(y)');
grid on;

subplot(2,1,2);
plot(xi, f_y, 'r-', 'LineWidth', 1.5);
xlabel('y'); ylabel('f(y)');
title('推算出的 score function f(y)');
grid on;

% ==== 创建插值函数 ====
fprintf('创建基于数据的插值非线性函数...\n');
nlfunc = @(yy) interp1(xi, f_y, yy, 'linear', 'extrap');

fprintf('新的非线性函数已生成！可以直接用于 dynamicOrica.\n');

end
