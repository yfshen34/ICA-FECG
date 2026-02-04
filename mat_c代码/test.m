% ====================== 胎动检测完整代码 ======================
clear; clc; close all; warning off;  % 清空环境并关闭警告

% ---------------------- 参数配置区 ----------------------
% 数据读取参数
file_path = 'F:\Data_2025\EEG_log_20250916_163101.csv';  % 替换为实际数据路径
targetVars = {'ax_g', 'ay_g', 'az_g', 'gx_dps', 'gy_dps', 'gz_dps'};  % 目标列名

% 特征计算参数
fs = 50;                          % IMU采样率（Hz，需与数据实际采样率一致）
win_size_motion = 0.2;            % 运动检测窗口大小（秒，适配短时间胎动）
step_motion = 0.05;               % 运动检测步长（秒，提高时间分辨率）
win_size_rest = 1.0;              % 静息基线窗口大小（秒，用于动态阈值更新）

% 多特征阈值参数
acc_threshold_factor = 2.0;       % 加速度方差阈值系数（均值+系数*标准差）
kurtosis_threshold = 4.0;         % 峰度阈值（>4视为尖锐分布，胎动特征）
zcr_threshold = 0.25;             % 过零率阈值（次/秒，高频波动特征）
freq_band = [1, 5];               % 胎动主频率带（Hz，对应胎儿踢腿/转身）

% Z轴主导验证参数
z_axis_ratio_threshold = 0.6;     % Z轴方差占比阈值（>60%视为垂直方向主导）
min_segment_variance = 0.01;      % 最小有效段方差（过滤噪声段，单位：g²）
verbose = true;                   % 是否输出调试信息（true=开启）

% ---------------------- 数据读取与预处理 ----------------------
% 读取CSV数据并转换为表格
try
    T = readtable(file_path);  % 读取表格数据
catch ME
    error('数据读取失败！请检查文件路径是否正确：\n%s\n错误信息：%s', file_path, ME.message);
end

% 匹配目标列（加速度+角速度）
varNames = T.Properties.VariableNames;  
matchedNames = {}; matchedIdx = []; 
for i = 1:length(varNames)
    if ismember(varNames{i}, targetVars)
        matchedNames{end+1} = varNames{i};
        matchedIdx(end+1) = i;
    end
end
if isempty(matchedNames)
    error('未找到目标列！请检查targetVars与CSV列名是否匹配。');
end

% 提取IMU数据并降采样（可选，根据实际需求调整）
six_axis_table = T(:, matchedIdx); 
imu_data = table2array(six_axis_table);  
imu = imu_data(1:4:end, :);  % 降采样（4倍下采样，减少计算量）
ax = imu(:,1); ay = imu(:,2); az = imu(:,3);  % 三轴加速度（g）
gx = imu(:,4); gy = imu(:,5); gz = imu(:,6);  % 三轴角速度（dps）

% 生成时间轴（基于原始采样率）
n_samples = numel(ax);
t = (0:n_samples-1)/fs;  % 时间轴（秒）

% ---------------------- 特征计算（多特征融合） ----------------------
fprintf('开始计算多特征...\n');
acc_norm = sqrt(ax.^2 + ay.^2 + az.^2);
% 初始化特征存储数组
num_win = floor((n_samples - win_size_motion*fs) / (step_motion*fs)) + 1;
acc_var = zeros(num_win,1);       % 加速度方差
kurtosis_acc = zeros(num_win,1);  % 加速度峰度
zcr_acc = zeros(num_win,1);       % 加速度过零率
freq_energy = zeros(num_win,1);   % 主频率带能量占比

 
% 滑动窗口计算特征
for idx = 1:num_win
    % 窗口索引（基于原始数据）
    start_idx = (idx-1)*step_motion*fs + 1;
    end_idx = start_idx + win_size_motion*fs - 1;
    if end_idx > n_samples
        break;  % 防止越界
    end
    % 注意：acc_norm需在循环前计算
    win_acc = acc_norm(start_idx:end_idx); 

    % 1. 加速度模长（提前计算，避免重复）
    % 注意：原代码中acc_norm未定义，此处补充计算
    acc_norm_win = sqrt(ax(start_idx:end_idx).^2 + ay(start_idx:end_idx).^2 + az(start_idx:end_idx).^2);
    
    % 2. 方差（适配短时间波动）
    acc_var(idx) = var(acc_norm_win);
    
    % 3. 峰度（捕捉尖锐冲击）
    kurtosis_acc(idx) = kurtosis(acc_norm_win);  % MATLAB kurtosis函数计算无偏峰度
    
    % 4. 过零率（高频波动次数）
    zcr_acc(idx) = sum(diff(sign(acc_norm_win - mean(acc_norm_win))) ~= 0) / (length(acc_norm_win)-1) / fs;
    
    % 5. 频域能量（主频率带占比）
    fft_vals = fft(acc_norm_win - mean(acc_norm_win));
    freqs = (0:length(fft_vals)-1)*fs/length(fft_vals);
    valid_freqs = freqs >= freq_band(1) & freqs <= freq_band(2);
    if any(valid_freqs)
        freq_energy(idx) = sum(abs(fft_vals(valid_freqs)).^2) / sum(abs(fft_vals).^2);
    else
        freq_energy(idx) = 0;  % 无有效频率成分
    end
end

% ---------------------- 动态阈值计算（基于静息段统计） ----------------------
fprintf('计算动态阈值...\n');

% 检测静息段（方差低于初始阈值的区域）
init_th_var = mean(acc_var(1:min(100,num_win))) + acc_threshold_factor*std(acc_var(1:min(100,num_win)));
rest_mask = acc_var < init_th_var;  % 静息段掩码（方差低于初始阈值）

% 分割静息段（避免过长静息段干扰）
rest_segments = detect_rest_segments(t, rest_mask, win_size_rest);  % 自定义静息段分割函数

% 基于静息段计算动态阈值（取静息段统计量的均值+标准差）
rest_var = acc_var([rest_segments.start_idx; rest_segments.end_idx]);  % 静息段方差
if ~isempty(rest_var)
    th_var = mean(rest_var) + acc_threshold_factor*std(rest_var);  % 加速度方差阈值
    th_kurtosis = mean(kurtosis_acc([rest_segments.start_idx; rest_segments.end_idx])) + (kurtosis_threshold - mean(kurtosis_acc([rest_segments.start_idx; rest_segments.end_idx])));  % 动态峰度阈值
    th_zcr = mean(zcr_acc([rest_segments.start_idx; rest_segments.end_idx])) + (zcr_threshold - mean(zcr_acc([rest_segments.start_idx; rest_segments.end_idx])));  % 动态过零率阈值
    th_freq_energy = mean(freq_energy([rest_segments.start_idx; rest_segments.end_idx])) + (mean(freq_energy([rest_segments.start_idx; rest_segments.end_idx])) - mean(freq_energy([rest_segments.start_idx; rest_segments.end_idx])));  % 动态频域能量阈值（示例，需根据实际调整）
else
    % 数据不足时使用保底阈值
    th_var = 0.1; th_kurtosis = 4.0; th_zcr = 0.25; th_freq_energy = 0.3;
    warning('静息段数据不足，使用保底阈值！');
end

% ---------------------- 多特征融合运动检测 ----------------------
fprintf('检测运动段...\n');

% 综合多特征判定运动点
is_motion = (acc_var > th_var) & ...
            (kurtosis_acc > th_kurtosis) & ...
            (zcr_acc > th_zcr) & ...
            (freq_energy > th_freq_energy);

% 合并运动段（过滤短时间噪声）
min_segment_duration=0.1
motion_segments = merge_motion_segments(t, is_motion, min_segment_duration);  % 自定义合并函数

% ---------------------- Z轴主导验证（关键后处理） ----------------------
fprintf('验证Z轴主导段...\n');

valid_segments = [];  % 初始化有效段列表

for i = 1:length(motion_segments)
    % 提取当前段的时间范围
    seg = motion_segments(i);
    seg_start = seg.start_time;
    seg_end = seg.end_time;
    
    % 调试信息输出
    if verbose
        fprintf('验证第%d段：起始=%.2fs，结束=%.2fs，原始持续时间=%.2fs\n',...
                i, seg_start, seg_end, seg_end - seg_start);
    end
    
    % 筛选当前段对应的加速度数据索引
    seg_mask = (t >= seg_start) & (t <= seg_end);
    if sum(seg_mask) == 0
        warning(sprintf('第%d段无有效数据，跳过！', i));
        continue;
    end
    
    % 提取当前段的三轴加速度数据
    ax_seg = ax(seg_mask);
    ay_seg = ay(seg_mask);
    az_seg = az(seg_mask);
    
    % 计算各轴加速度的方差
    ax_var_seg = var(ax_seg);
    ay_var_seg = var(ay_seg);
    az_var_seg = var(az_seg);
    total_var = ax_var_seg + ay_var_seg + az_var_seg;
    
    % 过滤总方差过小的段（噪声）
    if total_var < min_segment_variance
        if verbose
            fprintf('第%d段总方差=%.4f < 最小阈值=%.4f，跳过！\n',...
                    i, total_var, min_segment_variance);
        end
        continue;
    end
    
    % 计算Z轴方差占比（避免除以零）
    if total_var == 0
        warning(sprintf('第%d段总方差为零，跳过！', i));
        continue;
    end
    z_ratio = az_var_seg / total_var;
    
    % 验证Z轴主导条件
    if z_ratio > z_axis_ratio_threshold
        % 记录有效段详细信息
        valid_seg = struct('start_time', seg_start,'end_time', seg_end,'duration', seg_end - seg_start,'z_ratio', z_ratio, 'ax_var', ax_var_seg,'ay_var', ay_var_seg,'az_var', az_var_seg);
        valid_segments(end+1) = valid_seg;
        
        if verbose
            fprintf('第%d段通过验证！Z轴占比=%.1f%%\n', i, z_ratio*100);
        end
    else
        if verbose
            fprintf('第%d段未通过验证！Z轴占比=%.1f%% < 阈值=%.1f%%\n',...
                    i, z_ratio*100, z_axis_ratio_threshold*100);
        end
    end
end

% ---------------------- 结果可视化 ----------------------
if ~isempty(valid_segments)
    fprintf('检测到%d个有效胎动段！\n', length(valid_segments));
    plot_results(t, ax, acc_norm, valid_segments, th_var, th_kurtosis, th_zcr);
else
    warning('未检测到有效胎动段！');
end

% ====================== 辅助函数 ======================
% 检测静息段（方差低于阈值的连续区域）
function rest_segments = detect_rest_segments(t, mask, max_duration)
    diff_mask = diff([0; mask(:); 0]);
    rest_starts = find(diff_mask == 1);
    rest_ends = find(diff_mask == -1) - 1;
    rest_segments = struct('start_time', {}, 'end_time', {}, 'start_idx', {}, 'end_idx', {});
    
    for i = 1:length(rest_starts)
        s_time = t(rest_starts(i));
        e_time = t(rest_ends(i));
        if (e_time - s_time) <= max_duration
            rest_segments(end+1).start_time = s_time;
            rest_segments(end).end_time = e_time;
            rest_segments(end).start_idx = rest_starts(i);
            rest_segments(end).end_idx = rest_ends(i);
        end
    end
end

% 合并运动段（基于最小持续时间）
function merged_segments = merge_motion_segments(t, is_motion, min_duration)
    diff_flag = diff([0; is_motion(:); 0]);
    raw_starts = find(diff_flag == 1);
    raw_ends = find(diff_flag == -1) - 1;
    merged_segments = struct('start_time', {}, 'end_time', {});
    
    if isempty(raw_starts)
        return;
    end
    
    current_start = t(raw_starts(1));
    current_end = t(raw_ends(1));
    
    for k = 2:length(raw_starts)
        next_start = t(raw_starts(k));
        next_end = t(raw_ends(k));
        
        % 允许小间隔合并（间隔≤min_duration/2）
        if (next_start - current_end) <= (min_duration/2)
            current_end = max(current_end, next_end);
        else
            if (current_end - current_start) >= min_duration
                merged_segments(end+1).start_time = current_start;
                merged_segments(end).end_time = current_end;
            end
            current_start = next_start;
            current_end = next_end;
        end
    end
    
    % 保存最后一个段（检查持续时间）
    if (current_end - current_start) >= min_duration
        merged_segments(end+1).start_time = current_start;
        merged_segments(end).end_time = current_end;
    end
end

% 结果可视化函数
function plot_results(t, ax, acc_norm, valid_segments, th_var, th_kurtosis, th_zcr)
    figure('Position', [100 100 1200 1000], 'Name', '胎动检测结果');
    
    % 子图1：加速度模长与有效段
    subplot(4,1,1);
    plot(t, acc_norm, 'b-', 'LineWidth', 1.2); hold on;
    y_lim = ylim;
    for i = 1:length(valid_segments)
        seg = valid_segments(i);
        patch([seg.start_time seg.end_time seg.end_time seg.start_time],...
              [y_lim(1) y_lim(1) y_lim(2) y_lim(2)],...
              [0.9 0.9 0.7], 'FaceAlpha', 0.3, 'EdgeColor', 'r');
        text(seg.start_time + (seg.end_time - seg.start_time)/2, y_lim(2)*0.9,...
             sprintf('Z占比=%.1f%%', seg.z_ratio*100),...
             'HorizontalAlignment', 'center', 'Color', 'r', 'FontSize', 8);
    end
    title('加速度模长与有效胎动段（Z轴主导）');
    xlabel('时间 (s)');
    ylabel('加速度模长 (g)');
    legend('加速度模长', '有效段', 'Location', 'best');
    grid on; hold off;
    
    % 子图2：加速度方差与阈值
    subplot(4,1,2);
    win_size_plot = 0.5;  % 绘图窗口大小（秒）
    [acc_var_plot, t_plot] = moving_average(acc_var, round(win_size_plot*fs));
    plot(t_plot, acc_var_plot, 'r-', 'LineWidth', 1.2); hold on;
    yline(th_var, '--r', sprintf('Th_Var=%.3f', th_var), 'LineWidth', 1.5);
    title('加速度方差与阈值');
    xlabel('时间 (s)');
    ylabel('方差 (g²)');
    legend('方差', '阈值', 'Location', 'best');
    grid on; hold off;
    
    % 子图3：峰度与阈值
    subplot(4,1,3);
    [kurtosis_plot, t_plot] = moving_average(kurtosis_acc, round(win_size_plot*fs));
    plot(t_plot, kurtosis_plot, 'm-', 'LineWidth', 1.2); hold on;
    yline(th_kurtosis, '--m', sprintf('Th_Kurtosis=%.1f', th_kurtosis), 'LineWidth', 1.5);
    title('加速度峰度与阈值');
    xlabel('时间 (s)');
    ylabel('峰度');
    legend('峰度', '阈值', 'Location', 'best');
    grid on; hold off;
    
    % 子图4：过零率与阈值
    subplot(4,1,4);
    [zcr_plot, t_plot] = moving_average(zcr_acc, round(win_size_plot*fs));
    plot(t_plot, zcr_plot, 'c-', 'LineWidth', 1.2); hold on;
    yline(th_zcr, '--c', sprintf('Th_ZCR=%.2f', th_zcr), 'LineWidth', 1.5);
    title('过零率与阈值');
    xlabel('时间 (s)');
    ylabel('过零率 (次/秒)');
    legend('过零率', '阈值', 'Location', 'best');
    grid on; hold off;
end

% 滑动平均函数（用于平滑绘图）
function [smooth_data, smooth_t] = moving_average(data, window_size)
    if window_size < 1
        smooth_data = data;
        smooth_t = (0:length(data)-1)/fs;
        return;
    end
    half_window = floor(window_size/2);
    smooth_data = zeros(size(data));
    smooth_t = (0:length(data)-1)/fs;
    for i = 1:length(data)
        start_idx = max(1, i - half_window);
        end_idx = min(length(data), i + half_window);
        smooth_data(i) = mean(data(start_idx:end_idx));
    end
end