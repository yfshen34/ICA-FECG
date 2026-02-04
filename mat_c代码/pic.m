offset=10.0;
qrsAf=qrsAf(qrsAf >= offset);
qrsF_all_global=qrsF_all+ offset;
% qrsF_all_global=qrsF_all+ offset;
%验证matlab内的online前10s数据丢失

% offset=10.0;
% qrsAf=qrsAf(qrsAf <= (60-offset));
% % qrsF_all_global=qrsF_all(qrsF_all <= (60-offset));
% 
% qrsF_all_global=qrsF_all;

% offset=2.0;
% qrsAf=qrsAf(qrsAf <= (60-offset));
% qrsF_all_global=qrsF_all(qrsF_all <= (60-offset));



%qrsF_all_global=[];


%验证matlab内的尾部2s数据丢失


% plotSignalComparisonClose(qrsAf,qrsF_all_global); % 使用之前定义的绘图函数

tolerance=0.05;
recall = calculateRecallRate(qrsAf, qrsF_all_global, tolerance);
precision = calculatePrecision(qrsAf, qrsF_all_global, tolerance);
plotSignalComparisonClose(qrsAf, qrsF_all_global, precision, recall);

 text(4.2, 1.9, sprintf('召回率: %.1f%%', recall), ...
        'FontSize', 12, 'FontWeight', 'bold', ...
        'BackgroundColor', [1 1 1 0.7], 'EdgeColor', 'k');
 hold off;

function plotSignalComparisonClose(PFP_times, APFP_times, precision, recall)
    % 创建图形窗口
    figure('Color', 'white', 'Position', [100, 100, 800, 500]);
    
    % 设置轴范围和标签
    xlim([0 65]);
    ylim([0 2]);
    xlabel('时间 (秒)');
    title('信号位置对比图');
    set(gca, 'YTick', [0.5 1.2], 'YTickLabel', {'qrsF all', 'qrsAf'});
    grid on;
    hold on;
    
    % 绘制PFP数据 (蓝色竖线，在上部)
    for i = 1:length(PFP_times)
        x = PFP_times(i);
        line([x x], [1.55 0.85], 'Color', [0 0.4470 0.7410], 'LineWidth', 1.5);
    end
    
    % 绘制APFP数据 (橙色竖线，在下部)
    for i = 1:length(APFP_times)
        x = APFP_times(i);
        line([x x], [0.15 0.85], 'Color', [0.8500 0.3250 0.0980], 'LineWidth', 1.5);
    end
    
    % 添加黑色标记（示例位置）
    plot(0.2, 1.2, 'k_', 'MarkerSize', 10, 'LineWidth', 2); % PFP标记
    plot(0.2, 0.5, 'k_', 'MarkerSize', 10, 'LineWidth', 2); % APFP标记
    
    % 添加参考线
    plot([0 5], [1.2 1.2], 'k--', 'LineWidth', 0.5); % PFP参考线
    plot([0 5], [0.5 0.5], 'k--', 'LineWidth', 0.5); % APFP参考线
    
    % ===== 新增精确率标注 =====
    % 假设 precision 和 recall 已经在调用前计算好
    text(4.2, 1.7, sprintf('精确率: %.1f%%', precision), ...
        'FontSize', 12, 'FontWeight', 'bold', ...
        'BackgroundColor', [1 1 1 0.7], 'EdgeColor', 'k');
    
    % 原有的召回率标注（保持不变）
    text(4.2, 1.9, sprintf('召回率: %.1f%%', recall), ...
        'FontSize', 12, 'FontWeight', 'bold', ...
        'BackgroundColor', [1 1 1 0.7], 'EdgeColor', 'k');
    
    hold off;
end



% function plotSignalComparisonClose(PFP_times, APFP_times)
%     % 创建图形窗口
%     figure('Color', 'white', 'Position', [100, 100, 800, 500]);
% 
%     % 设置轴范围和标签
%     xlim([0 65]);
%     ylim([0 2]);
%     xlabel('时间 (秒)');
%     title('信号位置对比图');
%     set(gca, 'YTick', [0.5 1.2], 'YTickLabel', {'qrsF all', 'qrsAf'});
%     grid on;
%     hold on;
% 
%     % 绘制PFP数据 (蓝色竖线，在上部)
%     for i = 1:length(PFP_times)
%         x = PFP_times(i);
%         line([x x], [1.55 0.85], 'Color', [0 0.4470 0.7410], 'LineWidth', 1.5);
%     end
% 
%     % 绘制APFP数据 (橙色竖线，在下部)
%     for i = 1:length(APFP_times)
%         x = APFP_times(i);
%         line([x x], [0.15 0.85], 'Color', [0.8500 0.3250 0.0980], 'LineWidth', 1.5);
%     end
% 
%     % 添加黑色标记（示例位置）
%     plot(0.2, 1.2, 'k_', 'MarkerSize', 10, 'LineWidth', 2); % PFP标记
%     plot(0.2, 0.5, 'k_', 'MarkerSize', 10, 'LineWidth', 2); % APFP标记
% 
%     % 添加参考线
%     plot([0 5], [1.2 1.2], 'k--', 'LineWidth', 0.5); % PFP参考线
%     plot([0 5], [0.5 0.5], 'k--', 'LineWidth', 0.5); % APFP参考线
% 
% end

function recall = calculateRecallRate(ref_times, detected_times, tolerance)

    % 初始化匹配计数器
    true_positives = 0;
    
    % 遍历所有参考信号
    for i = 1:length(ref_times)
        ref_time = ref_times(i);
        
        % 检查是否有检测信号在容差范围内
        diffs = abs(detected_times - ref_time);
        if any(diffs <= tolerance)
            true_positives = true_positives + 1;
        end
    end
    
    % 计算召回率
    recall = (true_positives / length(ref_times)) * 100;
    
    % 显示结果
    fprintf('召回率计算结果:\n');
    fprintf('参考信号总数: %d\n', length(ref_times));
    fprintf('正确识别的信号数: %d\n', true_positives);
    fprintf('召回率: %.2f%%\n', recall);
end

function precision = calculatePrecision(ref_times, detected_times, tolerance)
    % 初始化计数器
    true_positives = 0;
    false_positives = 0;
    
    % 遍历所有检测信号
    for i = 1:length(detected_times)
        detected_time = detected_times(i);
        
        % 检查是否有参考信号在容差范围内
        diffs = abs(ref_times - detected_time);
        if any(diffs <= tolerance)
            true_positives = true_positives + 1;
        else
            false_positives = false_positives + 1;
        end
    end
    
    % 计算精确率
    precision = (true_positives / (true_positives + false_positives)) * 100;

    
    
    % 显示结果
    fprintf('精确率计算结果:\n');
    fprintf('检测信号总数: %d\n', length(detected_times));
    fprintf('正确检测的信号数: %d\n', true_positives);
    fprintf('误检信号数: %d\n', false_positives);
    fprintf('精确率: %.2f%%\n', precision);
end