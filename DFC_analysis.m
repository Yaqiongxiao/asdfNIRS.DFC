%%
%%%%%%%%%%%%%%%%%%%%%%%% DFC indexes %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear;
clc;
close all;
datapath = 'xxxx';
para = load(fullfile(datapath, 'cluster_stats.mat'));
vec_dFC = zeros(44, 20);
for id = 1:44
    vec_dFC(id, 1) = para.F{1, id}(1, 1);  % PO-state1
    vec_dFC(id, 2) = para.F{1, id}(1, 2);  % PO-state2
    vec_dFC(id, 3) = para.F{1, id}(1, 3);  % PO-state3
    vec_dFC(id, 4) = para.F{1, id}(1, 4);  % PO-state4
    
    vec_dFC(id, 5) = para.MDT{1, id}(1, 1); % MDT-state1
    vec_dFC(id, 6) = para.MDT{1, id}(1, 2); % MDT-state2
    vec_dFC(id, 7) = para.MDT{1, id}(1, 3); % MDT-state3
    vec_dFC(id, 8) = para.MDT{1, id}(1, 4); % MDT-state4
    
    vec_dFC(id, 9) = para.TM{1, id}(1, 2); % state1 -> state2 prob
    vec_dFC(id, 10) = para.TM{1, id}(1, 3); % state1 -> state3 prob
    vec_dFC(id, 11) = para.TM{1, id}(1, 4); % state1 -> state4 prob
    vec_dFC(id, 12) = para.TM{1, id}(2, 1); % state2 -> state1 prob
    vec_dFC(id, 13) = para.TM{1, id}(2, 3); % state2 -> state3 prob
    vec_dFC(id, 14) = para.TM{1, id}(2, 4); % state2 -> state4 prob
    vec_dFC(id, 15) = para.TM{1, id}(3, 1); % state3 -> state1 prob
    vec_dFC(id, 16) = para.TM{1, id}(3, 2); % state3 -> state2 prob
    vec_dFC(id, 17) = para.TM{1, id}(3, 4); % state3 -> state4 prob
    vec_dFC(id, 18) = para.TM{1, id}(4, 1); % state4 -> state1 prob
    vec_dFC(id, 19) = para.TM{1, id}(4, 2); % state4 -> state2 prob
    vec_dFC(id, 20) = para.TM{1, id}(4, 3); % state4 -> state3 prob   
end

% %%%%%%%%% NOTES:
% The brain states were reordered based on their occurrence frequency (from highest to lowest). 
% The mapping between original and new state labels for the main analysis
% (window length = 60s) is as follows:
% Original State 3 (33.29%) → State 1 (most frequent)
% Original State 4 (28.61%) → State 2
% Original State 2 (25.45%) → State 3
% Original State 1 (12.66%) → State 4 (least frequent)
%%%%%%%%%

%% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%% group_infor %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Group assignment correction: SubjectASD002 should be reassigned to TD023
vec_dFC(45,:)= vec_dFC(2,:);
vec_dFC(2,:)=[];
% For window length = 60 and step size = 1 analysis, one outlier was removed
vec_dFC(16,:)=[];
asd_group = vec_dFC(1:20, :);  % 20 ASD
td_group = vec_dFC(21:43, :);  % 23 TD

asd_mean = mean(asd_group);
asd_std = std(asd_group);
asd_sem = std(asd_group) ./ sqrt(20);
td_mean = mean(td_group);
td_std = std(td_group);
td_sem = std(td_group) ./ sqrt(23);

%%%%%%%%%%%%%%%%%%%% Mann-Whitney U test %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
p_values = zeros(1, size(vec_dFC, 2));
U_value= zeros(1, size(vec_dFC, 2));
z_value= zeros(1, size(vec_dFC, 2));
for i = 1:size(vec_dFC, 2)
    [p_values(i), ~, stats] = ranksum(asd_group(:, i), td_group(:, i));  
    U_value(i) = stats.ranksum;  
    z_value(i) = stats.zval;
end

% NOTE: these state labels refer to the original labeling scheme
% State labels reported in the main text were reorganized as follows:
% Original State 3 → State 1
% Original State 4 → State 2
% Original State 2 → State 3
% Original State 1 → State 4
varNames = {'Variable', 'ASD_Mean', 'ASD_SD', 'TD_Mean', 'TD_SD', 'Z_value', 'P_value'};
charNames = {
    'State1_PO', 'State2_PO', 'State3_PO', 'State4_PO',...
    'State1_MDT','State2_MDT','State3_MDT','State4_MDT',...
    'State1 -> State2', 'State1 -> State3', 'State1 -> State4', ...
    'State2 -> State1', 'State2 -> State3', 'State2 -> State4',...
    'State3 -> State1', 'State3 -> State2', 'State3 -> State4',...
    'State4 -> State1', 'State4 -> State2', 'State4 -> State3'};
numVars = size(vec_dFC, 2);
tableData = cell(numVars, length(varNames));
for i = 1:numVars
    tableData{i,1} = charNames{i};
    tableData{i,2} = sprintf('%.3f', asd_mean(i));
    tableData{i,3} = sprintf('%.3f', asd_std(i));
    tableData{i,4} = sprintf('%.3f', td_mean(i));
    tableData{i,5} = sprintf('%.3f', td_std(i));
    tableData{i,6} = sprintf('%.3f', z_value(i));
    tableData{i,7} = sprintf('%.3f', p_values(i));
end
resultTable = cell2table(tableData, 'VariableNames', varNames);
disp(resultTable)


%%
%%%%%%%%%%%%%%%%%%%%%%%%%%% plot bar %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% set color
asd_color = [1, 0.4, 0.4]; % warm for ASD
td_color = [0.4, 0.6, 1]; % cold for TD
bar_width = 0.4;
offset = 0.2;

% plot PO
figure('Name', 'PO', 'NumberTitle', 'off');
hold on;
x_positions_asd = (1:4) - offset;
x_positions_td = (1:4) + offset;
bar(x_positions_asd, asd_mean(1:4), bar_width, 'FaceColor', asd_color, 'EdgeColor', 'none', 'DisplayName', 'ASD');
bar(x_positions_td, td_mean(1:4), bar_width, 'FaceColor', td_color, 'EdgeColor', 'none', 'DisplayName', 'TD');
errorbar(x_positions_asd, asd_mean(1:4), asd_sem(1:4), 'k.', 'LineWidth', 1.5); 
errorbar(x_positions_td, td_mean(1:4), td_sem(1:4), 'k.', 'LineWidth', 1.5);
xticks(1:4);
xticklabels({'o_State 1', 'o_State2', 'o_State3', 'o_State4'});
ylabel('PO','fontsize',16,'FontName','Arial');
set(gca, 'FontSize', 14,'FontName','Arial');
grid off;
hold off;

% plot MDT
figure('Name', 'MDT', 'NumberTitle', 'off');
hold on;
x_positions_asd = (1:4) - offset;
x_positions_td = (1:4) + offset;
bar(x_positions_asd, asd_mean(5:8), bar_width, 'FaceColor', asd_color, 'EdgeColor', 'none', 'DisplayName', 'ASD');
bar(x_positions_td, td_mean(5:8), bar_width, 'FaceColor', td_color, 'EdgeColor', 'none', 'DisplayName', 'TD');
errorbar(x_positions_asd, asd_mean(5:8), asd_sem(5:8), 'k.', 'LineWidth', 1.5);
errorbar(x_positions_td, td_mean(5:8), td_sem(5:8), 'k.', 'LineWidth', 1.5);
xticks(1:4);
xticklabels({'o_State 1', 'o_State2', 'o_State3', 'o_State4'});
set(gca, 'FontSize', 14,'FontName','Arial');
ylabel('MDT','fontsize',16,'FontName','Arial');
grid off;
hold off;

% plot state transition
figure('Name', 'Transition probability', 'NumberTitle', 'off');
hold on;
x_positions_asd = (1:3) - offset;
x_positions_td = (1:3) + offset;
bar(x_positions_asd, asd_mean(10:12), bar_width, 'FaceColor', asd_color, 'EdgeColor', 'none', 'DisplayName', 'ASD');
bar(x_positions_td, td_mean(10:12), bar_width, 'FaceColor', td_color, 'EdgeColor', 'none', 'DisplayName', 'TD');
errorbar(x_positions_asd, asd_mean(10:12), asd_sem(10:12), 'k.', 'LineWidth', 1.5);
errorbar(x_positions_td, td_mean(10:12), td_sem(10:12), 'k.', 'LineWidth', 1.5);
xticks(1:3);
xticklabels({'State4 -> State1', 'State4 -> State2', 'State4 -> State3'});
set(gca, 'FontSize', 14,'FontName','Arial');
xlim([0.5 3.5]);
grid off;
hold off;

%%     
%%%%%%%%%%%%%%%%%%%% Spearman's partial correlation  %%%%%%%%%%%%%%%%%%%%%
clear;
clc;
close all;
datapath = 'C:\Users\such\Documents\2024-rs-fNIRS\2024_final_results\data_for_github\w60_s1';
para = load(fullfile(datapath, 'cluster_stats.mat'));
sig_dFC = zeros(44, 2);
for id = 1:44
    sig_dFC(id, 1) = para.MDT{1, id}(1, 1); % State4 in the main text
    sig_dFC(id, 2) = para.TM{1, id}(1, 3); % state4 -> state1 
end
% Revision of group info and outlier removal
sig_dFC(45,:)= sig_dFC(2,:);
sig_dFC(2,:)=[];
sig_dFC(16,:)=[];

% loading covar and interested var
cova_path = 'C:\Users\such\Documents\2024-rs-fNIRS\2024_final_results\data_for_github\matlab_code';
load(fullfile(cova_path, 'cova_for_spr.mat')); % age and gender
load(fullfile(cova_path,'interested_var_for_spr.mat')); % ABC, G-adapt, task ACC and RT
r_values = NaN(2,4);
p_values = NaN(2,4);

%  Convert 'NAN' strings to numeric NaN values
for col = 1:size(var_data, 2)
    for row = 1:size(var_data, 1)
        if ischar(var_data{row, col}) && strcmp(var_data{row, col}, 'NAN')
            var_data{row, col} = NaN; 
        end
    end
end

sig_dFC_labels = {'State 4 MDT','State 4 -> State 1'};
var_name =  {'ABC','Gesell adaptive score','ACC','RT'};
% set scatter color for two groups
asd_color = [1, 0.4, 0.4]; % warm for ASD
td_color = [0.4, 0.6, 1]; % cold for TD

for i = 1:2
    for j = 1:4
        sig_var = sig_dFC(:, i);
        interest_var = cell2mat(var_data(:, j));
        valid_idx = ~isnan(interest_var); % filter for valid_data
        valid_sig_var = sig_var(valid_idx);
        valid_interest_var = interest_var(valid_idx);
        valid_cova = cova(valid_idx, :);

        [rho, p] = partialcorr(valid_sig_var, valid_interest_var, valid_cova, 'Type', 'Spearman');
        % saving the results
        r_values(i, j) = rho;
        p_values(i, j) = p;
        
        % Calculate residuals
        mdl_y = fitlm(valid_cova, valid_sig_var);
        resid_y = mdl_y.Residuals.Raw;
        mdl_x = fitlm(valid_cova, valid_interest_var);
        resid_x = mdl_x.Residuals.Raw;
        % get group id
        valid_subjects = find(valid_idx);
        asd_idx = valid_subjects <= 20; 
        td_idx = valid_subjects > 20;
    
        % linear model
        mdl = fitlm(resid_x, resid_y);
        x_fit = linspace(min(resid_x), max(resid_x), 100)';
        [y_fit, CI] = predict(mdl, x_fit, 'Alpha', 0.05, 'Prediction', 'curve');
        
        figure;
        scatter(resid_x(asd_idx), resid_y(asd_idx), 70, asd_color, 'filled', 'DisplayName', 'ASD');
        hold on;
        scatter(resid_x(td_idx), resid_y(td_idx), 70, td_color, 'filled','DisplayName', 'TD');
        plot(x_fit, y_fit, 'k-', 'LineWidth', 1.5, 'DisplayName', 'Regression Line');
        fill([x_fit; flipud(x_fit)],[CI(:,1); flipud(CI(:,2))],[0.8, 0.8, 0.8], 'FaceAlpha', 0.3, ...
            'EdgeColor', 'none', 'DisplayName', '95% CI');
        
        % Set axis
        xlabel([var_name{j}], 'FontSize', 18, 'FontName', 'Arial','FontWeight', 'bold');
        ylabel([sig_dFC_labels{i}], 'FontSize', 18, 'FontName', 'Arial','FontWeight', 'bold',...
            'Units', 'normalized', 'Position', [-0.11, 0.5, 0]);
        title(sprintf('rho = %.2f, p = %.3f', rho, p), 'FontSize', 16, 'FontName', 'Arial');
        set(gca, 'Box', 'off');
        set(gca, 'FontName', 'Arial', 'FontSize', 16);
        set(gca, 'TickLength', [0.02 0.02]); % 设置刻度线长度
        set(gca, 'LineWidth', 1.5); % 设置坐标轴线宽
        grid off;
        hold off;
    end
end