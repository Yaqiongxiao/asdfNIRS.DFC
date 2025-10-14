%%
%%%%%%%%%%%%%%%%%%%%%%  --- PART ONE ---- %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear;
clc;
close all;
datapath = 'xxx'; % change here
para = load(fullfile(datapath, 'cluster_stats.mat'));
dirContent = dir(datapath);
folders = dirContent([dirContent.isdir] & ~ismember({dirContent.name}, {'.', '..'}));
nSubj = length(folders);
for id = 1:nSubj
    dFC(id, 1) = para.MDT{1, id}(1, 2); % MDT
end
dFC_label = {'State4 MDT'};
cova_path = 'xxx';
load(fullfile(cova_path, 'var_cova_24subj.mat'));
r_values = NaN(size(dFC,2),size(var,2));
p_values = NaN(size(dFC,2),size(var,2));
td_color = [0.4, 0.6, 1]; 
for i = 1:size(dFC,2)
    for j = 1:size(var,2)
        sig_var = dFC(:, i); 
        interest_var =var(:, j);
        valid_idx = ~isnan(interest_var); 
        valid_sig_var = sig_var(valid_idx);
        valid_interest_var = interest_var(valid_idx);
        valid_cova = cova(valid_idx, :); 
        [rho, p] = partialcorr(valid_sig_var, valid_interest_var, valid_cova, 'Type', 'Spearman');
        r_values(i, j) = rho;
        p_values(i, j) = p;
        mdl_y = fitlm(valid_cova, valid_sig_var);
        resid_y = mdl_y.Residuals.Raw;
        mdl_x = fitlm(valid_cova, valid_interest_var);
        resid_x = mdl_x.Residuals.Raw;
        mdl = fitlm(resid_x, resid_y);
        x_fit = linspace(min(resid_x), max(resid_x), 100)';
        [y_fit, CI] = predict(mdl, x_fit, 'Alpha', 0.05, 'Prediction', 'curve');
        
        figure;
        scatter(resid_x, resid_y, 70, td_color, 'filled', 'DisplayName', 'TD');
        hold on;
        plot(x_fit, y_fit, 'k-', 'LineWidth', 1.5, 'DisplayName', 'Regression Line');
        fill([x_fit; flipud(x_fit)], ...
            [CI(:,1); flipud(CI(:,2))], ...
            [0.8, 0.8, 0.8], 'FaceAlpha', 0.3, ...
            'EdgeColor', 'none', 'DisplayName', '95% CI');
        xlabel([var_name{j}], 'FontSize', 18, 'FontName', 'Arial','FontWeight', 'bold');
        ylabel([dFC_label{i}], 'FontSize', 18, 'FontName', 'Arial','FontWeight', 'bold',...
            'Units', 'normalized', 'Position', [-0.11, 0.5, 0]);
        title(sprintf('rho = %.2f, p = %.3f', rho, p), 'FontSize', 16, 'FontName', 'Arial');
        
        set(gca, 'Box', 'off');
        set(gca, 'FontName', 'Arial', 'FontSize', 16);
        set(gca, 'TickLength', [0.02 0.02]);
        set(gca, 'LineWidth', 1.5);
        grid off;
        hold off;
    end
end

%%
%%%%%%%%%%%%%%%%%%%%%%  --- PART TWO ---- %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear;
clc;
close all;
datapath = 'xxx'; % change here
load(fullfile(datapath, 'MvpaResults.mat'));
color = [0.4, 0.6, 1]; 
actual_RT = MvpaResults.predLabel(:,1);
predicted_RT = MvpaResults.predLabel(:,2);
resid_x = actual_RT;
resid_y = predicted_RT;
mdl = fitlm(resid_x, resid_y);
x_fit = linspace(min(resid_x), max(resid_x), 100)';
[y_fit, CI] = predict(mdl, x_fit, 'Alpha', 0.05, 'Prediction', 'curve');
figure;
scatter(resid_x, resid_y, 70, color, 'filled', 'DisplayName', 'TD');
hold on;
plot(x_fit, y_fit, 'k-', 'LineWidth', 1.5, 'DisplayName', 'Regression Line');
fill([x_fit; flipud(x_fit)], ...
    [CI(:,1); flipud(CI(:,2))], ...
    [0.8, 0.8, 0.8], 'FaceAlpha', 0.3, ...
    'EdgeColor', 'none', 'DisplayName', '95% CI');
xlabel('Actual RT', 'FontSize', 18, 'FontName', 'Arial', 'FontWeight', 'bold');
ylabel('Predicted RT', 'FontSize', 18, 'FontName', 'Arial', 'FontWeight', 'bold',...
    'Units', 'normalized', 'Position', [-0.11, 0.5, 0]);
[r, p] = corr(resid_x, resid_y);
title(sprintf('r = %.3f, p = %.3f', r, p), 'FontSize', 16, 'FontName', 'Arial');
set(gca, 'Box', 'off');
set(gca, 'FontName', 'Arial', 'FontSize', 16);
set(gca, 'TickLength', [0.02 0.02]); 
set(gca, 'LineWidth', 1.5);
grid off;
hold off;
% saveas(gcf, 'Actual_vs_Predicted_RT.png');