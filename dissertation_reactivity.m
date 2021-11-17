function out = dissertation_reactivity(physio)
if nargin < 1, physio = 'sbp'; end

%% load reactivity spreadsheet
cd /home/tek31/ProjectDrive/Users/Thomas/dissertation/datasets/
dat = readtable('analytic_sample_raw.csv', 'TreatAsEmpty', 'NA');

%% isolate to reactivity variable of interest
ids = dat.id;
react = dat(:, contains(dat.Properties.VariableNames, physio));

%% calculate baseline as the mean of the last 3 readings
clear pretask
for i = 1:height(react)
    tmp = rmmissing(table2array(react(i, 1:5)));
    pretask(i, :) = tmp(end-2:end);
end
baseline = mean(pretask, 2);
pretask_change = pretask - baseline;
    

%% create change scores
taskchange = table2array(react(:, 6:end)) - baseline;

%% find within-participant outliers (>3SD away from mean)
[a, b] = find(abs(scale(taskchange')) > 3);
for i = 1:length(a)
    taskchange(b(i), :)
end
length(a)
% 2 SBP readings need rescaled
% 5 HR readings need rescaled


%% find between-participant outliers (not as imporant here)
[a, b] = find(abs(scale(taskchange)) > 3);
for i = 1:length(a)
    taskchange(b(i), :)
end

%% use filloutliers to interpolate within-participant outliers
taskchange_rescaled = filloutliers(taskchange, 'nearest', 'mean', 2);

%% get mean and SD
subdat = [pretask_change taskchange_rescaled]';
submean = mean(subdat');
substd = std(subdat');

%% plot 3 pretask and 16 task readings as change scores - 


%% show all timeseries
figure; 
plot(1:3, subdat(1:3, :), 'Color', [.5 .5 .5])
hold on; plot(4:11, subdat(4:11, :), 'Color', [.5 .5 .5])
plot(12:19,subdat(12:19, :), 'Color', [.5 .5 .5])


%% plot means and SD's
figure; set(gcf, 'Position', [ 200 200 700 350])
line(1:3, submean(1:3), 'Color', [0 0 0], 'LineWidth', 3)
hold on; line(4:11, submean(4:11), 'Color', [0 0 0], 'LineWidth', 3)
line(12:19, submean(12:19), 'Color', [0 0 0], 'LineWidth', 3)
%errorbar(1:3, submean(1:3), substd(1:3), 'LineStyle', 'none', 'Color', [0 0 0])
errorbar(1:19, submean, substd, substd, 'LineStyle', 'none', 'Color', [0 0 0], 'LineWidth', 3)
scatter(1:3, submean(1:3),  200,  'o', 'MarkerFaceColor', [0 0 0], 'CData', [0 0 0])
scatter(4:2:19, submean(4:2:end),  200,  'o', 'MarkerFaceColor', 'r', 'CData', [0 0 0])
scatter(5:2:19, submean(5:2:end),  200,  'o', 'MarkerFaceColor', 'b', 'CData', [0 0 0])
g= gca; set(g, 'XTick', 0:20, 'XTickLabel', [], 'FontSize', 14)


%%


%% another way to plot
redblue = repmat([rgb('red') ;rgb('blue')], 4, 1);
figure; 
l = lineplot_columns([pretask_change taskchange_rescaled])
delete(l.line_han)
line(1:3, l.m(1:3), 'Color', [0 0 0], 'LineWidth', 3, 'Marker', 'o', 'MarkerSize', 8, 'MarkerFaceColor', [0.5000 0.5000 0.5000])
scatter(4:12, l.m(4:12), 'Color', [0 0 0], 'LineWidth', 3, 'Marker', 'o', 'MarkerSize', 8, 'MarkerFaceColor', redblue)


% change error bar

%% write file
csvwrite(sprintf('analytic_sample_reactivity_%s.csv', physio'), [ids taskchange_rescaled])


