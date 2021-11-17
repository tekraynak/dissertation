

load('/home/tek31/ProjectDrive/Users/Thomas/dissertation/behavior/all_behav_stats.mat')



%% correlations between ACC/RT and predicted reactivity

for p = 1:2
    react = readtable(sprintf('datasets/analytic_sample_reactivity_%s.csv', phys{p}))
    for i = 1:size(behav, 1)
     load(sprintf('predict/idiographic/lassopcr/%s/%d_lassopcr_%s.mat', phys{p}, ids(i), phys{p}));
     acc_rho(i, p) = corr(out.yfit(out.cv_folds.fold_indicator), acc(i, :)', 'type', 'Spearman', 'rows', 'complete');
     rt_rho(i, p) = corr(out.yfit(out.cv_folds.fold_indicator), rt(i, :)', 'type', 'Spearman', 'rows', 'complete');
    end
end

save('behavior/predicted_behavior_associations.mat', 'acc_rho', 'rt_rho')



%% correlations between ACC/RT and observed reactivity

for p = 1:2
    react = readtable(sprintf('datasets/analytic_sample_reactivity_%s.csv', phys{p}))
    for i = 1:size(behav, 1)
     load(sprintf('predict/idiographic/lassopcr/%s/%d_lassopcr_%s.mat', phys{p}, ids(i), phys{p}));
     acc_rho(i, p) = corr(react{i, 2:end}', acc(i, :)', 'type', 'Spearman', 'rows', 'complete');
     rt_rho(i, p) = corr(react{i, 2:end}', rt(i, :)', 'type', 'Spearman', 'rows', 'complete');
    end
end
     
save('behavior/observed_behavior_associations.mat', 'acc_rho', 'rt_rho')


    
%% PARTIAL CORRELATIONS



for p = 1:2
    react = readtable(sprintf('datasets/analytic_sample_reactivity_%s.csv', phys{p}))
    for i = 1:size(behav, 1)
     load(sprintf('predict/idiographic/lassopcr/%s/%d_lassopcr_%s.mat', phys{p}, ids(i), phys{p}));
     acc_rho(i, p) = partialcorr(out.yfit(out.cv_folds.fold_indicator), acc(i, :)', repmat([1;2], 8, 1), 'type', 'Spearman', 'rows', 'complete');
     rt_rho(i, p) = partialcorr(out.yfit(out.cv_folds.fold_indicator), rt(i, :)', repmat([1;2], 8, 1), 'type', 'Spearman', 'rows', 'complete');
    end
end
save('behavior/predicted_behavior_associations_partialcorr.mat', 'acc_rho', 'rt_rho')


for p = 1:2
    react = readtable(sprintf('datasets/analytic_sample_reactivity_%s.csv', phys{p}))
    for i = 1:size(behav, 1)
     load(sprintf('predict/idiographic/lassopcr/%s/%d_lassopcr_%s.mat', phys{p}, ids(i), phys{p}));
     acc_rho(i, p) = partialcorr(react{i, 2:end}', acc(i, :)', repmat([1;2], 8, 1), 'type', 'Spearman', 'rows', 'complete');
     rt_rho(i, p) = partialcorr(react{i, 2:end}', rt(i, :)', repmat([1;2], 8, 1), 'type', 'Spearman', 'rows', 'complete');
    end
end

save('behavior/observed_behavior_associations_partialcorr.mat', 'acc_rho', 'rt_rho')
