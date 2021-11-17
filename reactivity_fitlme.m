cd /home/tek31/ProjectDrive/Users/Thomas/dissertation

physio = {'sbp' 'hr'}

for p = 1:2
    disp(physio{p})
    dat = readtable(sprintf('datasets/analytic_sample_reactivity_%s.csv', physio{p}));
    disp('Overall reactivity')
    [~,P,~,STATS] = ttest(mean(dat{:, 2:end}, 2))
    disp('mixed effects model')
    % mixed effects model
    t = table;
    t.phys = reshape(dat{:, 2:end}', [size(dat, 1)*(size(dat, 2)-1), 1]);
    t.id = repelem(dat{:, 1}, 16);
    t.condition = repmat([2; 1], size(dat, 1)*8, 1);
    t.task = repmat([zeros(8, 1); ones(8, 1)], size(dat, 1), 1);
    writetable(t, sprintf('datasets/analytic_sample_reactivity_%s_long.csv', physio{p}));
    fitlme(t, 'phys ~ condition * task + (1 | id)', 'FitMethod', 'REML')
    
    
end




[h, p, ci, stats] = ttest(mean(dat{:, 2:2:end}, 2), mean(dat{:, 3:2:end}, 2))