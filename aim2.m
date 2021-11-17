
cd /home/tek31/ProjectDrive/Users/Thomas/dissertation
dat = readtable('datasets/analytic_sample_raw.csv');
ids = dat.id;

%% set whether the similarity measure was cosine or dot
similarity_method = 'cosine'

%% apply signatures
for i = 1:length(ids)
    apply_signature(ids(i), similarity_method)
end

%% load output from above step
clear dat;
for i = 1:length(ids)
    load(sprintf('predict/signatures/%s/%d.mat', similarity_method, ids(i)));
    dat(i) = out;
end
out = dat;

%% dat.rho is a 2x2 matrix correlating the predicted and observed values
% where rows = physio and colums = pattern
% physio measures are HR & SBP
% patterns are Eisenbarth (HR) & Gianaros (SBP)
% therefore:
%   rho(1, 1) is correlation of HR and Eisenbarth map
%   rho(1, 2) is correlation of HR and Gianaros map
%   rho(2, 1) is correlation of SBP and Eisenbarth map
%   rho(2, 2) is correlation of SBP and Gianaros map
% same patterna applies for stats below
physio = {'HR' 'SBP'}; signature = {'Eisenbarth' 'Gianaros'};

for s = 1:2
    for p = 1:2
        fprintf(1, '\n Predicting %s from %s\n', physio{p}, signature{s});
        rho = nan(size(ids)); rho_pval = rho; mae = rho; bf10 = rho; r2 = rho;
        for i = 1:length(dat)
            rho(i) = dat(i).rho(p, s);
            %rho_pval(i) = dat(i).rho_pval(p, s);
            r2(i) = dat(i).r2ss(p, s);
            mae(i) = dat(i).mae(p, s);
            bf10(i) = dat(i).bf(p, s);   
        end
        bf10 = filloutliers(bf10, 'nearest');
        bf01 = 1./bf10;
        t = table(ids, rho, r2, mae, bf10, bf01);
        barplot_columns(t(:, 2:end))
        writetable(t, sprintf('predict/signatures/%s_%s_%s_stats.csv', similarity_method, physio{p}, signature{s}));
    end
end

    

%% concatenate spearman rho values
clear rho
for i = 1:length(dat)
    rho(:, :, i) = dat(i).rho;
end

%% run ttest on each signature / physiology combination
for p = 1:2
    for s = 1:2
        [~, pval(p, s), ~, stats(p, s)] = ttest(squeeze(rho(p, s, :)));        
    end
end

for p = 1:2
    for s = 1:2
        fprintf(1, '\nPrediction of %s from %s signature\n is min = %.2f max = %.2f\n\n', ...
            	physio{p}, signature{s}, min(squeeze(rho(p, s, :))), max(squeeze(rho(p, s, :))))
        barplot_columns(squeeze(rho(p, s, :)));
    end
end


