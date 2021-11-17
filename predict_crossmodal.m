function predict_crossmodal(method)
if nargin< 1, method = 'lassopcr'; end
    
%% load data
cd /home/tek31/ProjectDrive/Users/Thomas/dissertation
dat = readtable('datasets/analytic_sample_raw.csv', 'TreatAsEmpty', 'NA');
ids = dat.id;

phys = {'sbp' 'hr'};


%% predict HR from SBP and vice versa
train = [1 2]; test = [2 1];
for tt = 1:2
    
    fprintf(1, '\nUsing %s to predict %s\n', phys{train(tt)}, phys{test(tt)})
    
    rho = nan(size(ids)); rho_pval = rho; mae = rho; bf10 = rho; r2 = rho;

    for i = 1:length(ids)
         %% load both sets of predictions
        for p = 1:2
            load(sprintf('predict/idiographic/%s/%s/%d_%s_%s.mat', method, phys{p}, ids(i), method, phys{p}));
            cdat{p} = out;
        end
        yfit = double(cdat{train(tt)}.yfit);
        Y = cdat{test(tt)}.Y;
        
        %% calculate predicted-observed Spearman rho
        [rho(i), rho_pval(i)] = corr(Y, yfit, 'type', 'Spearman');

        %% calculate MAE
        mae(i) = mean(abs(Y - yfit));

        %% calculate Bayes Factor
        % https://klabhub.github.io/bayesFactor/
        bf10(i) = bf.regression(Y, yfit);

        %% calculate R2 - sum of squares
        Y = scale(Y); yfit = scale(yfit);
        rss = sum((Y - yfit).^2);
        tss = sum((Y - mean(Y)).^2);
        r2(i) = 1 - (rss/tss);
    end
    
    bf10 = filloutliers(bf10, 'nearest');
    bf01 = 1./bf10;
    t = table(ids, rho, rho_pval, r2, mae, bf10, bf01);
    barplot_columns(t(:, 2:end))
    writetable(t, sprintf('predict/idiographic/%s_crossmodal_%s-to-%s.csv', method, phys{train(tt)}, phys{test(tt)}));
    
end

        
    
        
