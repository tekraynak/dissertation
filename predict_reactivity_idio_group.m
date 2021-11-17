function out = predict_reactivity_idio_group(physio, method)
if nargin< 1, physio = 'sbp'; end
if nargin < 2, method = 'lassopcr'; end
    
%% load data
cd /home/tek31/ProjectDrive/Users/Thomas/dissertation
dat = readtable(sprintf('datasets/analytic_sample_reactivity_%s.csv', physio), 'TreatAsEmpty', {'NA'});
ids = dat{:, 1};


%% load all idiographic maps
clear wdat
for i = 1:length(ids)
    load(sprintf('predict/idiographic/%s/%s/%d_%s_%s.mat', method, physio, ids(i), method, physio));
    wdat(:, i) = out.weight_obj.dat;
end

%% 10-fold cv, stratify on mean reactivity
k = 10;
CV = stratified_holdout_set(mean(dat{:, 2:end}, 2), k);


rho = nan(size(ids)); rho_pval = rho; mae = rho; bf10 = rho; r2 = rho;

%% run cv
for kk = 1:k
    % calculate mean weightmap for participants in training set
    trainwdat = mean(wdat(:, CV.trIdx{kk}), 2);
    
    testids = ids(CV.teIdx{kk});
    % load beta images for each person in test set
    for i = 1:length(testids)
        
        %% find their index in the overall ID list
        whID = find(ids == testids(i));
        
        %% name images
        for b = 1:8
            img(b) = filenames(sprintf('betaseries/%d/msit/beta_000%d.nii', testids(i), b), 'absolute');
            img(b+8) = filenames(sprintf('betaseries/%d/stroop/beta_000%d.nii', testids(i), b), 'absolute');
        end
        
        %% load images
        f = fmri_data(img, 'masks/grey.nii');
        
        %% calculate predicted reactivity
        yfit = f.dat'*trainwdat;
        
        %% get reactivity
        Y = dat{whID, 2:end}';
        
        %% calculate predicted-observed Spearman rho
        [rho(whID), rho_pval(whID)] = corr(Y, yfit, 'type', 'Spearman');

        %% calculate MAE
        mae(whID) = mean(abs(Y - yfit));

        %% calculate Bayes Factor
        % https://klabhub.github.io/bayesFactor/
        bf10(whID) = bf.regression(Y, yfit);

        %% calculate R2 - sum of squares
        Y = scale(Y); yfit = scale(yfit);
        rss = sum((Y - yfit).^2);
        tss = sum((Y - mean(Y)).^2);
        r2(whID) = 1 - (rss/tss);
        
    end
    
    fprintf(1, '\nDone with fold %s', kk)
end
    
%bf10 = filloutliers(bf10, 'nearest');
bf01 = 1./bf10;
t = table(ids, rho, rho_pval, r2, mae, bf10, bf01);
barplot_columns(t(:, 2:end))
writetable(t, sprintf('predict/idiographic/%s_%s_group-pred_stats.csv', method, physio));
    

