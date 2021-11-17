function out = predict_reactivity_idiographic(id, physio, method, bycondition)
if nargin < 4, bycondition = 0; end
if nargin < 3, method = 'lassopcr'; end % lassopcr using my code (not predict)
if nargin < 2, physio = 'sbp'; end
if nargin < 1, id = 1003; end

%%
wdir = cdtodrive;
dissdir = [wdir '/Users/Thomas/dissertation'];
cd(dissdir)

%% load reactivity
% updated to allow for predicting behavior (accuracy, RT)
if any(strcmp(physio, {'sbp' 'hr'}))
    dat = readtable(sprintf('datasets/analytic_sample_reactivity_%s.csv', physio), 'TreatAsEmpty', {'NA'});
    dat = dat(dat{:, 1} == id, :);
    react = dat{:, 2:end};
    % baseline = dat.([physio '_baseline']);
    % react = table2array(dat(:, 6:21)) - baseline';
elseif any(strcmp(physio, {'acc' 'rt'}))
    load('behavior/all_behav_stats.mat');
    tmp = [behav(find(id == ids), :).blockstats];
    react = horzcat(tmp.(char(['mean_' physio])))';
end


%% check for images
if ~exist(sprintf('betaseries/%d/msit', id)) || ~exist(sprintf('betaseries/%d/stroop', id))
    out.pred_obs_rho = NaN; 
    out.numblocks = NaN;
    return
end

%% name images
for i = 1:8
    img(i) = filenames(sprintf('betaseries/%d/msit/beta_000%d.nii', id, i), 'absolute');
    img(i+8) = filenames(sprintf('betaseries/%d/stroop/beta_000%d.nii', id, i), 'absolute');
end

%% remove NaN data
img = img(~isnan(react))';
react = react(~isnan(react))';

%% load images
f = fmri_data(img, 'masks/grey.nii');
f.Y = react;

%% set up folds - either LOOCV or by condition
if bycondition
    whfolds = repmat([1;2], 8, 1);
else
    whfolds = length(react); %LOOCV
end

%% lasso-pcr 
% shown are the pred-outcome-r for 1002 sbp
switch method
    case 'lassopcr'
        out = lassopcr_cv(f, [], 'nfolds', whfolds, 'rho'); % .6024
    case 'lassopcr-alpha05'
        out = lassopcr_cv(f, [], 'nfolds', whfolds, 'rho', 'alpha', .5); % .6024
    case 'predict-pcr'
        [~, out] = predict(f, 'algorithm_name', 'cv_pcr', 'nfolds', whfolds); % .6033
    case 'predict-lassopcr'
        [~, out] = predict(f, 'algorithm_name', 'cv_lassopcr', 'nfolds', whfolds, 'estimateparams');
        %[~, out] = predict(f, 'algorithm_name', 'cv_lassopcr', 'nfolds', length(react)); % .6033 - this produces the same results as PCR
    case 'predict-elasticnet' % think this works now after changing line 1618
        [~, out] = predict(f, 'algorithm_name', 'cv_lassopcrmatlab', 'Alpha', .5, 'nfolds', whfolds, 'EstimateParams', 'Lambda'); % .6033
end


% 1002 sbp yields the following pred-outcome-r
% [cverr, stats, optout] = predict(f, 'algorithm_name', 'cv_lassopcr', 'nfolds', 16, 'error_type', 'mse', 'estimateparams')
% .5613
% [~, out] = predict(f, 'algorithm_name', 'cv_lassopcrmatlab', 'Alpha', .5, 'nfolds', length(react), 'EstimateParams', 'Lambda')
% .6556
% [~, out] = predict(f, 'algorithm_name', 'cv_lassopcrmatlab', 'Alpha', 1, 'nfolds', length(react), 'EstimateParams', 'Lambda')
% .6556

%% calculate predicted-observed Spearman rho
out.rho = corr(out.Y, out.yfit, 'type', 'Spearman');

%% calculate MAE
out.mae = mean(abs(out.Y - out.yfit));

%% calculate Bayes Factor
% https://klabhub.github.io/bayesFactor/
out.bf = bf.regression(out.Y, out.yfit);

%% calculate R2 - sum of squares
rss = sum((out.Y - out.yfit).^2);
tss = sum((out.Y - mean(out.Y)).^2);
out.r2ss = 1 - (rss/tss);


%% save number of blocks
out.numblocks = length(react);

%% generate file name 
if bycondition
    fname = sprintf('predict/idiographic_bycondition/%s/%s/%d_%s_%s', method, physio, id, method, physio);
else
    fname = sprintf('predict/idiographic/%s/%s/%d_%s_%s', method, physio, id, method, physio);
end

%% save output
save(sprintf('%s.mat', fname), 'out');
if strcmp(method, 'original')
    saveas(gcf, sprintf('%s.png', fname));
end
