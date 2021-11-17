function out = predict_reactivity_idiographic_bycondition(id, physio)
if nargin < 2, physio = 'sbp'; end
if nargin < 1, id = 1002; end

%%
wdir = cdtodrive;
dissdir = [wdir '/Users/Thomas/dissertation'];
cd(dissdir)

% %% load reactivity
% dat = readtable('datasets/pip_reactivity_2020.csv', 'TreatAsEmpty', {'NA'});
% if ~any (id == dat.id); return; end
% dat = dat(dat.id == id, contains(dat.Properties.VariableNames, physio));
% 
% baseline = dat.([physio '_baseline']);
% react = table2array(dat(:, 6:21)) - baseline';

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
    tmp.pred_obs_rho = NaN; 
    tmp.numblocks = NaN;
    return
end

%% name images
for i = 1:8
    img(i) = filenames(sprintf('betaseries/%d/msit/beta_000%d.nii', id, i), 'absolute');
    img(i+8) = filenames(sprintf('betaseries/%d/stroop/beta_000%d.nii', id, i), 'absolute');
end

%% load images and run LASSO-PCR according to condition
cond = {'incongruent' 'congruent'};
clear f out tmp

for c = 1:2
    f(c) = fmri_data(img(c:2:16), 'masks/grey.nii');
    f(c).Y = react(c:2:16);

    %% remove NaN data
    f(c).dat = f(c).dat(:, ~isnan(f(c).Y));
    f(c).Y = f(c).Y(~isnan(f(c).Y));

    %% lasso-pcr on that specific condition
    tmp(c) = lassopcr_cv(f(c), [], 'nfolds', length(f(c).Y), 'rho');

    %% calculate predicted-observed Spearman rho
    rho(c) = corr(tmp(c).Y, tmp(c).yfit, 'type', 'Spearman');

    %% save number of blocks
   %out(c).numblocks = length(out(c).Y);
    
end

close all;

%% cross-validate across conditions to assess generalizability
for c = 1:2
    switch c
        case 1
            train = 1; test = 2;
        case 2
            train = 2; test = 1;
    end
    pred(c).vals = tmp(train).constant + tmp(train).weight_obj.dat'*f(test).dat;
end

%% combine across conditions and plot
out.yfit = horzcat(pred.vals)';
out.Y = horzcat(f.Y)';

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
    
%% generate file name 
method = 'lassopcr';
fname = sprintf('predict/idiographic_bycondition/%s/%s/%d_%s_%s', method, physio, id, method, physio);

%% save output
save(sprintf('%s.mat', fname), 'out');
%saveas(gcf, sprintf('%s.png', fname));


