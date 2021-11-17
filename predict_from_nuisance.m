function out2 = predict_from_nuisance(id, physio)
% runs LASSOPCR predicting physiology from nuisance variable
% first output is WM + CSF + motion
% second output is WM + CSF

wdir = cdtodrive;
dissdir = [wdir '/Users/Thomas/dissertation'];
cd(dissdir)

%% predict reactivity from nuisance variables

if nargin < 1, id = 1002; end
if nargin < 2, physio = 'sbp'; end


tasks = {'msit' 'stroop'};
onsets_sec = [10 80 150 220 290 360 430 500]; % onsets in seconds
onsets = onsets_sec/2; % onsets in scans
dur = 30; % duration in scans


%% load nuisance variables and divide nuisance variables into block-related averages
clear dat
for t = 1:length(tasks)
    dat{t} = dlmread(sprintf('betaseries/%d/%s/%d_%s_nuisance.txt', id, tasks{t}, id, tasks{t}));
    for b = 1:8
        dat2{t}(b, :) = mean(dat{t}(onsets(b):onsets(b)+30, :));
    end
end

%% concatenate across 2 tasks
nv = vertcat(dat2{1}, dat2{2}); % 16 x 32 array of nuisance variables
wmcsf = nv(:, [7 8 15 16 23 31]);

%% load reactivity

dat = readtable(sprintf('datasets/analytic_sample_reactivity_%s.csv', physio), 'TreatAsEmpty', {'NA'});
dat = dat(dat{:, 1} == id, :);
react = dat{:, 2:end};

%% LASSO-PCR on all nuiscance variables
out = lassopcr_cv(nv, react', 'nfolds', 16)

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

out2(1) = out;

%% LASSO-PCR on WMCSF
out = lassopcr_cv(wmcsf, react', 'nfolds', 16)

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

out2(2) = out;
close all;

