

if nargin < 1, physio = 'sbp'; end
if nargin < 2, n = 10; end

nfolds = 10;

%%
wdir = cdtodrive;
dissdir = [wdir '/Users/Thomas/dissertation'];
cd(dissdir)

%% load reactivity
dat = readtable(sprintf('datasets/analytic_sample_reactivity_%s.csv', physio), 'TreatAsEmpty', {'NA'});

%% restrict to sample size (total n = 242; testing on smaller samples)
dat = dat(1:n, :);
ids = dat{:,1};
react = dat{:, 2:end};

%% set up image file paths
for s = 1:n
    for i = 1:8
        img(s, i) = filenames(sprintf('betaseries/%d/msit/beta_000%d.nii', ids(s), i), 'absolute');
        img(s, i+8) = filenames(sprintf('betaseries/%d/stroop/beta_000%d.nii', ids(s), i), 'absolute');
    end
end

%% reshape image paths and reactivity variables into arrays
react1 = reshape(react', [16*n, 1]);
img1 = reshape(img', [16*n, 1]);
ids1 = repelem(ids, 16)';
f = fmri_data(img1, which('grey.nii'))
f.Y = react1;

%% mlpcr2
[B, Bb, Bw, pc_b, sc_b, pc_w, sc_w] = mlpcr2(f.dat', f.Y, 'subjIDs', ids1);

%% mlpcr2 with cross-validation using predict
[~, stats] = predict(f, 'algorithm_name', 'cv_mlpcr', 'subjIDs', ids1, 'nfolds', 10)
% i think the between and within maps are stored in stats.other_output - B Bb Bw

[~, stats] = predict(f, 'algorithm_name', 'cv_lassopcr', 'subjIDs', ids1, 'nfolds', 10)

%% try to bootstrap the mlpcr2 method

opt = statset('UseParallel', true);


x = 1;
for i = 1:n
    for b = 1:16
        bootdat(i, b, :) = f.dat(:, x);
        x = x+1;
    end
end
        
bootfun = @(x, y) mlpcr2_bootfun(x, y);
boot =  bootstrp(10, bootfun, bootdat, react);%, 'Options', opt);


function Bout = mlpcr2_bootfun(f)
    Y = f(:, 1); X = f(:, 2:end);
    [B, Bb, Bw, pc_b, sc_b, pc_w, sc_w] = mlpcr2(X,Y,varargin);
    Bout = [Bb' Bw'];
end


%% old mlpcr (doesn't work)
X = f.dat';
Y = react1;
n_subj = n;
subjid = ids1;
 fit_lme_options = {'FitMethod','REML','CovariancePattern','isotropic'}
   [w,I,lme] = mlpcr(X,Y,'topLvl', {ones(length(Y),1),n_subj-1}, ...
                    'withinSubj', {subjid,8}, 'fitlmeoptions', fit_lme_options);
           