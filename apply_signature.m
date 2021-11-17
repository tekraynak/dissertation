function out = apply_signature(id, method, doplot)
if nargin < 1, id = 1004; end
if nargin < 2, method = 'dot'; end
if nargin < 3, doplot = 0; end

%%
wdir = cdtodrive;
dissdir = [wdir '/Users/Thomas/dissertation'];
cd(dissdir)
out.id = id;

%% %% load reactivity
sbp = readtable('datasets/analytic_sample_reactivity_sbp.csv');
sbp = sbp{sbp{:, 1} == id, 2:end};
hr = readtable('datasets/analytic_sample_reactivity_hr.csv');
hr = hr{hr{:, 1} == id, 2:end};
out.physio = {'hr' 'sbp'};
out.Y = [hr' sbp'];

%% check for images
if ~exist(sprintf('betaseries/%d/msit', id)) || ~exist(sprintf('betaseries/%d/stroop', id))
    out.rho = NaN;
    out.rho_incongruent = NaN;
    out.rho_congruent = NaN;
    out.numblocks = NaN;
    return
end

%% name images
for i = 1:8
    img(i) = filenames(sprintf('betaseries/%d/msit/beta_000%d.nii', id, i), 'absolute');
    img(i+8) = filenames(sprintf('betaseries/%d/stroop/beta_000%d.nii', id, i), 'absolute');
end

%% remove NaN data
%img = img(~isnan(react))';
%react = react(~isnan(react))';

%% load images
f = fmri_data(img, 'masks/grey.nii');

%% load signatures
out.signatures = {'Eisenbarth 2016 (HR)' 'Gianaros 2017 (SBP)'};
sig(1) = fmri_data('signatures/ANS_Eisenbarth_JN_2016_HR_pattern.img', 'masks/grey.nii');
sig(2) = fmri_data('signatures/avgChgSBP_onAvgTaskHE_weight.nii', 'masks/grey.nii');

%% apply signatures (dot-product) to beta series images
for i = 1:2
    switch method
        case 'dot'
            out.yfit(:, i) = canlab_pattern_similarity(f.dat, sig(i).dat);
        case 'cosine'
            out.yfit(:, i) = canlab_pattern_similarity(f.dat, sig(i).dat, 'cosine_similarity');
    end
end
% also can try cosine similarity;

%% calculate metrics for each pair of Y and yfit
for YY = 1:2
    for yy = 1:2
        tmp_Y = out.Y(:, YY);
        tmp_yfit = out.yfit(:, yy);
        
        %% calculate predicted-observed correlations 
        out.rho(YY, yy) =  corr(tmp_Y, tmp_yfit, 'type', 'Spearman');
        % rho is a matrix showing Spearman's row between each predicted and observed pair
        % rows = observed variable (top = HR, bottom = SBP)
        % columsn = signature (left = Eisenbarth, right = Gianaros)

        %% calculate MAE
        out.mae(YY, yy) = mean(abs(tmp_Y - tmp_yfit));

        %% calculate Bayes Factor
        % https://klabhub.github.io/bayesFactor/
        out.bf(YY, yy) = bf.regression(tmp_Y, tmp_yfit);

        %% calculate R2 - sum of squares
        tmp_Y = scale(tmp_Y);
        tmp_yfit = scale(tmp_yfit);
        rss = sum((tmp_Y - tmp_yfit).^2);
        tss = sum((tmp_Y - mean(tmp_Y)).^2);
        out.r2ss(YY, yy) = 1 - (rss/tss);
    end
end


%% scatterplots showing predicted & observed, color according to condition
if doplot
    figure; 
    subplot(1, 2, 1);
    scatter(out(1).predicted)
    hold on; scatter(predicted(2:2:end), out.observed(2:2:end), 'b', '*'); lsline;
    title(sprintf('overall rho = %.2f\n incongruent rho = %.2f\ncongruent rho = %.2f', ...
        out.rho, out.rho_incongruent, out.rho_congruent));
end

%% save
save(sprintf('predict/signatures/%s/%d.mat', method, id), 'out');

