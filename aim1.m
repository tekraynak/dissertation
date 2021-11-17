
wdir = cdtodrive;
dissdir = [wdir '/Users/Thomas/dissertation'];
cd(dissdir)


%% generate list of IDs in analytic sample
dat = readtable('datasets/analytic_sample_raw.csv');
ids = dat.id;
length(ids)
phys = {'sbp' 'hr'};

%% run idiographic prediction routines on each participant 

%% primary approach - LASSO-PCR

parfor i = 1:length(ids)
    for p = 1:2
        predict_reactivity_idiographic(ids(i), phys{p});
        close all;
    end
end

%% secondary approaches, because the primary approaches were unsatisfactory
methods = {'lassopcr-alpha05' 'predict-elasticnet'  'predict-lassopcr' 'predict-pcr'};
parfor i = 1:length(ids)
     for p = 1:2
         for m = 1:length(methods)
             predict_reactivity_idiographic(ids(i), phys{p}, methods{m});
         end
     end
end

%% collect output 
% matrix is method(4) x physio(2) x subject(247)
clear rho mae bf r2
methods = {'lassopcr' 'predict-elasticnet' 'predict-pcr'};
for m = 1:length(methods)
     for p = 1:2
         rho = nan(size(ids)); rho_pval = rho; mae = rho; bf10 = rho; r2 = rho;
         for i = 1:length(ids)
            load(sprintf('predict/idiographic/%s/%s/%d_%s_%s.mat', methods{m}, phys{p}, ids(i), methods{m}, phys{p}));
            [rho(i), rho_pval(i)] = corr(out.yfit, out.Y, 'type', 'Spearman');
            mae(i) = out.mae;
            bf10(i) = out.bf;
            r2(i) = out.r2ss;
            disp(i)
         end
        bf10 = filloutliers(bf10, 'nearest');% tends to get outliers
        bf01 = 1./bf10;
        t = table(ids, rho, rho_pval, r2, mae, bf10, bf01);
        barplot_columns(t(:, 2:end));
        writetable(t, sprintf('predict/idiographic/%s_%s_stats.csv', methods{m}, phys{p}));
        
    end
    disp(i)
end
   



%% plots of overall performance
for m = 1:length(methods)
    barplot_columns(squeeze(rho(m, :, :))')
    g = gca; 
    set(g, 'XTickLabels', {'SBP' 'HR'}, 'FontSize', 18)
    set(g.YLabel, 'String', 'Spearman rho')
    saveas(gcf, sprintf('predict/idiographic/%s_rho_distribution.png', methods{m}))
    disp('MAE')
    barplot
end

%% correspondence between reactivity
% how much SBP and HR covary across time
% right panel of Eisenbarth Fig 3
cp{1} = readtable('datasets/analytic_sample_reactivity_hr.csv');
cp{2} = readtable('datasets/analytic_sample_reactivity_sbp.csv');
for i = 1:height(cp{1})
    %physio_r(i) = corr(cp{1}{i, 2:end}', cp{2}{i, 2:end}');
    physio_rho(i) = corr(cp{1}{i, 2:end}', cp{2}{i, 2:end}', 'type', 'Spearman');
end
%barplot_columns(physio_r')
barplot_columns(physio_rho')

%% correspondence of performance between physio measures
corr(squeeze(rho(1, 1, :)), squeeze(rho(1, 2, :)), "rows", "pairwise")
figure; scatter(squeeze(rho(1, 1, :)), squeeze(rho(1, 2, :)), [], 'k.'); lsline;

%% plot showing how PCR does as well as LASSOPCR for participants >0 
% and does better for participants <0
figure; subplot(1, 2, 1); title('HR')
scatter(squeeze(rho(1, 1, :)), squeeze(rho(5, 1, :)), [], 'k.');
hold on; plot([-1:.1:1],[-1:.1:1]);
xlabel('Predicted-Observed rho from LASSO-PCR'); 
ylabel('Predicted-Observed rho from PCR');

subplot(1, 2, 2); title('HR')
scatter(squeeze(rho(1, 2, :)), squeeze(rho(5, 2, :)), [], 'k.');
hold on; plot([-1:.1:1],[-1:.1:1], 'k');
xlabel('Predicted-Observed rho from LASSO-PCR'); 
ylabel('Predicted-Observed rho from PCR');


%% what factors relate to whether the model predicts?
% baseline CV, mean reactivity, var reactivity, motion, behavior (ACC RT)

% load performance
clear rho numcomponents
for i = 1:length(ids)
    for p = 1:2
            clear out; 
            load(sprintf('predict/idiographic/lassopcr/%s/%d_lassopcr_%s.mat', phys{p}, ids(i), phys{p}));
            rho(i, p) = out.rho;
            numcomponents(i, p) = out.num_components_from_best_lambda;
    end
end
sbp_components = numcomponents(:, 1);
hr_components = numcomponents(:, 2);

desc = readtable('datasets/analytic_sample_descriptives.csv', 'TreatAsEmpty', 'NA');

load('predict/idiographic/all_output.mat');
sbp_rho = squeeze(rho(1, 1, :));
hr_rho = squeeze(rho(1, 2, :));
id = ids;
pred_t = table(id, sbp_rho, hr_rho);

pred_t = innerjoin(pred_t, desc);
qc_stats = readtable('qc/qc_stats.csv')
pred_t = innerjoin(pred_t, qc_stats);
pred_t.fd = mean([pred_t.stroop_fdmean pred_t.msit_fdmean], 2);

pred_t = innerjoin(pred_t, table(id, sbp_components, hr_components))

% compare people with 0 components retained to thos with >0 components retained
vars = {'sbp_baseline' 'sbp_all_chg' 'sbp_all_sd' 'ACC' 'RT' 'fd'}
for v = 1:length(vars)
    disp(vars{v})
    tmp = pred_t.(vars{v});
    [~,P,~,STATS] = ttest2(tmp(pred_t.sbp_components==0), tmp(pred_t.sbp_components>0))
end

vars = {'hr_baseline' 'hr_all_chg' 'hr_all_sd' 'ACC' 'RT' 'fd'}
for v = 1:length(vars)
    disp(vars{v})
    tmp = pred_t.(vars{v});
    [~,P,~,STATS] = ttest2(tmp(pred_t.hr_components==0), tmp(pred_t.hr_components>0))
end

    
% linear regression of individual differences
fitlm(pred_t, 'sbp_rho ~ sbp_baseline + sbp_all_chg + sbp_all_sd + ACC + RT + fd ')
fitlm(pred_t, 'hr_rho ~ hr_baseline + hr_all_chg + hr_all_sd + ACC + RT + fd')


%% load maps, ttest, threshold, write nii
clear t
for m = 1:length(methods)   
    for p = 1:2
        load('predict/idiographic/predict-elasticnet/hr/1002_predict-elasticnet_hr.mat'); % load to set up fmri_data object
        w = out.weight_obj;
        for i = 1:length(ids)
            load(sprintf('predict/idiographic/%s/%s/%d_%s_%s.mat', methods{m}, phys{p}, ids(i), methods{m}, phys{p}));
            w.dat(:, i) = out.weight_obj.dat;
        end
        t(m, p) = ttest(w); orthviews(t(m, p))
        write(t(m, p), 'fname', sprintf('predict/idiographic/%s_%s_ttest.nii', methods{m}, phys{p}), 'overwrite');
        tfdr = threshold(t(m, p), .05, 'fdr', 'k', 50); orthviews(tfdr)
        write(tfdr, 'thresh', 'fname', sprintf('predict/idiographic/%s_%s_ttest_fdr05_k50.nii', methods{m}, phys{p}), 'overwrite'); 
    end
end

%% compare maps


% load pairs of maps for everyone
clear ww
for i = 1:length(ids)
    for p = 1:2
        load(sprintf('predict/idiographic/lassopcr/%s/%d_lassopcr_%s.mat', phys{p}, ids(i),phys{p}));
        ww(i, p, :) = out.weight_obj.dat;
    end
    disp(i)
end

% ttest on maps
for i = 1:size(ww, 3) % loop through voxels
    [~, map_p(i), ~, stats]  = ttest(squeeze(ww(:, 1, i)), squeeze(ww(:, 2, i)));
    map_t(i) = stats.tstat;
end

st = statistic_image('type', 'T', 'dfe', size(ww, 1) - 1, 'volInfo', fdat.volInfo, 'dat', map_t', 'p', map_p' )

% do searchlight correlation for each participant
%[sc_r(i), sc_dat(i)] = searchlight_correlation(t(1, 1), t(1, 2), 5) % 5 voxel radius

% load pairs of maps for everyone
clear ww
for i = 1:length(ids)
    for p = 1:2
        load(sprintf('predict/idiographic/lassopcr/%s/%d_lassopcr_%s.mat', phys{p}, ids(i),phys{p}));
        ww(i, p) = out.weight_obj;
    end
    disp(i)
end
% run searchlight using parallel (time-consuming)
clear sc_r sc_dat
load('masks/volInfo.mat')
for i = 1:size(ww, 1)
   tmp1 = ww(i, 1); tmp2 = ww(i, 2);
   if ~(all(tmp1.dat == 0) | all(tmp2.dat == 0))
        tmp1.volInfo = volInfo; tmp2.volInfo = volInfo;
        t = tic;
        [sc_r{i}, sc_dat{i}] = searchlight_correlation(tmp1, tmp2, 5) % 5 voxel radius
        toc(t)
        write(sc_dat{i}, 'overwrite', 'fname', sprintf('predict/idiographic/lassopcr/sbp_hr_searchlight/%d.nii', ids(i)));
    end
end

s