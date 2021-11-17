%% routines to check quality in dissertation project
%   these routines are to be run on the n=247 participamts with complete fmri and reactivity data
%   this sample may be reduced further if problematic participants are identified

%% check signal coverage of Stroop and MSIT functional data
% do this in 2 ways:
% 1 (as proposed) calculate temporal variance map of smoothed images 
%       and calculate proportion of voxels with nonzero variation
% 2 (additional) calculate number of voxels in mask image 
%       produced by SPM first-level analysis

%% change directory
[drive, rawdir]  = cdtodrive; 
preprocdir = [drive '/PIP/SPM12/Preprocessed'];
cd Users/Thomas/dissertation; dissdir = pwd;

addpath('masks')

cd betaseries

%% get ids (n=247)
ids = cellfun(@str2num, filenames('*'));
tasks = {'Stroop' 'MSIT'};

%% construct mask of temporal variance, compute ratio of voxels with nonzero variance

for i = 1:length(ids)
    for t = 1:2
        func = fmri_data(filenames(sprintf('%s/%d/visit_baseline/%s/swr*.nii', preprocdir, ids(i), tasks{t})), which('grey.nii'));
        func_var_perc(i, t) = mean(var(func.dat') ~= 0);
    end
end

%% save this file
save('../qc/func_var_perc.mat', 'ids', 'func_var_perc')

%%  get number of voxels in SPM mask

for i = 1:length(ids)
    for t = 1:2
        r = region(fmri_data(sprintf('%d/%s/mask.nii', ids(i), tasks{t}), which('grey.nii')));
        spmmask_numvoxels(i, t )= r.numVox;
    end
end

%% compute framewise displacement for each task
tasks = {'Stroop' 'MSIT'};
for i = 1:length(ids)
    for t = 1:2
        fd = framewisedisplacement(textread(char(filenames(sprintf('%s/PIP/SPM12/Preprocessed/%d/visit_baseline/%s/rp*.txt', ...
            drive, ids(i), tasks{t})))));
        fdmean(i, t) = fd.fdmean;
    end
end

%% combine key variables into 1 table
qc = table;
qc.id = ids;
qc.stroop_variance = func_var_perc(:, 1);
qc.msit_variance = func_var_perc(:, 2);
qc.stroop_mask = spmmask_numvoxels(:, 1);
qc.msit_mask = spmmask_numvoxels(:, 2);
qc.stroop_fdmean = fdmean(:, 1);
qc.msit_fdmean = fdmean(:, 2);
writetable(qc, 'qc_stats.csv')

%% identify outliers
for c = 2:width(qc)
    fprintf(1, '\n\n %s \n', qc.Properties.VariableNames{c});
    tmp = scale(qc{:, c});
    out = find(abs(tmp) > 3);
    for o = 1:length(out)
        fprintf(1, '\n%d z = %.2f', ids(out(o)), tmp(out(o)));
    end
end

%% get behavioral performance
for i = 1:length(ids)
    for t = 1:2
        behav(i, t) = pip_readeprime(char(filenames(sprintf('%s/PIP_*/E-Prime/RawData/Scanner/%d*/%s*.txt', ...
            rawdir, ids(i), tasks{t}))));    
    end
end


%% flag IDs with lower-than-chance accuracy (25% Stroop, 33% MSIT)
for i = 1:size(behav, 1)
    if behav(i, 1).Congruent_ACC < .25 | behav(i, 2).Congruent_ACC < .33
        fprintf(1, '\n%d - %.2f Stroop Cong ACC; %.2f MSIT Con ACC\n', ids(i), behav(i, 1).Congruent_ACC, behav(i, 2).Congruent_ACC);
    end
end

%6425 - 0.92 Stroop ACC; 0.29 MSIT ACC
%6429 - 0.02 Stroop ACC; 0.99 MSIT ACC

