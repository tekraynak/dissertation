%% check signal coverage of Stroop and MSIT functional data
% do this in 2 ways:
% 1 (as proposed) calculate temporal variance map of smoothed images 
%       and calculate proportion of voxels with nonzero variation
% 2 (additional) calculate number of voxels in mask image 
%       produced by SPM first-level analysis

%% change directory
drive = cdtodrive; 
preprocdir = [drive '/PIP/SPM12/Preprocessed'];
cd Users/Thomas/dissertation; dissdir = pwd;

addpath('masks')

cd betaseries

%% get ids
ids = filenames('*');
tasks = {'Stroop' 'MSIT'};

%% construct mask of temporal variance, compute ratio of voxels with nonzero variance

for i = 1:length(ids)
    for t = 1:2
        func = fmri_data(filenames(sprintf('%s/%d/visit_baseline/%s/swr*.nii', preprocdir, str2num(ids{i}), tasks{t})), which('grey.nii'));
        func_var_perc(i, t) = mean(var(func.dat') ~= 0);
    end
end

%% save this file
save('func_var_perc.mat', 'ids', 'func_var_perc')

%%  get number of voxels in SPM mask

for i = 1:length(ids)
    for t = 1:2
        r = region(fmri_data(sprintf('%d/%s/mask.nii', str2num(ids{i}), tasks{t}), which('grey.nii')));
        spmmask_numvoxels(i, t )= r.numVox;
    end
end