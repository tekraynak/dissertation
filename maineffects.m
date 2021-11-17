


drive = cdtodrive;

cd /home/tek31/ProjectDrive/Users/Thomas/dissertation


%% generate list of IDs in analytic sample
dat = readtable('datasets/analytic_sample_raw.csv');
ids = dat.id;
length(ids)
phys = {'sbp' 'hr'};


%% load data
tasks = {'Stroop' 'MSIT'}

clear tdat
for t = 1:2
    clear fils
    for i = 1:length(ids)
        fils{i} = sprintf('%s/PIP/SPM12/First_Level/%d/%s/con_0001.nii', drive, ids(i), tasks{t});
    end
    tdat = ttest(fmri_data(fils, 'masks/grey.nii'))
    write(tdat, 'fname', sprintf('maps/MainEffects_%s.nii', tasks{t}), 'overwrite')
    thdat = threshold(tdat, .05, 'fdr', 'k', 50);
    write(thdat, 'fname', sprintf('maps/MainEffects_%s_fdr05_k50.nii', tasks{t}), 'thresh', 'overwrite')
end