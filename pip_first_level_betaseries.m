function matlabbatch = pip_first_level_betaseries(id)
if nargin < 1, id = 1003; end
% calculates first level GLM for Stroop and MSIT that incorporates block-specific regressors 
% 8 blocks per task
% GLM does not model condition effects, per committee request
% produces a separate GLM for each task, per committee request

%% CD to drive and set paths
wdir = cdtodrive; 
preprocdir = sprintf('%s/PIP/SPM12/Preprocessed/%d/visit_baseline/', wdir, id);

dissdir = '/home/tek31/ProjectDrive/Users/Thomas/dissertation'; % where analyses will be conducted and stored

%% loop through Stroop and MSIT separately

Tasks = {'Stroop' 'MSIT'}; tasks = {'stroop' 'msit'};

for t = 1:2

    %% get block length and start times
    times = load(sprintf('%s/PIP/SummaryFiles/block_length_%s.mat', wdir, tasks{t}));

    if ~any(id == [times.dat.id]); return; end % exit if we cannot find timing for this participant

    blocklength = times.dat([times.dat.id] == id).blocklength;
    
    blockstart(1) = 10;
    for i = 2:8
        blockstart(i) = blockstart(i-1) + 70;
    end

    fixstart(1) = 0; fixlength(1) = 10;
    for i = 2:8
        fixstart(i) = blockstart(i-1) + blocklength(i-1);
        fixlength(i) = 60 - blocklength(i-1) + 10;
    end


    %% get functional and realignment data
    funcs = filenames(sprintf('%s/%s/*swr*', preprocdir, Tasks{t}));

    %% exit if functional data is not correct size
    if size(funcs, 1) ~= 280
            fprintf(1, '\nMISSING OR INCORRECT NUMBER OF SCANS - ABORT\n')
        continue;
    end

    %% make new first level directory
    savedir = sprintf('%s/betaseries/%d/%s', dissdir, id, tasks{t});
    mkdir(savedir)

    %% get realignment data, extract WM+CSF signal, and make new nuisance data file
    rpfile = filenames(sprintf('%s/%s/rp*txt', preprocdir, Tasks{t}));
    rp = load(rpfile{1});
    tmp = extract_from_rois(char(funcs), [dissdir '/masks/espm_white_prob_90.nii']);
    wm = tmp.average_data;
    tmp = extract_from_rois(char(funcs), [dissdir '/masks/espm_csf_prob_75.nii']);
    csf = tmp.average_data;
    nuisance = [rp wm csf];
    nuisance = [nuisance [zeros(1, size(nuisance, 2)); diff(nuisance)]]; % add derivative
    nuisance = [nuisance (nuisance.^2)];
    nuisance_file = sprintf('%s/betaseries/%d/%s/%d_%s_nuisance.txt', dissdir, id, tasks{t}, id, tasks{t});
    dlmwrite(nuisance_file, nuisance);

    %% make SPM batch
    clear matlabbatch
    matlabbatch{1}.spm.stats.fmri_spec.dir = {savedir};
    matlabbatch{1}.spm.stats.fmri_spec.timing.units = 'secs';
    matlabbatch{1}.spm.stats.fmri_spec.timing.RT = 2;
    matlabbatch{1}.spm.stats.fmri_spec.timing.fmri_t = 16;
    matlabbatch{1}.spm.stats.fmri_spec.timing.fmri_t0 = 8;


    matlabbatch{1}.spm.stats.fmri_spec.sess(1).scans = funcs;
    % matlabbatch{1}.spm.stats.fmri_spec.sess(1).cond(1).name = 'Incongruent';
    % matlabbatch{1}.spm.stats.fmri_spec.sess(1).cond(1).onset = blockstart(1:2:8)';
    % matlabbatch{1}.spm.stats.fmri_spec.sess(1).cond(1).duration = stroop_blocklength(1:2:8)';
    % matlabbatch{1}.spm.stats.fmri_spec.sess(1).cond(1).tmod = 0;
    % matlabbatch{1}.spm.stats.fmri_spec.sess(1).cond(1).pmod = struct('name', {}, 'param', {}, 'poly', {});
    % matlabbatch{1}.spm.stats.fmri_spec.sess(1).cond(1).orth = 1;
    % matlabbatch{1}.spm.stats.fmri_spec.sess(1).cond(2).name = 'Congruent';
    % matlabbatch{1}.spm.stats.fmri_spec.sess(1).cond(2).onset = blockstart(2:2:8)';
    % matlabbatch{1}.spm.stats.fmri_spec.sess(1).cond(2).duration = stroop_blocklength(2:2:8)';
    % matlabbatch{1}.spm.stats.fmri_spec.sess(1).cond(2).tmod = 0;
    % matlabbatch{1}.spm.stats.fmri_spec.sess(1).cond(2).pmod = struct('name', {}, 'param', {}, 'poly', {});
    % matlabbatch{1}.spm.stats.fmri_spec.sess(1).cond(2).orth = 1;
    % matlabbatch{1}.spm.stats.fmri_spec.sess(1).cond(3).name = 'Fixation';
    % matlabbatch{1}.spm.stats.fmri_spec.sess(1).cond(3).onset = stroop_fixstart';
    % matlabbatch{1}.spm.stats.fmri_spec.sess(1).cond(3).duration = stroop_fixlength';
    % matlabbatch{1}.spm.stats.fmri_spec.sess(1).cond(3).tmod = 0;
    % matlabbatch{1}.spm.stats.fmri_spec.sess(1).cond(3).pmod = struct('name', {}, 'param', {}, 'poly', {});
    % matlabbatch{1}.spm.stats.fmri_spec.sess(1).cond(3).orth = 1;

    % block-specific regressors
    for i = 1:8
        matlabbatch{1}.spm.stats.fmri_spec.sess(1).cond(i).name = sprintf('Block%d', i);
        matlabbatch{1}.spm.stats.fmri_spec.sess(1).cond(i).onset = blockstart(i);
        matlabbatch{1}.spm.stats.fmri_spec.sess(1).cond(i).duration = blocklength(i);
        matlabbatch{1}.spm.stats.fmri_spec.sess(1).cond(i).tmod = 0;
        matlabbatch{1}.spm.stats.fmri_spec.sess(1).cond(i).pmod = struct('name', {}, 'param', {}, 'poly', {});
        matlabbatch{1}.spm.stats.fmri_spec.sess(1).cond(i).orth = 1;
    end


    matlabbatch{1}.spm.stats.fmri_spec.sess(1).multi = {''};
    matlabbatch{1}.spm.stats.fmri_spec.sess(1).regress = struct('name', {}, 'val', {});
    matlabbatch{1}.spm.stats.fmri_spec.sess(1).multi_reg = {nuisance_file};
    matlabbatch{1}.spm.stats.fmri_spec.sess(1).hpf = 128;

    matlabbatch{2}.spm.stats.fmri_est.spmmat(1) = cfg_dep('fMRI model specification: SPM.mat File', substruct('.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','spmmat'));
    matlabbatch{2}.spm.stats.fmri_est.write_residuals = 0;
    matlabbatch{2}.spm.stats.fmri_est.method.Classical = 1;

    %% save batch
    save(sprintf('%s/batch.mat', savedir), 'matlabbatch');

    %% run batch
    spm_jobman('run', matlabbatch);

end
