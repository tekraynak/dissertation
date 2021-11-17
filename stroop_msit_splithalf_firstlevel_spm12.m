function stroop_msit_splithalf_firstlevel_spm12(id, task, whichhalf)
% Run split-level first-level design matrix and contrast analysis for Stroop or MSIT task
% This divides the run in half and runs design matrices separately
% Does this by loading the original design matrix and editing file locations, onset times accordingly
% NEED TO SPECIFY: subject ID, which task, to run on, and which half to run on
if nargin < 2, task = 'Stroop'; end
if nargin < 3, whichhalf = 'both'; end

% can say whichhalf = 'both' to run both (sets recursive loop)
if strcmp(whichhalf, 'both')
    stroop_msit_splithalf_firstlevel_spm12(id, task, 1)
    stroop_msit_splithalf_firstlevel_spm12(id, task, 2)
    return;
end
    
% % Need robust WLS toolbox before running this script
% % http://www.diedrichsenlab.org/imaging/robustWLS_spm12.html
% % http://www.diedrichsenlab.org/download/rWLS_v4.0_SPM12.zip
% % Since ER is an event-related design, and we applied slice-timing correction
% % in the preprocess step, microtime resolution and microtime onset need to be change.
% % microtime resolution = st_nslices (# slices acquired)
% % microtime onset = slice_Ref_ER (slice_reference)

%% cd to subject directory
mywd = cdtodrive; cd PIP/SPM12

%% Create directory for saving output
switch whichhalf; case 1; halfstring = '1st'; case 2; halfstring = '2nd'; end
mkdir(sprintf('First_Level/%d/SplitHalf/%s_%sHalf', id, task, halfstring));

%% Display stats
fprintf(1, '\n\nRunning Split-Half First-Level Design for %d - %s Half\n\n', id, halfstring);

%% Load original design matrix (not SplitHalf)
load(sprintf('First_Level/%d/%s/%d_%s_level1_batch.mat', id, task, id, task));

%% Load functional data, realignment parameters
% file location differs between studies

funcfiles = filenames(sprintf('Preprocessed/%d/visit_baseline/%s/swr*.nii', id, task), 'absolute');
rpfile = char(filenames(sprintf('Preprocessed/%d/visit_baseline/%s/rp*.txt', id, task), 'absolute'));


%% Change design parameters depending on whether it's the 1st or 2nd half

% Load rp file and cut in half
% Cut functional image list in half

% Design matrix changes
%   Trim design matrix onsets such that only trials within the time frame are allowed
%   'whichtrials' variable says which trials are to be retained for each condition



rp = load(rpfile);
design = matlabbatch{1}.spm.tools.rwls.fmri_rwls_spec.sess.cond;

switch whichhalf
    case 1
        rp = rp(1:140, :);
        funcfiles = funcfiles(1:140);
       % whichtrials = horzcat(design(1:3).onset) <= 172;
    case 2
        rp = rp(141:end, :);
        funcfiles = funcfiles(141:end);
        %whichtrials = horzcat(design(1:3).onset) >= 172;
end

% for Stroop and MSIT, the 1st and 2nd halves are identical to each other,
%   so can simply trim latter half of design and use it for either the 1st
%   half or 2nd half data
design(1).onset = design(1).onset(1:4); % fixation - from 8 blocks to 4
design(2).onset = design(2).onset(1:2); % incongruent - from 4 blocks to 2
design(3).onset = design(3).onset(1:2); % congruent - from 4 blocks to 2

% 
% Trim trials from each condition
% for c = 1:3 % for each of the 3 conditions
%     design(c).onset = design(c).onset(whichtrials(:, c));
%     design(c+3).onset = design(c+3).onset(whichtrials(:, c));
%     design(c+3).duration = design(c+3).duration(whichtrials(:, c));
%     design(c+6).onset = design(c+6).onset(whichtrials(:, c));
%     design(c+9).onset = design(c+9).onset(whichtrials(:, c));            
% end
% 
% %   For 2nd half, subtract 172sec from each retained onset
% if whichhalf == 2
%     for i = 1:length(design)
%         design(i).onset = design(i).onset - 140;
%     end
% end
    
% save new rp.txt file in new First_Level folder
save(sprintf('First_Level/%d/SplitHalf/%s_%sHalf/rp.txt', id, task, halfstring), '-ascii', '-double', '-tabs', 'rp')

%% Make changes to matlabbatch for new SplitHalf design

% link to new rp.txt file
matlabbatch{1}.spm.tools.rwls.fmri_rwls_spec.sess.multi_reg = cellstr(sprintf('%s/First_Level/%d/SplitHalf/%s_%sHalf/rp.txt', pwd, id, task, halfstring));

% new directory
matlabbatch{1}.spm.tools.rwls.fmri_rwls_spec.dir = cellstr(sprintf('%s/First_Level/%d/SplitHalf/%s_%sHalf', pwd, id, task, halfstring));

% updated list of functional images
matlabbatch{1}.spm.tools.rwls.fmri_rwls_spec.sess.scans = funcfiles;

% new design matrix
matlabbatch{1}.spm.tools.rwls.fmri_rwls_spec.sess.cond = design;

% update mask location
matlabbatch{1}.spm.tools.rwls.fmri_rwls_spec.mask = cellstr(which('brainMask_MNI152+csf.nii'));

%% save new batch file
save(sprintf('First_Level/%d/SplitHalf/%s_%sHalf/batch.mat', id, task, halfstring), 'matlabbatch');

%% run batch
spm_jobman('run', matlabbatch);
