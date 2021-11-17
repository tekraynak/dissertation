%% make dataset for dissertation study
% data quality check and construct dataset

%% participants are to be included if they:
%  1. have at least 3 presetressor physio (SBP & HR) readings and all task readings
%  2. have complete and useable fMRI data for both Stroop and MSIT
%       (check signal coverage and behavioral accuracy)

drive = cdtodrive; cd Users/Thomas/dissertation

%% first, see who is missing reactivity data
%% load reactivity spreadsheet
dat = readtable('datasets/pip_bp_raw.csv', 'TreatAsEmpty', 'NA');
dat = sortrows(dat,'id','ascend');
ids331 = dat.id;

%% remove DBP and MAP variables
dat(:, contains(dat.Properties.VariableNames, {'dbp' 'map'})) = [];


%% remove 2098 -  later reported taking antihypertensives
dat(dat.id == 2098, :) = [];

%% now look to see who has useable fmri data
missing_funcs = zeros(size(dat.id))
for i = 1:length(dat.id)
    if ~exist(sprintf('betaseries/%d/msit/beta_0001.nii', dat.id(i))) & ...
            ~exist(sprintf('betaseries_incomplete/%d/msit/beta_0001.nii', dat.id(i))) 
        fprintf(1, '\n%d missing msit data', dat.id(i));
        missing_funcs(i) = 1;
    end
    if ~exist(sprintf('betaseries/%d/stroop/beta_0001.nii', dat.id(i))) & ...
            ~exist(sprintf('betaseries_incomplete/%d/stroop/beta_0001.nii', dat.id(i)))
        fprintf(1, '\n%d missing stroop data', dat.id(i));
        missing_funcs(i) = 1;
    end
end
dat.id(find(missing_funcs)); sum(missing_funcs)
        % 1005 missing stroop data - incomplete volumes 
        
        % 2073 missing stroop data - incomplete volumes, cannot reconstruct
        
        % 2101 missing msit data - FIXED
        
        % 2107 missing msit data - incomplete volumes, cannot reconstruct
        
        % 2108 missing msit data - incomplete volumes, cannot reconstruct
        % 2108 missing stroop data - incomplete volumes, cannot reconstruct
        
        % 2121 missing stroop data - FIXED
        
        % 3150 missing msit data - incomplete volumes, cannot reconstruct
        
        % 4293 missing msit data - FIXED
        % 4293 missing stroop data - FIXED
        
        % 6547 missing msit data - FIXED
        % 6547 missing stroop data - FIXED

        % 7581 missing msit data - missing behavioral data, possibly did
        %       not complete both tasks - exclude
        
        % 7614 missing msit data 
        % 7614 missing stroop data - missing behavioral data, possibly did
        %       not complete both tasks - exclude

%% remove the above participants with missing fmri data
dat = dat(~missing_funcs, :);
size(dat)

%% review signal coverage, motion, and behavioral performance
%   see qc.m and related output in qc folder

% tasks = {'Stroop' 'MSIT'};
% for i = 1:height(dat)
%     for t = 1:2
%         fd = framewisedisplacement(textread(char(filenames(sprintf('%s/PIP/SPM12/Preprocessed/%d/visit_baseline/%s/rp*.txt', drive, dat.id(i), tasks{t})))));
%         fdmean(i, t) = fd.fdmean;
%     end
% end

%% remove 6425 and 6429 for poor behavioral performance (less than chance on easy trials)
dat(dat.id == 6425 | dat.id == 6429, :) = []; size(dat)

%% remove for poor signal coverage
dat(dat.id == 5384 | dat.id == 6548 | dat.id == 3166, :) = []; size(dat)

%% identify participants who do not have at least 3 prestressor physio readings or do not have all task readings
% do for SBP and HR separately
sbp = dat(:, contains(dat.Properties.VariableNames, 'sbp'));
sbp.numbaseline = sum(~isnan(sbp{:, 1:5}), 2);
sbp.numtask = sum(~isnan(sbp{:, 6:21}), 2);

hr = dat(:, contains(dat.Properties.VariableNames, 'hr'));
hr.numbaseline = sum(~isnan(hr{:, 1:5}), 2);
hr.numtask = sum(~isnan(hr{:, 6:21}), 2);

missing_cv = sbp.numbaseline < 3 | sbp.numtask ~= 16 | hr.numbaseline < 3 | hr.numtask ~= 16;
% turns out everyone has at least 3 baseline readings, but many have missing task readings

%% report how many participants were removed at this stage and 
fprintf(1, '\n%d participants removed for having missing reactivity task values\n', sum(missing_cv)); 
% 75

%%   get counts of participants who have varying levels of missingness
[cnt_unique, unique_a] = hist(sbp.numtask, unique(sbp.numtask));
table(unique_a, cnt_unique', 'VariableNames', {'Num_SBP_readings' 'Num_participants'})

%     Num_SBP_readings    Num_participants
%     ________________    ________________
% 
%             9                   1       
%            11                   2       
%            12                   5       
%            13                   7       
%            14                  24       
%            15                  36       
%            16                 242       

[cnt_unique, unique_a] = hist(hr.numtask, unique(hr.numtask));
table(unique_a, cnt_unique', 'VariableNames', {'Num_HR_readings' 'Num_participants'})

%     Num_HR_readings    Num_participants
%     _______________    ________________
% 
%            9                   1       
%           11                   2       
%           12                   5       
%           13                   7       
%           14                  24       
%           15                  35       
%           16                 243       

%% write list of participants who have missing reactivity values (to check paper files)
writetable(dat(missing_cv, :), 'datasets/participants_missing_reactivity.csv')

%% reduce to participants with complete reactivity data
dat = dat(~missing_cv, :);



%% write temporary dataset
size(dat)
writetable(dat, 'datasets/analytic_sample_raw.csv')


