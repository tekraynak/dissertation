function idiographic_maps_searchlight(id, varargin)

wdir = cdtodrive;
dissdir = [wdir '/Users/Thomas/dissertation'];
cd(dissdir)
load('masks/volInfo.mat')
phys = {'sbp' 'hr'};

%% load data
for i = 1:length(varargin)
    switch(varargin{i})
        case 'sbp'
            load(sprintf('predict/idiographic/lassopcr/sbp/%d_lassopcr_sbp.mat', id));
            tmp{i} = out.weight_obj;
            if all(tmp{i}.dat == 0)
                disp('map contains all zeros'); return; 
            end
        case 'hr'
            load(sprintf('predict/idiographic/lassopcr/hr/%d_lassopcr_hr.mat', id));
            tmp{i} = out.weight_obj;
            if all(tmp{i}.dat == 0)
                disp('map contains all zeros'); return; 
            end
        case 'eisen'
            tmp{i} = fmri_data('signatures/ANS_Eisenbarth_JN_2016_HR_pattern.img', 'masks/grey.nii');
        case 'jaha'
            tmp{i} = fmri_data('signatures/avgChgSBP_onAvgTaskHE_weight.nii', 'masks/grey.nii');
    end
end
        

tmp{1}.volInfo = volInfo;
tmp{2}.volInfo = volInfo;

%% run searchlight
t = tic;
[searchlight_r, searchlight_dat] = searchlight_correlation(tmp{1}, tmp{2}, 5); % 5 voxel radius
toc(t)

%% write NII image
write(searchlight_dat, 'overwrite', 'fname', ...
    sprintf('predict/idiographic/lassopcr/sbp_hr_searchlight/%s_%s_%d.nii', varargin{1}, varargin{2}, id));
