wdir = cdtodrive; cd Users/Thomas/dissertation
addpath('masks');

level1dir = [wdir '/PIP/SPM12/First_Level'];

dat = readtable('datasets/analytic_sample_raw.csv');
ids = dat.id;

Tasks = {'Stroop' 'MSIT'};

%% load split-level data for each task and calculate internal consistency
% code adapted from IMT_reliability.m

for t = 1:2 % do the same for each task
    clear sl corrAB 
    for i = 1:length(ids)
        sl{i, 1} = filenames(sprintf('%s/%d/SplitHalf/%s_1stHalf/con_0001.nii', level1dir, ids(i), Tasks{t}), 'char');
        sl{i, 2} = filenames(sprintf('%s/%d/SplitHalf/%s_2ndHalf/con_0001.nii', level1dir, ids(i), Tasks{t}), 'char');
    end
    blockA = fmri_data(sl(:, 1), which('grey.nii'));
    blockB = fmri_data(sl(:, 2), which('grey.nii'));
    
    %% generate structure
   clear dat
   dat(:,:,1) = blockA.dat; dat(:,:,2) = blockB.dat;  
   
     %% convert zeros to NaN so they are ignored in the correlations
    dat(dat == 0 ) = NaN;

    %% mask according to how many people have all data (not NaN) in a voxel
    for i = 1:size(dat, 1)
        %voxelmask(i) = sum(all(squeeze(dat(i, : , :)),2));
        voxelmask(i) = sum(any(~isnan(squeeze(dat(i,:,:))), 2));
    end

   %% number of people who must have a value in a voxel for that voxel to be included
    % go with 50% for now
    thresh = round(size(dat,2)/2); 
    
    %% at each voxel, rescale outliers across people
    dat = filloutliers(dat, 'nearest',  2);

    %% at each voxel, rescale outliers and calculate correlation 
    for i = 1:size(dat, 1)
        corrAB(i) = corr(dat(i, :, 1)', dat(i, :, 2)', 'rows', 'pairwise');
    end

    %% calculate SBreliability
    SBreliability = (2*corrAB)./(1+corrAB);

    %% apply mask to ignore bad SB values
    % thresh based on number/proportion of people who must have data in a voxel
    % go with 50% for now
    thresh = round(size(dat,2)/2); 
    SBreliability(voxelmask< thresh) = 0;
        
     %% generate map
    SB(t) = blockA; SB(t).dat = SBreliability';

    %% view and save
    orthviews(SB(t))
    write(SB(t), 'fname', sprintf('reliability/Spearman_Brown_reliability_map_%s.nii', Tasks{t}));


end