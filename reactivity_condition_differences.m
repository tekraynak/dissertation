function out = reactivity_condition_differences(id, physio)
if nargin < 1, id = 1002; end
if nargin < 2, physio = 'hr'; end

%%
wdir = cdtodrive;
dissdir = [wdir '/Users/Thomas/dissertation'];
cd(dissdir)

%% load reactivity
dat = readtable('datasets/pip_reactivity_2020.csv', 'TreatAsEmpty', {'NA'});
if ~any (id == dat.id); return; end
dat = dat(dat.id == id, contains(dat.Properties.VariableNames, physio));

baseline = dat.([physio '_baseline']);
react = table2array(dat(:, 6:21)) - baseline';

%% calculate ttest between conditions
[~, p, ci, stats] = ttest(react(1:2:end), react(2:2:end));
out.tstat = stats.tstat;
out.p = p;
