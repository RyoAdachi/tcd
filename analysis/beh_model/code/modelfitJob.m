%%  Fit the behavioral model using fmmincon

function modelfitJob(model_number)

if ~ischar(model_number)
    error('Input must be a string.');
end

% You might want to change the mount point
rp = '/home/radachi/research/CD/';

% Load the subject list
load([rp 'analysis/subjectList.mat']);

matlab_script = 'modelfit';
walltime = '47:59:59';

for iSub = 1:length(subjects)

text=['/home/radachi/research/CD/analysis/beh_model/cd_modelfit_tolmanJob_all.sh ',... 
matlab_script, ' ', walltime, ' ', model_number, ' ', subjects{iSub}];
system(text);
out=1;
    
end
