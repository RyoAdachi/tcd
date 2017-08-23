%%
% glm54 TESTS THE REGIONS CORRELATING WITH
% (1) IMAGE PRESENTATION
% (2) BUTTON PRESS
% (3) every 250ms
% (4) BOCPD CPP (demean)
% (5) DBM, diffMean (abs)
% (6) Heu DV (|4-5|)
% (7) Decision period (control)
%%

function glm54(subjectName)

clearvars -except subjectName;

spm('defaults', 'fmri');

% You might want to change the mount point
rp = '/home/radachi/research/CD/';

% Load the subject list
load([rp 'analysis/subjectList.mat']);
iSub = find(ismember(subjects,subjectName));

% Load the regressor values
load([rp 'analysis/beh/fmri_parameter_create/result/parameter_concat.mat']);

%% 1st level

spm_jobman('initcfg');

% Create the directories for SPM.mat to be saved
if ~exist([rp 'analysis/fmri/glm/glm54/result'], 'dir')
    mkdir([rp 'analysis/fmri/glm/glm54/'], 'result');
end

if ~exist([rp 'analysis/fmri/glm/glm54/result/1stlevel'], 'dir')
    mkdir([rp 'analysis/fmri/glm/glm54/result/'], '1stlevel');
end

if ~exist([rp 'analysis/fmri/glm/glm54/result/1stlevel/' subjectName], 'dir')
    mkdir([rp 'analysis/fmri/glm/glm54/result/1stlevel/'],subjectName);
else
    rmdir([rp 'analysis/fmri/glm/glm54/result/1stlevel/' subjectName],'s');
    mkdir([rp 'analysis/fmri/glm/glm54/result/1stlevel/'],subjectName);
end

if strcmp(subjectName,'012')
    vecSess = [1,3,4];
elseif strcmp(subjectName,'023')
    vecSess = [2,3,4];
else
    vecSess = [1,2,3,4];
end

iJob = 1;

matlabbatch{iJob}.spm.stats.fmri_spec.dir = ...
    {[rp 'analysis/fmri/glm/glm54/result/1stlevel/' subjectName '/']};
matlabbatch{iJob}.spm.stats.fmri_spec.timing.units = 'secs';
matlabbatch{iJob}.spm.stats.fmri_spec.timing.RT = 1;
matlabbatch{iJob}.spm.stats.fmri_spec.timing.fmri_t = 14;
matlabbatch{iJob}.spm.stats.fmri_spec.timing.fmri_t0 = 1;
% Put scanning data
countImage = 0;
for iRun = 1:length(vecSess)
    nii_files = dir([rp 'data/scan/' subjectName ...
        '/rsBOLD_MB_1_000' num2str(vecSess(iRun)+5) '/', 'swra*.nii']);
    for j = 1:size(nii_files,1)
        countImage = countImage + 1;
        matlabbatch{iJob}.spm.stats.fmri_spec.sess(1).scans{countImage} = ...
            [rp 'data/scan/' subjectName ...
            '/rsBOLD_MB_1_000' num2str(vecSess(iRun)+5) '/' nii_files(j).name];
    end
end

%% Image presentation (main/control pooled)

iCond = 1;

matlabbatch{iJob}.spm.stats.fmri_spec.sess(1).cond(iCond).name = ...
    parameter.subjects(iSub).run(1).onset(1).name;
matlabbatch{iJob}.spm.stats.fmri_spec.sess(1).cond(iCond).onset = ...
    parameter.subjects(iSub).run(1).onset(1).onset;
matlabbatch{iJob}.spm.stats.fmri_spec.sess(1).cond(iCond).duration = ...
    parameter.subjects(iSub).run(1).onset(1).duration;
matlabbatch{iJob}.spm.stats.fmri_spec.sess(1).cond(iCond).tmod = ...
    parameter.subjects(iSub).run(1).onset(1).tmod;

matlabbatch{iJob}.spm.stats.fmri_spec.sess(1).cond(iCond).pmod = ...
    struct('name', {}, 'param', {}, 'poly', {});
matlabbatch{iJob}.spm.stats.fmri_spec.sess(1).cond(iCond).orth = 0;

%% Button press (main/control pooled)

iCond = iCond + 1;

matlabbatch{iJob}.spm.stats.fmri_spec.sess(1).cond(iCond).name = ...
    parameter.subjects(iSub).run(1).onset(2).name;
matlabbatch{iJob}.spm.stats.fmri_spec.sess(1).cond(iCond).onset = ...
    parameter.subjects(iSub).run(1).onset(2).onset;
matlabbatch{iJob}.spm.stats.fmri_spec.sess(1).cond(iCond).duration = ...
    parameter.subjects(iSub).run(1).onset(2).duration;
matlabbatch{iJob}.spm.stats.fmri_spec.sess(1).cond(iCond).tmod = ...
    parameter.subjects(iSub).run(1).onset(2).tmod;

matlabbatch{iJob}.spm.stats.fmri_spec.sess(1).cond(iCond).pmod = ...
    struct('name', {}, 'param', {}, 'poly', {});
matlabbatch{iJob}.spm.stats.fmri_spec.sess(1).cond(iCond).orth = 0;

%% DV

iCond = iCond + 1; iPmod = 1;

matlabbatch{iJob}.spm.stats.fmri_spec.sess(1).cond(iCond).name = ...
    parameter.subjects(iSub).run(1).onset(7).name;
matlabbatch{iJob}.spm.stats.fmri_spec.sess(1).cond(iCond).onset = ...
    parameter.subjects(iSub).run(1).onset(7).onset;
matlabbatch{iJob}.spm.stats.fmri_spec.sess(1).cond(iCond).duration = ...
    parameter.subjects(iSub).run(1).onset(7).duration;
matlabbatch{iJob}.spm.stats.fmri_spec.sess(1).cond(iCond).tmod = ...
    parameter.subjects(iSub).run(1).onset(7).tmod;

matlabbatch{iJob}.spm.stats.fmri_spec.sess(1).cond(iCond).pmod(iPmod) = ...
    parameter.subjects(iSub).run(1).pmod(9);

iPmod = iPmod + 1;

matlabbatch{iJob}.spm.stats.fmri_spec.sess(1).cond(iCond).pmod(iPmod) = ...
    parameter.subjects(iSub).run(1).pmod(22);

iPmod = iPmod + 1;

matlabbatch{iJob}.spm.stats.fmri_spec.sess(1).cond(iCond).pmod(iPmod) = ...
    parameter.subjects(iSub).run(1).pmod(43);

matlabbatch{iJob}.spm.stats.fmri_spec.sess(1).cond(iCond).orth = 0;

%% Decision period (control)

iCond = iCond + 1;

matlabbatch{iJob}.spm.stats.fmri_spec.sess(1).cond(iCond).name = ...
    parameter.subjects(iSub).run(1).onset(35).name;
matlabbatch{iJob}.spm.stats.fmri_spec.sess(1).cond(iCond).onset = ...
    parameter.subjects(iSub).run(1).onset(35).onset;
matlabbatch{iJob}.spm.stats.fmri_spec.sess(1).cond(iCond).duration = ...
    parameter.subjects(iSub).run(1).onset(35).duration;
matlabbatch{iJob}.spm.stats.fmri_spec.sess(1).cond(iCond).tmod = ...
    parameter.subjects(iSub).run(1).onset(35).tmod;

matlabbatch{iJob}.spm.stats.fmri_spec.sess(1).cond(iCond).pmod = ...
    struct('name', {}, 'param', {}, 'poly', {});
matlabbatch{iJob}.spm.stats.fmri_spec.sess(1).cond(iCond).orth = 0;

%% Movement regressor

matlabbatch{iJob}.spm.stats.fmri_spec.sess(1).multi = {''};
matlabbatch{iJob}.spm.stats.fmri_spec.sess(1).regress = ...
    struct('name', {}, 'val', {});

matlabbatch{iJob}.spm.stats.fmri_spec.sess(1).multi_reg{1} = ...
    [rp 'analysis/fmri/glm/nuisance/result/' subjectName '/nuisance.mat'];

% HPF
matlabbatch{iJob}.spm.stats.fmri_spec.sess(1).hpf = 128;

%% Derivatives, autoregressive model

matlabbatch{iJob}.spm.stats.fmri_spec.fact = struct('name', {}, 'levels', {});
matlabbatch{iJob}.spm.stats.fmri_spec.bases.hrf.derivs = [0 0];
matlabbatch{iJob}.spm.stats.fmri_spec.volt = 1;
matlabbatch{iJob}.spm.stats.fmri_spec.global = 'None';
matlabbatch{iJob}.spm.stats.fmri_spec.mask = ...
    {[rp 'analysis/fmri/mask/image/wholebrain_binary.nii']};
matlabbatch{iJob}.spm.stats.fmri_spec.cvi = 'AR(1)';

%% Estimation

iJob = iJob + 1;

matlabbatch{iJob}.spm.stats.fmri_est.spmmat = ...
    {[rp 'analysis/fmri/glm/glm54/result/1stlevel/' subjectName '/SPM.mat']};
matlabbatch{iJob}.spm.stats.fmri_est.write_residuals = 0;
matlabbatch{iJob}.spm.stats.fmri_est.method.Classical = 1;

%% Run batch

spm_jobman('run', matlabbatch);

clear matlabbatch;

end