%% Contrast creator for glm2
% glm2 TESTS THE REGIONS CORRELATING WITH
% (1) IMAGE PRESENTATION
% (2) BUTTON PRESS
% (3) every 250ms
% (4) BOCPD Var (demean)
% (5) BOCPD CPP (demean)
% (6) Decision period (control)
%%

function glm2_contrast()

clear;

spm('defaults', 'fmri');
spm_jobman('initcfg');

% You might want to change the mount point
rp = '/home/radachi/research/CD/';

% Load the subject list
load([rp 'analysis/subjectList.mat']);

% Contrast names to be saved
contrast = {'BOCPD-CPP(+)','BOCPD-CPP(-)','BOCPD-Var(+)','BOCPD-Var(-)',...
    'every250ms(+)','every250ms(-)'};
    
%% 1st level

for iSub = 1:length(subjects)
    
    %% Contrast
    
    iJob = 1; iCont = 1;
    
    matlabbatch{iJob}.spm.stats.con.spmmat = ...
        {[rp 'analysis/fmri/glm/glm2/result/1stlevel/' subjects{iSub} '/SPM.mat']};
    
    matlabbatch{iJob}.spm.stats.con.consess{iCont}.tcon.name = contrast{1};
    matlabbatch{iJob}.spm.stats.con.consess{iCont}.tcon.weights = [0 0 0 0 1];
    matlabbatch{iJob}.spm.stats.con.consess{iCont}.tcon.sessrep = 'none';
    
    iCont = iCont + 1;
    
    matlabbatch{iJob}.spm.stats.con.consess{iCont}.tcon.name = contrast{2};
    matlabbatch{iJob}.spm.stats.con.consess{iCont}.tcon.weights = [0 0 0 0 -1];
    matlabbatch{iJob}.spm.stats.con.consess{iCont}.tcon.sessrep = 'none';
    
    iCont = iCont + 1;
    
    matlabbatch{iJob}.spm.stats.con.consess{iCont}.tcon.name = contrast{3};
    matlabbatch{iJob}.spm.stats.con.consess{iCont}.tcon.weights = [0 0 0 1];
    matlabbatch{iJob}.spm.stats.con.consess{iCont}.tcon.sessrep = 'none';
    
    iCont = iCont + 1;
    
    matlabbatch{iJob}.spm.stats.con.consess{iCont}.tcon.name = contrast{4};
    matlabbatch{iJob}.spm.stats.con.consess{iCont}.tcon.weights = [0 0 0 -1];
    matlabbatch{iJob}.spm.stats.con.consess{iCont}.tcon.sessrep = 'none';
    
    iCont = iCont + 1;
    
    matlabbatch{iJob}.spm.stats.con.consess{iCont}.tcon.name = contrast{5};
    matlabbatch{iJob}.spm.stats.con.consess{iCont}.tcon.weights = [0 0 1];
    matlabbatch{iJob}.spm.stats.con.consess{iCont}.tcon.sessrep = 'none';
    
    iCont = iCont + 1;
    
    matlabbatch{iJob}.spm.stats.con.consess{iCont}.tcon.name = contrast{6};
    matlabbatch{iJob}.spm.stats.con.consess{iCont}.tcon.weights = [0 0 -1];
    matlabbatch{iJob}.spm.stats.con.consess{iCont}.tcon.sessrep = 'none';
    
    matlabbatch{iJob}.spm.stats.con.delete = 1;
    
    %% Run batch
    
    spm_jobman('run', matlabbatch);
    
    clear matlabbatch;
    
end

%% 2nd level

% Create the directories for SPM.mat to be saved
if ~exist([rp 'analysis/fmri/glm/glm2/result/2ndlevel'], 'dir')
    mkdir([rp 'analysis/fmri/glm/glm2/result/'], '2ndlevel');
else
    rmdir([rp 'analysis/fmri/glm/glm2/result/2ndlevel'],'s');
    mkdir([rp 'analysis/fmri/glm/glm2/result/'], '2ndlevel');
end

for iCon = 1:length(contrast)
    
    if ~exist([rp 'analysis/fmri/glm/glm2/result/2ndlevel/' ...
            contrast{iCon}], 'dir')
        mkdir([rp 'analysis/fmri/glm/glm2/result/2ndlevel/' ...
            contrast{iCon}]);
    end
    
    %% SPECIFICATION
    
    iJob = 1;
    
    matlabbatch{iJob}.spm.stats.factorial_design.dir = ...
        {[rp 'analysis/fmri/glm/glm2/result/2ndlevel/' contrast{iCon}]};
    
    for iSub = 1:numel(subjects)
        if iCon <= 9
            matlabbatch{iJob}.spm.stats.factorial_design.des.t1.scans{iSub} = ...
                [rp 'analysis/fmri/glm/glm2/result/1stlevel/' ...
                subjects{iSub} '/con_000' num2str(iCon) '.nii'];
        else
            matlabbatch{iJob}.spm.stats.factorial_design.des.t1.scans{iSub} = ...
                [rp 'analysis/fmri/glm/glm2/result/1stlevel/' ...
                subjects{iSub} '/con_00' num2str(iCon) '.nii'];
        end
    end
    
    matlabbatch{iJob}.spm.stats.factorial_design.cov = ...
        struct('c', {}, 'cname', {}, 'iCFI', {}, 'iCC', {});
    matlabbatch{iJob}.spm.stats.factorial_design.masking.tm.tm_none = 1;
    matlabbatch{iJob}.spm.stats.factorial_design.masking.im = 1;
    matlabbatch{iJob}.spm.stats.factorial_design.masking.em = {''};
    matlabbatch{iJob}.spm.stats.factorial_design.globalc.g_omit = 1;
    matlabbatch{iJob}.spm.stats.factorial_design.globalm.gmsca.gmsca_no = 1;
    matlabbatch{iJob}.spm.stats.factorial_design.globalm.glonorm = 1;
    
    %% ESTIMATION
    
    iJob = iJob + 1;
    
    matlabbatch{iJob}.spm.stats.fmri_est.spmmat = ...
        {[rp 'analysis/fmri/glm/glm2/result/2ndlevel/' contrast{iCon} '/SPM.mat']};
    matlabbatch{iJob}.spm.stats.fmri_est.write_residuals = 0;  
    matlabbatch{iJob}.spm.stats.fmri_est.method.Classical = 1;
    
    %% CONTRAST
    
    iJob = iJob + 1;
    
    matlabbatch{iJob}.spm.stats.con.spmmat = ...
        {[rp 'analysis/fmri/glm/glm2/result/2ndlevel/' contrast{iCon} '/SPM.mat']};
    matlabbatch{iJob}.spm.stats.con.consess{1}.tcon.name = contrast{iCon};
    matlabbatch{iJob}.spm.stats.con.consess{1}.tcon.convec = 1;
    matlabbatch{iJob}.spm.stats.con.consess{1}.tcon.sessrep = 'none';
    matlabbatch{iJob}.spm.stats.con.delete = 1;
    
    
    %% RUN SPM
    spm_jobman('run', matlabbatch);
    
    clear matlabbatch;
    
end

end