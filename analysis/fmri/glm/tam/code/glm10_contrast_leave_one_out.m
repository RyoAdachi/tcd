%% Contrast creator for glm10
% glm10 TESTS THE REGIONS CORRELATING WITH
% (1) IMAGE PRESENTATION
% (2) BUTTON PRESS
% (3) every 250ms
% (4) Heu DV (-X~0) always
% (5) Heu DV (-Y-X~-X) always
% (6) Heu DV (|4-5|)
% (7) Decision period (control)

%%

function glm10_contrast_leave_one_out(subjectName)

clearvars -except subjectName;

spm('defaults', 'fmri');
spm_jobman('initcfg');

% You might want to change the mount point
rp = '/home/radachi/research/CD/';

% Load the subject list
load([rp 'analysis/subjectList.mat']);
iSub = find(ismember(subjects,subjectName));

% Contrast names to be saved
contrast = {'X(+)','X(-)','Y(+)','Y(-)','DV(+)','DV(-)',...
    'every250ms(+)','every250ms(-)'};

% Create the directories for SPM.mat to be saved
if ~exist([rp 'analysis/fmri/glm/glm10/result/2ndlevel_leave_one_out'], 'dir')
    mkdir([rp 'analysis/fmri/glm/glm10/result/'], '2ndlevel_leave_one_out');
end


subject_out = subjects;
subject_out(iSub) = [];

if ~exist([rp 'analysis/fmri/glm/glm10/result/2ndlevel_leave_one_out/' subjectName], 'dir')
    mkdir([rp 'analysis/fmri/glm/glm10/result/2ndlevel_leave_one_out'], subjectName);
else
    rmdir([rp 'analysis/fmri/glm/glm10/result/2ndlevel_leave_one_out/' subjectName],'s');
    mkdir([rp 'analysis/fmri/glm/glm10/result/2ndlevel_leave_one_out'], subjectName);
end

for iCon = 1:length(contrast)
    
    if ~exist([rp 'analysis/fmri/glm/glm10/result/2ndlevel_leave_one_out/' ...
            subjectName '/' contrast{iCon}], 'dir')
        mkdir([rp 'analysis/fmri/glm/glm10/result/2ndlevel_leave_one_out/' ...
            subjectName '/' contrast{iCon}]);
    end
    
    %% SPECIFICATION
    
    iJob = 1;
    
    matlabbatch{iJob}.spm.stats.factorial_design.dir = ...
        {[rp 'analysis/fmri/glm/glm10/result/2ndlevel_leave_one_out/' ...
        subjectName '/' contrast{iCon}]};
    
    for j = 1:numel(subject_out)
        if iCon <= 9
            matlabbatch{iJob}.spm.stats.factorial_design.des.t1.scans{j} = ...
                [rp 'analysis/fmri/glm/glm10/result/1stlevel/' ...
                subject_out{j} '/con_000' num2str(iCon) '.nii'];
        else
            matlabbatch{iJob}.spm.stats.factorial_design.des.t1.scans{j} = ...
                [rp 'analysis/fmri/glm/glm10/result/1stlevel/' ...
                subject_out{j} '/con_00' num2str(iCon) '.nii'];
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
        {[rp 'analysis/fmri/glm/glm10/result/2ndlevel_leave_one_out/' ...
        subjectName '/' contrast{iCon} '/SPM.mat']};
    matlabbatch{iJob}.spm.stats.fmri_est.write_residuals = 0;
    matlabbatch{iJob}.spm.stats.fmri_est.method.Classical = 1;
    
    %% CONTRAST
    
    iJob = iJob + 1;
    
    matlabbatch{iJob}.spm.stats.con.spmmat = ...
        {[rp 'analysis/fmri/glm/glm10/result/2ndlevel_leave_one_out/' ...
        subjectName '/' contrast{iCon} '/SPM.mat']};
    matlabbatch{iJob}.spm.stats.con.consess{1}.tcon.name = contrast{iCon};
    matlabbatch{iJob}.spm.stats.con.consess{1}.tcon.convec = 1;
    matlabbatch{iJob}.spm.stats.con.consess{1}.tcon.sessrep = 'none';
    matlabbatch{iJob}.spm.stats.con.delete = 1;
    
    
    %% RUN SPM
    spm_jobman('run', matlabbatch);
    
    clear matlabbatch;
    
end

end