%%
% CALCULATES THE MEAN T1 IMAGE
%%

function meanT1()

spm('defaults', 'fmri');
spm_jobman('initcfg');

% Load the subject IDs
load('../analysis/subjectList.mat');

expMean = '(';
for iSub = 1:length(subjects)
    nii_files = dir(['../data/scan/' subjects{iSub} '/T1_3D_32ch_7min_0001/wms*.nii']);
    matlabbatch{1}.spm.util.imcalc.input(iSub) = ...
        {['../data/scan/' subjects{iSub} '/T1_3D_32ch_7min_0001/' nii_files(1).name]};
    
    if iSub == length(subjects)
        expMean = strcat(expMean,'i',num2str(iSub),')/',num2str(length(subjects)));
    else
        expMean = strcat(expMean,'i',num2str(iSub),'+');
    end
    
end
matlabbatch{1}.spm.util.imcalc.output = 'meanT1';
matlabbatch{1}.spm.util.imcalc.outdir = {'../data/scan/meanT1/'};
matlabbatch{1}.spm.util.imcalc.expression = expMean;
matlabbatch{1}.spm.util.imcalc.var = struct('name', {}, 'value', {});
matlabbatch{1}.spm.util.imcalc.options.dmtx = 0;
matlabbatch{1}.spm.util.imcalc.options.mask = 0;
matlabbatch{1}.spm.util.imcalc.options.interp = 1;
matlabbatch{1}.spm.util.imcalc.options.dtype = 4;

spm_jobman('run', matlabbatch);

end