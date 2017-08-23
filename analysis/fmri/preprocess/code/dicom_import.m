%%
% DICOM IMPORT
% Input is subject ID (e.g. 002)
%%

function dicom_import(id)

spm('defaults', 'fmri');
spm_jobman('initcfg');

% DICOM data
dicom_files = dir(['../data/scan/' id '/dicom/*.dcm']);

for j = 1 :size(dicom_files,1)
    matlabbatch{1}.spm.util.import.dicom.data{j}  = ...
        ['../data/scan/' id '/dicom/' dicom_files(j).name];
end

matlabbatch{1}.spm.util.import.dicom.root = 'series';
matlabbatch{1}.spm.util.import.dicom.outdir = {['../data/scan/' id '/']};
matlabbatch{1}.spm.util.import.dicom.protfilter = '.*';
matlabbatch{1}.spm.util.import.dicom.convopts.format = 'nii';
matlabbatch{1}.spm.util.import.dicom.convopts.icedims = 0;

spm_jobman('run', matlabbatch);

end
