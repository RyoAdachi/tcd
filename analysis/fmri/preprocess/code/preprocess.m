%%
% PREPROCESSING
% Input is subject ID (e.g. 002)
%%

function preprocess(id)

if nargin == 0
    error('Use argument of subject ID (e.g. 002)');
end
    
spm('defaults', 'fmri');
spm_jobman('initcfg');

%rp_spm = '/Applications/MATLAB_R2013a.app/toolbox/spm12';
rp_spm = '/home/radachi/Documents/MATLAB/toolbox/spm12';

nSlice = 56; % Number of axial slices
TR = 1; % Repetition time
TA = 0; % For multiband, TA is set to 0 here since it is not used
refTime = 0; % In ms
nRun = 4; % Four sessions
fVoxelSize = [2.5 2.5 2.5]; % Voxel size for the functional images (mm)
sVoxelSize = [1 1 1]; % Voxel size for the functional images (mm)
fwhmSize = [8 8 8]; % Smoothing size

% Load the image acquisition time information
load('slice_timing_multiband.mat');

%% SLICE TIMING CORRECTION

for iRun = 1:nRun    
    nii_files = dir(['../data/scan/' id '/rsBOLD_MB_1_000' num2str(iRun+5) '/f*.nii']);
    for j = 1 :size(nii_files,1)
        matlabbatch{1}.spm.temporal.st.scans{iRun}(j,1)  = ...
            {['../data/scan/' id '/rsBOLD_MB_1_000' num2str(iRun+5) '/' nii_files(j).name]};
    end    
end
matlabbatch{1}.spm.temporal.st.nslices = nSlice;
matlabbatch{1}.spm.temporal.st.tr = TR;
matlabbatch{1}.spm.temporal.st.ta = TA;
matlabbatch{1}.spm.temporal.st.so = slice_timing_multiband;
matlabbatch{1}.spm.temporal.st.refslice = refTime;
matlabbatch{1}.spm.temporal.st.prefix = 'a';

%% REALIGN AND RESLICE

matlabbatch{2}.spm.spatial.realign.estwrite.data{1}(1) = ...
    cfg_dep('Slice Timing: Slice Timing Corr. Images (Sess 1)', ...
    substruct('.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('()',{1}, '.','files'));
matlabbatch{2}.spm.spatial.realign.estwrite.data{2}(1) = ...
    cfg_dep('Slice Timing: Slice Timing Corr. Images (Sess 2)', ...
    substruct('.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('()',{2}, '.','files'));
matlabbatch{2}.spm.spatial.realign.estwrite.data{3}(1) = ...
    cfg_dep('Slice Timing: Slice Timing Corr. Images (Sess 3)', ...
    substruct('.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('()',{3}, '.','files'));
matlabbatch{2}.spm.spatial.realign.estwrite.data{4}(1) = ...
    cfg_dep('Slice Timing: Slice Timing Corr. Images (Sess 4)', ...
    substruct('.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('()',{4}, '.','files'));
matlabbatch{2}.spm.spatial.realign.estwrite.eoptions.quality = 0.9;
matlabbatch{2}.spm.spatial.realign.estwrite.eoptions.sep = 4;
matlabbatch{2}.spm.spatial.realign.estwrite.eoptions.fwhm = 5;
matlabbatch{2}.spm.spatial.realign.estwrite.eoptions.rtm = 0; % Default realign to the mean but here to the first by convention
matlabbatch{2}.spm.spatial.realign.estwrite.eoptions.interp = 2;
matlabbatch{2}.spm.spatial.realign.estwrite.eoptions.wrap = [0 0 0];
matlabbatch{2}.spm.spatial.realign.estwrite.eoptions.weight = '';
matlabbatch{2}.spm.spatial.realign.estwrite.roptions.which = [2 1];
matlabbatch{2}.spm.spatial.realign.estwrite.roptions.interp = 4;
matlabbatch{2}.spm.spatial.realign.estwrite.roptions.wrap = [0 0 0];
matlabbatch{2}.spm.spatial.realign.estwrite.roptions.mask = 1;
matlabbatch{2}.spm.spatial.realign.estwrite.roptions.prefix = 'r';

%% COREGISTERATION

matlabbatch{3}.spm.spatial.coreg.estimate.ref(1) = ...
    cfg_dep('Realign: Estimate & Reslice: Mean Image', ...
    substruct('.','val', '{}',{2}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','rmean'));
nii_files = dir(['../data/scan/' id '/T1_3D_32ch_7min_0001/s*.nii']);
matlabbatch{3}.spm.spatial.coreg.estimate.source = ...
    {['../data/scan/' id '/T1_3D_32ch_7min_0001/' nii_files(1).name]};
matlabbatch{3}.spm.spatial.coreg.estimate.other = {''};
matlabbatch{3}.spm.spatial.coreg.estimate.eoptions.cost_fun = 'nmi';
matlabbatch{3}.spm.spatial.coreg.estimate.eoptions.sep = [4 2];
matlabbatch{3}.spm.spatial.coreg.estimate.eoptions.tol = [0.02 0.02 0.02 0.001 0.001 0.001 0.01 0.01 0.01 0.001 0.001 0.001];
matlabbatch{3}.spm.spatial.coreg.estimate.eoptions.fwhm = [7 7];

%% SEGMENTATION

matlabbatch{4}.spm.spatial.preproc.channel.vols(1) = ...
    cfg_dep('Coregister: Estimate: Coregistered Images', ...
    substruct('.','val', '{}',{3}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','cfiles'));
matlabbatch{4}.spm.spatial.preproc.channel.biasreg = 0.001;
matlabbatch{4}.spm.spatial.preproc.channel.biasfwhm = 60;
matlabbatch{4}.spm.spatial.preproc.channel.write = [0 1];
matlabbatch{4}.spm.spatial.preproc.tissue(1).tpm = {[rp_spm '/tpm/TPM.nii,1']};
matlabbatch{4}.spm.spatial.preproc.tissue(1).ngaus = 1;
matlabbatch{4}.spm.spatial.preproc.tissue(1).native = [1 0];
matlabbatch{4}.spm.spatial.preproc.tissue(1).warped = [0 0];
matlabbatch{4}.spm.spatial.preproc.tissue(2).tpm = {[rp_spm '/tpm/TPM.nii,2']};
matlabbatch{4}.spm.spatial.preproc.tissue(2).ngaus = 1;
matlabbatch{4}.spm.spatial.preproc.tissue(2).native = [1 0];
matlabbatch{4}.spm.spatial.preproc.tissue(2).warped = [0 0];
matlabbatch{4}.spm.spatial.preproc.tissue(3).tpm = {[rp_spm '/tpm/TPM.nii,3']};
matlabbatch{4}.spm.spatial.preproc.tissue(3).ngaus = 2;
matlabbatch{4}.spm.spatial.preproc.tissue(3).native = [1 0];
matlabbatch{4}.spm.spatial.preproc.tissue(3).warped = [0 0];
matlabbatch{4}.spm.spatial.preproc.tissue(4).tpm = {[rp_spm '/tpm/TPM.nii,4']};
matlabbatch{4}.spm.spatial.preproc.tissue(4).ngaus = 3;
matlabbatch{4}.spm.spatial.preproc.tissue(4).native = [1 0];
matlabbatch{4}.spm.spatial.preproc.tissue(4).warped = [0 0];
matlabbatch{4}.spm.spatial.preproc.tissue(5).tpm = {[rp_spm '/tpm/TPM.nii,5']};
matlabbatch{4}.spm.spatial.preproc.tissue(5).ngaus = 4;
matlabbatch{4}.spm.spatial.preproc.tissue(5).native = [1 0];
matlabbatch{4}.spm.spatial.preproc.tissue(5).warped = [0 0];
matlabbatch{4}.spm.spatial.preproc.tissue(6).tpm = {[rp_spm '/tpm/TPM.nii,6']};
matlabbatch{4}.spm.spatial.preproc.tissue(6).ngaus = 2;
matlabbatch{4}.spm.spatial.preproc.tissue(6).native = [0 0];
matlabbatch{4}.spm.spatial.preproc.tissue(6).warped = [0 0];
matlabbatch{4}.spm.spatial.preproc.warp.mrf = 1;
matlabbatch{4}.spm.spatial.preproc.warp.cleanup = 1;
matlabbatch{4}.spm.spatial.preproc.warp.reg = [0 0.001 0.5 0.05 0.2];
matlabbatch{4}.spm.spatial.preproc.warp.affreg = 'mni';
matlabbatch{4}.spm.spatial.preproc.warp.fwhm = 0;
matlabbatch{4}.spm.spatial.preproc.warp.samp = 3;
matlabbatch{4}.spm.spatial.preproc.warp.write = [0 1];

%% NORMALIZATION (FUNCTIONAL)

matlabbatch{5}.spm.spatial.normalise.write.subj.def(1) = ...
    cfg_dep('Segment: Forward Deformations', ...
    substruct('.','val', '{}',{4}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','fordef', '()',{':'}));
matlabbatch{5}.spm.spatial.normalise.write.subj.resample(1) =...
    cfg_dep('Realign: Estimate & Reslice: Resliced Images (Sess 1)', ...
    substruct('.','val', '{}',{2}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','sess', '()',{1}, '.','rfiles'));
matlabbatch{5}.spm.spatial.normalise.write.subj.resample(2) = ...
    cfg_dep('Realign: Estimate & Reslice: Resliced Images (Sess 2)', ...
    substruct('.','val', '{}',{2}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','sess', '()',{2}, '.','rfiles'));
matlabbatch{5}.spm.spatial.normalise.write.subj.resample(3) = ...
    cfg_dep('Realign: Estimate & Reslice: Resliced Images (Sess 3)', ...
    substruct('.','val', '{}',{2}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','sess', '()',{3}, '.','rfiles'));
matlabbatch{5}.spm.spatial.normalise.write.subj.resample(4) = ...
    cfg_dep('Realign: Estimate & Reslice: Resliced Images (Sess 4)', ...
    substruct('.','val', '{}',{2}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','sess', '()',{4}, '.','rfiles'));
matlabbatch{5}.spm.spatial.normalise.write.woptions.bb = [-78 -112 -70
    78 76 85];
matlabbatch{5}.spm.spatial.normalise.write.woptions.vox = fVoxelSize;
matlabbatch{5}.spm.spatial.normalise.write.woptions.interp = 4;

%% NORMALIZATION (STRUCTURAL)

matlabbatch{6}.spm.spatial.normalise.write.subj.def(1) = ...
    cfg_dep('Segment: Forward Deformations', ...
    substruct('.','val', '{}',{4}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','fordef', '()',{':'}));
matlabbatch{6}.spm.spatial.normalise.write.subj.resample(1) = ...
    cfg_dep('Segment: Bias Corrected (1)', ...
    substruct('.','val', '{}',{4}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','channel', '()',{1}, '.','biascorr', '()',{':'}));
matlabbatch{6}.spm.spatial.normalise.write.woptions.bb = [-78 -112 -70
    78 76 85];
matlabbatch{6}.spm.spatial.normalise.write.woptions.vox = sVoxelSize;
matlabbatch{6}.spm.spatial.normalise.write.woptions.interp = 4;

%% SMOOTHING

matlabbatch{7}.spm.spatial.smooth.data(1) = ...
    cfg_dep('Normalise: Write: Normalised Images (Subj 1)', ...
    substruct('.','val', '{}',{5}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('()',{1}, '.','files'));
matlabbatch{7}.spm.spatial.smooth.fwhm = fwhmSize;
matlabbatch{7}.spm.spatial.smooth.dtype = 0;
matlabbatch{7}.spm.spatial.smooth.im = 0;
matlabbatch{7}.spm.spatial.smooth.prefix = 's';

spm_jobman('run', matlabbatch);

end
