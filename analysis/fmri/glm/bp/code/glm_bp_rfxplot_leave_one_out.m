%% Non intractive mode of rfxplot

%%

function glm_bp_rfxplot_leave_one_out()

rp = '/home/radachi/research/CD/analysis/fmri/glm';

load(['/home/radachi/research/CD/analysis/subjectList.mat']);

% Define globalmax for each ROI (from 2ndlevel using all subjects)
ROINAME = {'dlPFCR','IPSR'};
MAXCOORDINATE = {[43,40,22],[33,-40,32]}; % Peak using all particiapnts within independent ROI

beta = zeros(21,2);

for iRoiName = 1:length(ROINAME)
    
    for iSub = 1:numel(subjects)
        
        %%
        % Get coordinates used for each subject
        % (nearest max to global max excluding subject)
        
        % Load SPM.mat
        load([rp '/glm31/result/2ndlevel_leave_one_out/' ...
            subjects{iSub} '/BP(main-control)/SPM.mat']);
        % Load template xSPM for 2nd level contrast
        load([rp '/glm_rfxplot/2ndlevel_xSPM.mat']);
        xSPM.swd = [rp '/glm31/result/2ndlevel_leave_one_out/' ...
            subjects{iSub} '/BP(main-control)'];
        xSPM.title = 'BP(main-control)';
        xSPM.thresDesc = 'none';
        
        [hreg, xSPM, SPM] = spm_results_ui('Setup', xSPM);
        
        spm_mip_ui('SetCoords', MAXCOORDINATE{iRoiName});
        
        nearestMax(iSub,:) = spm_mip_ui('Jump', spm_mip_ui('FindMIPax'), 'nrmax');
        
        clear xSPM SPM
        
        %% rfxplot
        
        for iCol = 1:2
            
            % Load RFX SPM.mat
            rfxdir = [rp '/glm31/result/2ndlevel/BP(main-control)'];
            load(fullfile(rfxdir,'SPM.mat'));
            
            % Create xSPM struct
            load([rp '/glm_rfxplot/rfx_bp_xSPM.mat']);
            xSPM.swd = rfxdir;
            xSPM.title = 'BP(main-control)';
            
            % define RFX coordinates
            rfxxyz = nearestMax(iSub,:);
            
            % Load a previously saved opts struct
            if iCol == 1
                load([rp '/glm_rfxplot/rfx_opts_bp_main_leave_one_out.mat']);
            else
                load([rp '/glm_rfxplot/rfx_opts_bp_control_leave_one_out.mat']);
            end
            opts.group = {iSub};
            
            % Call rfxplot non-interactively
            data = rfxplot(rfxxyz,opts,SPM,xSPM);
            
            % Save the result
            beta(iSub,iCol) = data(iSub).effect{:};
            
            clear xSPM SPM
            
        end
        
    end
    
    save([rp '/glm_rfxplot/result/glm31_BP(main-control)_' ...
        ROINAME{iRoiName} '_leave_one_out.mat'], 'beta', 'nearestMax');
    
end

end