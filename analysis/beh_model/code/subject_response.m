%%
% This function create the matrix of subject response (1/0) of each bin in 
% a matrix (72*120). 
% Row is ordered using the sequence identifier (1,...,100).

% Output is the following.
% (1) 21*1 cell, each containing 72*120 matrix of each participant's response.
%%

function subResponse = subject_response(subName)

cd(fileparts(mfilename('fullpath')));

nRun = 4; % Number of runs 
nTrial = 24; % Number of trials in one run
mainSeq = [1:18,25:42,49:66,73:90]; % Sequence number used for CD condition
nBin = 120; % Number of bins in one trial

if nargin == 1
    subjects = {subName};
else
    % Load the subject list used
    load('../subjectList.mat');
end

% Define the output size
subResponse = cell(length(subjects),1);

for iSub = 1:length(subjects)
    
    % seqOrder is loaded
    load(['../../data/beh/cdm' subjects{iSub} '/seqOrder.mat']);
    
    for iRun = 1:nRun
        
        % Load subject behavior
        data = load(['../../data/beh/cdm' subjects{iSub} ...
            '/cdm' subjects{iSub} '_run' num2str(iRun) '.mat']);
        
        for iTrial = 1:nTrial
            
            if data.cnd(iTrial,1) ~= 10 % if main trial
                if isempty(data.detectTime{iTrial})
                    response = zeros(1,nBin+1);
                else
                    response = histc(data.detectTime{iTrial} ...
                        - data.runStart(iTrial),0:0.25:30);
                end
                subResponse{iSub}(find(seqOrder(iRun,iTrial)==mainSeq),:) ...
                    = response(1,1:end-1);
            end
            
        end
        
    end
    
end

save('../subResponse.mat','subResponse');

end