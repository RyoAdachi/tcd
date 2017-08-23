%%
% This function evaluate the model's performance by calculating the 
% f measure between the model's response and participants' response.

% The inputs are the following.
% (1) model's response (2*1 cell each containing a 72*120 matrix of 1/0).
% (2) a model's BP within distanceTime[s] apart from actual subject
% response is a TP.
% (3) subject identifier (e.g. '002'). This is used in random model.
% For other Bayesian models, leave it blank.

%%

function [fIndividual, fIndividualEasy, fIndividualHard, ...
    precision, precisionEasy, precisionHard, ...
    sensitive, sensitiveEasy, sensitiveHard, ...
    fmeasureEven, fmeasureEvenEasy, fmeasureEvenHard, ...
    precisionEven, precisionEvenEasy, precisionEvenHard, ...
    sensitiveEven, sensitiveEvenEasy, sensitiveEvenHard, fIndividualCrossValidationH] = ...
    model_subject_response_evaluator(response, distanceTime, subName)

cd(fileparts(mfilename('fullpath')));

% Load the subjet list used
load('../subjectList.mat');

% Number of sequence in the main condition
nSeq = 72; 

% Get subject response matrix (21 cell, each contains 72*120 matrix)
load('../subResponse.mat');

% Load difficulty vector
load('../diffVec.mat');

% This is used in the random model simulation
if nargin > 2
    subResponse = {subResponse{find(strcmp(subjects,subName))}};
    subjects = {subName};
end

% Define the vectors
fIndividual = zeros(length(subjects),1);
fIndividualEasy = zeros(length(subjects),1);
fIndividualHard = zeros(length(subjects),1);
fIndividualCrossValidationH = zeros(length(subjects),2);

for iSub = 1:length(subjects)
        
    % Get the sequence order for the subject
    load(['../../data/beh/cdm' subjects{iSub} '/seqOrder.mat']);
    
    % Initialize (store TP, FP, FN for each block)
    tpAgg = zeros(4,2); fnAgg = zeros(4,2); fpAgg = zeros(4,2);
    tpAggOddEven = zeros(2,2); fnAggOddEven = zeros(2,2); fpAggOddEven = zeros(2,2);
    
    % Which of the two sequences did the subject use?
    % If 1 then seq14 was used; if 2 then seq95 was used
    seqIndex = -rem(str2num(subjects{iSub}),2)+2;
    
    for iSeq = 1:nSeq
        
        % Is this in the odd/even trial?
        [~, trialType] = find(reshape(seqOrder',2,[])'==iSeq);
        
        % Which block does the sequene belong to?
        [blockIndex, ~] = find(seqOrder==iSeq);
        
        % Location of the subject and model button press
        % by the bin number
        nSubBp = find(subResponse{iSub}(iSeq,:) == 1);
        if nargin > 2
            nModelBp = find(response{1}(iSeq,:) == 1);
        else
            nModelBp = find(response{seqIndex}(iSeq,:) == 1);
        end
        
        % Classify a models button press to TP, FP and get FN
        TP = 0; tmpModelBp = nModelBp;
        if ~isempty(nSubBp)
            for iSubBp = 1:length(nSubBp)
                tmp = tmpModelBp - nSubBp(iSubBp);
                neg_min = -min(abs(tmp(tmp<0)));
                pos_min = min(tmp(tmp>=0));
                if ~isempty(neg_min) && abs(neg_min) <= distanceTime/0.25
                    tmpModelBp(find(tmp==neg_min)) = [];
                    TP = TP + 1;
                elseif ~isempty(pos_min) && pos_min <= distanceTime/0.25
                    tmpModelBp(find(tmp==pos_min)) = [];
                    TP = TP + 1;
                end
            end
        end
        
        FN = length(nSubBp) - TP;
        FP = length(nModelBp) - TP;
        
        cnd = diffVec{seqIndex}(iSeq);
        
        % Aggeregate the TP, FN, FP
        tpAgg(blockIndex,cnd) = tpAgg(blockIndex,cnd) + TP;
        fnAgg(blockIndex,cnd) = fnAgg(blockIndex,cnd) + FN;
        fpAgg(blockIndex,cnd) = fpAgg(blockIndex,cnd) + FP;
        
        tpAggOddEven(trialType,cnd) = tpAggOddEven(trialType,cnd) + TP;
        fnAggOddEven(trialType,cnd) = fnAggOddEven(trialType,cnd) + FN;
        fpAggOddEven(trialType,cnd) = fpAggOddEven(trialType,cnd) + FP;
        
    end
    
    precision = sum(sum(tpAgg))/sum(sum(tpAgg+fpAgg));
    sensitive = sum(sum(tpAgg))/sum(sum(tpAgg+fnAgg));
    precisionEasy = sum(tpAgg(:,1))/sum(tpAgg(:,1)+fpAgg(:,1));
    sensitiveEasy = sum(tpAgg(:,1))/sum(tpAgg(:,1)+fnAgg(:,1));
    precisionHard = sum(tpAgg(:,2))/sum(tpAgg(:,2)+fpAgg(:,2));
    sensitiveHard = sum(tpAgg(:,2))/sum(tpAgg(:,2)+fnAgg(:,2));
    fIndividual(iSub,1) = 2/(1/precision+1/sensitive);
    fIndividualEasy(iSub,1) = 2/(1/precisionEasy+1/sensitiveEasy);
    fIndividualHard(iSub,1) = 2/(1/precisionHard+1/sensitiveHard);
    
    precisionOddEven = sum(tpAggOddEven,2)./sum((tpAggOddEven+fpAggOddEven),2);
    sensitiveOddEven = sum(tpAggOddEven,2)./sum((tpAggOddEven+fnAggOddEven),2);
    fIndividualCrossValidationH(iSub,:) = 2./(1./precisionOddEven+1./sensitiveOddEven);
    
    precisionEven = sum(tpAggOddEven(2,:))/sum(tpAggOddEven(2,:)+fpAggOddEven(2,:));
    precisionEvenEasy = tpAggOddEven(2,1)/(tpAggOddEven(2,1)+fpAggOddEven(2,1));
    precisionEvenHard = tpAggOddEven(2,2)/(tpAggOddEven(2,2)+fpAggOddEven(2,2));
    
    sensitiveEven = sum(tpAggOddEven(2,:))/sum(tpAggOddEven(2,:)+fnAggOddEven(2,:));
    sensitiveEvenEasy = tpAggOddEven(2,1)/(tpAggOddEven(2,1)+fnAggOddEven(2,1));
    sensitiveEvenHard = tpAggOddEven(2,2)/(tpAggOddEven(2,2)+fnAggOddEven(2,2));
    
    fmeasureEven = 2./(1./precisionEven+1./sensitiveEven);
    fmeasureEvenEasy = 2./(1./precisionEvenEasy+1./sensitiveEvenEasy);
    fmeasureEvenHard = 2./(1./precisionEvenHard+1./sensitiveEvenHard);
    
end

end
