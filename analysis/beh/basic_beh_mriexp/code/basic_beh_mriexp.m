%% THIS SCRIPT COMPUTES THE BASIC STATISTICS OF BEHAVIORAL DATA OF FMRI SUBJECTS

% Conditions [1,2]: easy, hard in change detection condition, 
% [10]: oddball condition
% For no change: pValue = 0.1,0.15 are easy and pValue=0.45,0.5 are hard
% For oddball condition, pValue = 0.3;

%% VARIABLES STORED
% tpMainAggregate: tp for main condition across all trials (row: nChange, col1: easy, col2: hard)
% tpOddAggregate: tp for odd
% fpMainAggregate: fp for main
% fpOddAggregate: fp for odd
% fnMainAggregate: fn for main
% fnOddAggregate: fn for odd
% cndMainAggregate: condition count for main
% cndOddAggregate: condition count for odd
% tpMainSession: tp for main condition for each session
% fpMainSession: fp for main condition for each session
% fnMainSession: fn for main condition for each session
% tpOddSession: tp for odd condition for each session
% fpOddSession: fp for odd condition for each session
% fnOddSession: fn for odd condition for each session
% nMainChangeSession: Number of changes for each main session
% nMainBpressSession: Number of button press for each main session
% nOddSession: Number of changes for each odd session
% nOddBpressSession: Number of button press for each odd session
% changeHL: Number of changes in main condition for L->H and H->L change
% vecTP: TP for each sequence ordered in the original sequence order
% vecFP: FP for each sequence ordered in the original sequence order
% vecFN: FN for each sequence ordered in the original sequence order
% fMainSession: : f measure for main condition for each session
% precisionMainSession: precision for main condition for each session
% sensitiveMainSession: sensitivity for main condition for each session
% mainBpress: Number of button press for main condition per trial
% precisionMain: precision for main condition across all trials
% precisionMainDifficulty: precision for main condition for each difficulty
% precisionHL: precision for main condition for L->H and H->L change
% sensitiveMain: sensitivity for main condition across all trials
% sensitiveMainDifficulty: sensitivity for main condition for each difficulty
% sensitiveHL: sensitivity for main condition for L->H and H->L change
% fMain: f measure for main condition across all trials
% fMainDifficulty: f measure for main condition for each difficulty
% fHL: f measure for main condition for L->H and H->L change
% oddBpress: Number of button press for oddball condition per trial
% precisionOdd: precision for oddball condition across all trials
% sensitiveOdd: sensitivity for oddball condition across all trials
% fOdd: f measure for oddball condition across all trials
% mainBpressDifficulty: Number of button press for main condition per trial for each difficulty
% nChangeDifficulty: number of change for each difficulty
% nOdd: number of oddballs 
% fOddSession: f measure for contol condition for each session
% precisionOddSession: precision for contol condition for each session
% sensitiveOddSession: sensitivity for contol condition for each session
% easyPrecisionMainSession: precision for easy condition for each session
% easySensitiveMainSession: sensitivity for easy condition for each session
% easyfMainSession: fmeasure for easy condition for each session
% hardPrecisionMainSession: precision for hard condition for each session
% hardSensitiveMainSession: sensitivity for hard condition for each session
% hardfMainSession: fmeasure for hard condition for each session
% precisionMainOddEven: precision for main condition for odd/even numbered trials
% sensitiveMainOddEven: sensitivity for main condition for odd/even numbered trials
% fMainOddEven: fmeasure for main condition for odd/even numbered trials
% easyPrecisionMainOddEven: precision for easy condition for odd/even numbered trials
% easySensitiveMainOddEven: sensitivity for easy condition for odd/even numbered trials
% easyfMainOddEven: fmeasure for easy condition for odd/even numbered trials
% hardPrecisionMainOddEven: precision for hard condition for odd/even numbered trials
% hardSensitiveMainOddEven: sensitivity for hard condition for odd/even numbered trials
% hardfMainOddEven: fmeasure for hard condition for odd/even numbered trials
%%

clear; clc; close all;

% Session information
nRun = 4; % Number of runs;
nMain = 18; % Number of main trials in one run
nOdd = 6; % Number of oddball trials in one run
nTrial = nMain + nOdd; % Total number of trials in one run
mainDuration = 5; % A BP within this time duration is regarded as a correct change detection
oddDuration = 1; % A BP within this time duration is regarded as a correct oddball detection
secImage = 1/8; % Duration of one image presentation

% Which numbered sequences are used for the main task
mainSeqNum = [1:18,25:42,49:66,73:90];

% Load the subject IDs used for the analysis
load('../../../subjectList.mat');

for iSub = 1:length(subjects)
    
    tpMainAggregate = zeros(4,2); fpMainAggregate = zeros(4,2); fnMainAggregate = zeros(4,2);
    tpOddAggregate = zeros(4,1); fpOddAggregate = zeros(4,1); fnOddAggregate = zeros(4,1);
    cndMainAggregate = zeros(4,2); cndOddAggregate = zeros(4,1);
    
    tpMainSession = zeros(4,2); fpMainSession = zeros(4,2); fnMainSession = zeros(4,2);    
    tpMainOddEven = zeros(2,2); fpMainOddEven = zeros(2,2); fnMainOddEven = zeros(2,2);
    tpOddSession = zeros(4,1); fpOddSession = zeros(4,1); fnOddSession = zeros(4,1);
    nMainChangeSession = zeros(4,1); nMainBpressSession = zeros(4,1);
    nOddSession = zeros(4,1); nOddBpressSession = zeros(4,1);
    
    changeHL = zeros(3,2); % 1col:LtoH, 2col:HtoL, 1row:actual, 2row:TP, 3row:FP
    
    vecTP = zeros(72,1); vecFN = zeros(72,1); vecFP = zeros(72,1);
    
    for iRun = 1:nRun
        
        data = load(['../../../../data/beh/cdm' subjects{iSub} '/cdm' subjects{iSub} '_run' num2str(iRun) '.mat']);
        
        for iTrial = 1:nTrial
            
            TP = 0; FN = 0; FP = 0;
            bpress = data.detectTime{iTrial} - data.runStart(iTrial); % Timing of button press
            
            switch data.cnd(iTrial,1)
                
                case {1,2} % Main condition
                    
                    nChange = length(data.changeLocation{iTrial}); % Number of actual changes
                    
                    if ~(nChange == 0)
                        diff_pValue = diff(data.pValue{iTrial}); % Difference in the p of Bernoulli trials
                        for iChange = 1:nChange
                            tmp = bpress - (secImage + 2*secImage*(data.changeLocation{iTrial}(iChange) - 1));
                            if sum(tmp >= 0) > 0
                                TP = TP + (min(tmp(tmp >= 0)) <= mainDuration);
                            end
                        end
                    end
                    
                    FN = nChange - TP;
                    FP = length(bpress) - TP;
                    
                    if nChange == 1
                        iCol = -0.5*sign(diff_pValue)+1.5; % Low to High or vice versa
                        % L to H is 1Col, H to L is 2Col
                        changeHL(:,iCol) = changeHL(:,iCol) + [1; TP; FP];
                    end
                    
                    clear iCol
                    
                    index = find(data.nSeq(iTrial) == mainSeqNum);
                    vecTP(index) = TP; vecFN(index) = FN; vecFP(index) = FP;
                    
                    %%
                    
                    iRow = nChange + 1;
                    
                    if nChange == 0
                        switch data.pValue{iTrial}
                            case {0.1,0.15}
                                iCol = 1;
                            case {0.45,0.5}
                                iCol = 2;
                        end
                    else
                        iCol = data.cnd(iTrial,1);
                    end
                    
                    % Result for the total sessions aggregated
                    % Summarize main condition (row:change number, col:difficulty)
                    tpMainAggregate(iRow,iCol) = tpMainAggregate(iRow,iCol) + TP;
                    fpMainAggregate(iRow,iCol) = fpMainAggregate(iRow,iCol) + FP;
                    fnMainAggregate(iRow,iCol) = fnMainAggregate(iRow,iCol) + FN;
                    cndMainAggregate(iRow,iCol) = cndMainAggregate(iRow,iCol) + 1;
                    
                    % Result for each session
                    tpMainSession(iRun,iCol) = tpMainSession(iRun,iCol) + TP;
                    fpMainSession(iRun,iCol) = fpMainSession(iRun,iCol) + FP;
                    fnMainSession(iRun,iCol) = fnMainSession(iRun,iCol) + FN;
                    
                    % Result for odd/even trials
                    index = 2 - mod(iTrial,2);
                    tpMainOddEven(index,iCol) = tpMainOddEven(index,iCol) + TP;
                    fpMainOddEven(index,iCol) = fpMainOddEven(index,iCol) + FP;
                    fnMainOddEven(index,iCol) = fnMainOddEven(index,iCol) + FN;
                    
                    % Check for any differences between sessions
                    nMainChangeSession(iRun) = nMainChangeSession(iRun) + nChange;
                    nMainBpressSession(iRun) = nMainBpressSession(iRun) + length(bpress);
                    
                case {10} % Control condition
                    
                    nOdd = length(data.oddballLocation{iTrial}); % Number of actual oddballs
                    
                    if ~(nOdd == 0)
                        for iOdd = 1:nOdd
                            tmp = bpress - (secImage + 2*secImage*(data.oddballLocation{iTrial}(iOdd) - 1));
                            if sum(tmp >= 0) > 0
                                TP = TP + (min(tmp(tmp >= 0)) <= oddDuration);
                            end
                        end
                    end
                    
                    FN = nOdd - TP;
                    FP = length(bpress) - TP;
                    
                    %%
                    
                    iRow = nOdd + 1;
                    iCol = 1;
                    
                    % Result for overall sessions
                    % Summarize oddball condition (row:oddball number, col:difficulty)
                    tpOddAggregate(iRow,iCol) = tpOddAggregate(iRow,iCol) + TP;
                    fpOddAggregate(iRow,iCol) = fpOddAggregate(iRow,iCol) + FP;
                    fnOddAggregate(iRow,iCol) = fnOddAggregate(iRow,iCol) + FN;
                    cndOddAggregate(iRow,iCol) = cndOddAggregate(iRow,iCol) + 1;
                    
                    % Result for each session
                    tpOddSession(iRun) = tpOddSession(iRun) + TP;
                    fpOddSession(iRun) = fpOddSession(iRun) + FP;
                    fnOddSession(iRun) = fnOddSession(iRun) + FN;
                    
                    % Check for any differences between sessions
                    nOddSession(iRun) = nOddSession(iRun) + nOdd;
                    nOddBpressSession(iRun) = nOddBpressSession(iRun) + length(bpress);
                    
            end
            
        end
        
    end
    
    % Compute other relevant variables (main condition)
    mainBpress = sum(sum(tpMainAggregate+fpMainAggregate))/sum(sum(cndMainAggregate));
    mainBpressDifficulty = sum(tpMainAggregate+fpMainAggregate,1)./sum(cndMainAggregate,1);
    nChangeDifficulty = sum(tpMainAggregate+fnMainAggregate,1)./sum(cndMainAggregate,1);
    precisionMain = sum(sum(tpMainAggregate))/sum(sum(tpMainAggregate + fpMainAggregate));
    precisionMainDifficulty = sum(tpMainAggregate,1)./sum(tpMainAggregate+fpMainAggregate,1);
    precisionHL = changeHL(2,:)./sum(changeHL(2:3,:),1);
    sensitiveMain = sum(sum(tpMainAggregate))/sum(sum(tpMainAggregate + fnMainAggregate));
    sensitiveMainDifficulty = sum(tpMainAggregate,1)./sum(tpMainAggregate+fnMainAggregate,1);
    sensitiveHL = changeHL(2,:)./changeHL(1,:);
    fMain = 2/(1/precisionMain + 1/sensitiveMain);
    fMainDifficulty = 2./(1./precisionMainDifficulty + 1./sensitiveMainDifficulty);
    fHL = 2./(1./precisionHL + 1./sensitiveHL);
    precisionMainSession = sum(tpMainSession,2)./sum((tpMainSession + fpMainSession),2);
    sensitiveMainSession = sum(tpMainSession,2)./sum((tpMainSession + fnMainSession),2);
    fMainSession = 2./(1./precisionMainSession + 1./sensitiveMainSession);
    easyPrecisionMainSession = tpMainSession(:,1)./(tpMainSession(:,1) + fpMainSession(:,1));
    easySensitiveMainSession = tpMainSession(:,1)./(tpMainSession(:,1) + fnMainSession(:,1));
    easyfMainSession = 2./(1./easyPrecisionMainSession + 1./easySensitiveMainSession);
    hardPrecisionMainSession = tpMainSession(:,2)./(tpMainSession(:,2) + fpMainSession(:,2));
    hardSensitiveMainSession = tpMainSession(:,2)./(tpMainSession(:,2) + fnMainSession(:,2));
    hardfMainSession = 2./(1./hardPrecisionMainSession + 1./hardSensitiveMainSession);
    
    precisionMainOddEven = sum(tpMainOddEven,2)./sum((tpMainOddEven + fpMainOddEven),2);
    sensitiveMainOddEven = sum(tpMainOddEven,2)./sum((tpMainOddEven + fnMainOddEven),2);
    fMainOddEven = 2./(1./precisionMainOddEven + 1./sensitiveMainOddEven);
    easyPrecisionMainOddEven = tpMainOddEven(:,1)./(tpMainOddEven(:,1) + fpMainOddEven(:,1));
    easySensitiveMainOddEven = tpMainOddEven(:,1)./(tpMainOddEven(:,1) + fnMainOddEven(:,1));
    easyfMainOddEven = 2./(1./easyPrecisionMainOddEven + 1./easySensitiveMainOddEven);
    hardPrecisionMainOddEven = tpMainOddEven(:,2)./(tpMainOddEven(:,2) + fpMainOddEven(:,2));
    hardSensitiveMainOddEven = tpMainOddEven(:,2)./(tpMainOddEven(:,2) + fnMainOddEven(:,2));
    hardfMainOddEven = 2./(1./hardPrecisionMainOddEven + 1./hardSensitiveMainOddEven);    
    
    % Compute other relevant variables (oddball condition)
    oddBpress = sum(tpOddAggregate + fpOddAggregate)/sum(cndOddAggregate);
    precisionOdd = sum(tpOddAggregate)/sum(tpOddAggregate + fpOddAggregate);
    sensitiveOdd = sum(tpOddAggregate)/sum(tpOddAggregate + fnOddAggregate);
    fOdd = 2/(1/precisionOdd + 1/sensitiveOdd);
    nOdd = sum(tpOddAggregate + fnOddAggregate)/sum(cndOddAggregate);
    precisionOddSession = tpOddSession./(tpOddSession + fpOddSession);
    sensitiveOddSession = tpOddSession./(tpOddSession + fnOddSession);
    fOddSession = 2./(1./precisionOddSession + 1./sensitiveOddSession);
    
    save(['../result/cdm' subjects{iSub} '.mat'],...
        'tpMainAggregate','tpOddAggregate','fpMainAggregate','fpOddAggregate','fnMainAggregate',...
        'fnOddAggregate','cndMainAggregate','cndOddAggregate','tpMainSession','fpMainSession',...
        'fnMainSession','tpOddSession','fpOddSession','fnOddSession','nMainChangeSession',...
        'nMainBpressSession','nOddSession','nOddBpressSession','changeHL','vecTP','vecFP','vecFN',...
        'fMainSession','precisionMainSession','sensitiveMainSession','mainBpress','precisionMain',...
        'precisionMainDifficulty','precisionHL','sensitiveMain','sensitiveMainDifficulty',...
        'sensitiveHL','fMain','fMainDifficulty','fHL','oddBpress','precisionOdd','sensitiveOdd',...
        'fOdd','mainBpressDifficulty','nChangeDifficulty','nOdd','fOddSession','precisionOddSession',...
        'sensitiveOddSession','easyPrecisionMainSession','easySensitiveMainSession','easyfMainSession',...
        'hardPrecisionMainSession','hardSensitiveMainSession','hardfMainSession',...
        'precisionMainOddEven','sensitiveMainOddEven','fMainOddEven',...
        'easyPrecisionMainOddEven','easySensitiveMainOddEven','easyfMainOddEven',...
        'hardPrecisionMainOddEven','hardSensitiveMainOddEven','hardfMainOddEven');
    
end

