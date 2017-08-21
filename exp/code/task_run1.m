%% TEMPORAL CHANGE DETECTION TASK FOR FMRI EXPERIMENT (2/15/2015)

%%

% Below is the information of the variables saved.
% mriStart: start time of the fmri scanner
% mriEnd: end time of the fmri scanner

% For all variable below, rows index trials.
% nSeq: which of the 96 sequences of data_mri_main_20150205...mat is used in each trial
% cnd: Property of the sequence
%      First column is the identity of the trial (1:easy, 2:hard or 10:oddball)
%      Second column indicate the type of change (1: low to high, 2: high to low freq)
% stimBinary: existence (1) or nonexistence(0) of stimulus in 120 bins
% changeLocation: changepoint location in bin number (1-120)
% pValue: p (Bernoulli prob) of the sequence (e.g. [0.1,0.5]: change from low to high freq)
% oddballLocation: timing of oddball presentation in bin number
% detectTime: timing of button press
% runStart: start time of each trial
% runEnd: end time of each trial (should be about 30s later of runStart)
% stimStart: onset of each stimulus
% stimEnd: offset of each stimulus
% jitter: ITI information

%%

function task_run1()

clear; clc;

rand('state',sum(100*clock));

% Subject ID (e.g. 001)
id = input('Input subject ID (e.g. 001): ','s');

global windowH windowW window

%% THESE SEQUENCES ARE USED

strName = [14,95]; % if id is odd number then 34, 45 when even
iRun = 1;

%%
strUsed = -rem(str2num(id),2)+2;

% Load the predefined sequence
dgp = load(['../raw_sequence/data_mri_main_20150205_' num2str(strName(strUsed)) '.mat']);
nTotal = dgp.nTotal; % Number of trials in one session (=24)
nRun = dgp.nRun;
whiteRGB = [255 255 255];
greenRGB = [0 255 0];
blackRGB = [0 0 0];

% Enable unified mode of KbName, so KbName accepts identical key names on
% all operating systems:
priorityLevel=MaxPriority(['GetSecs'],['KbCheck'],['KbWait'],['GetClicks']);
KbName('UnifyKeyNames');
ListenChar(2)
HideCursor

% Open full screen (external display in the scanner)
whichScreen = 1;
window = Screen(whichScreen, 'OpenWindow');
[windowW, windowH] = Screen('WindowSize', window);
% Enable alpha blending
Screen('BlendFunction', window, GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA) ;

% Black background
Screen('FillRect', window, blackRGB);

% Randomize the order of sequence 
% (number of easy, hard, oddball conditions across runs)
seqOrder = transpose(reshape(1:nRun*nTotal,nTotal,nRun));
for iCol=1:nTotal
    seqOrder(:,iCol) = seqOrder(randperm(nRun),iCol);
end
for iRow=1:nRun
    seqOrder(iRow,:) = seqOrder(iRow,randperm(nTotal));
end

% Store the order of sequences in the four runs
save(['../../data/beh/cdm' id '/seqOrder.mat'],'seqOrder');

% Jitter matrix (uniform, mean=4)
jitterVec = 2:4/(nTotal-1):6;

%% Task start

% Preassign variables
runStart = zeros(nTotal,1); runEnd = zeros(nTotal,1); 
stimStart = cell(nTotal,1); stimEnd = cell(nTotal,1);
stimBinary = cell(nTotal,1); changeLocation = cell(nTotal,1); detectTime = cell(nTotal,1);
pValue = cell(nTotal,1); oddballLocation = cell(nTotal,1); nSeq = zeros(nTotal,1);
cnd = zeros(nTotal,2); 

for iTrial = 1:nTotal    
    nSeq(iTrial) = seqOrder(iRun,iTrial);
    stimBinary{iTrial} = dgp.S{nSeq(iTrial)};
    changeLocation{iTrial} = dgp.changeLocation{nSeq(iTrial)};
    pValue{iTrial} = dgp.pValue{nSeq(iTrial)};
    oddballLocation{iTrial} = dgp.oddballLocation{nSeq(iTrial)};
    cnd(iTrial,:) = transpose(dgp.cnd(:,nSeq(iTrial)));
end

jitter = jitterVec(1,randperm(nTotal)); % pseudorandom ITI (2-6s UNIF)

% Ready screen
Screen(window,'TextSize',30);
disp_text = ['Get ready to start session ' num2str(iRun) '.'];
DrawFormattedText(window, disp_text, 'center', 'center', whiteRGB);
clear disp_text;
Screen('Flip',window);

% Task initiated by trigger pulse ('5')
while 1
    [keyIsDown, ~, keyCode ] = KbCheck([-1]);
    if keyIsDown
        if strcmp(KbName(keyCode),'5%')
            break;
        end
        while KbCheck; end
    end
end

% Start time of the scanner
mriStart = GetSecs;

for iTrial = 1:nTotal
    
    switch cnd(iTrial,1)
        
        case {10} % oddball condition
            
            % Fixation
            fixation(jitter(iTrial), greenRGB);
            
            % Image presentation
            [runStart(iTrial), runEnd(iTrial), detectTime{iTrial}, stimStart{iTrial}, stimEnd{iTrial}] = ...
                oddball_temporal_mri(stimBinary{iTrial}, oddballLocation{iTrial});
            
        case {1,2} % main condition
            
            % Fixation
            fixation(jitter(iTrial), whiteRGB);
            
            % Image presentation
            [runStart(iTrial), runEnd(iTrial), detectTime{iTrial}, stimStart{iTrial}, stimEnd{iTrial}] = ...
                temporal_change_mri(stimBinary{iTrial});
            
    end
    
end

% End time of the scanner
mriEnd = GetSecs;

% End of run screen
Screen(window,'TextSize',30);
disp_text = ['End of session ' num2str(iRun) '.'];
DrawFormattedText(window, disp_text, 'center', 'center', whiteRGB);
clear disp_text;
Screen('Flip',window);

save (['../../data/beh/cdm' id '/cdm' id '_run' num2str(iRun)],'mriStart','mriEnd','nSeq','cnd','stimBinary',...
    'changeLocation','pValue','oddballLocation','detectTime','runStart','runEnd',...
    'stimStart','stimEnd','jitter');

load('gong.mat'); sound(y);

% Exit 
ListenChar(0);
ShowCursor;
Screen(window,'Close');

end
