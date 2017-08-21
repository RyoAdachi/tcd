%% TEMPORAL CHANGE DETECTION TASK FOR FMRI EXPERIMENT (2/15/2015)

% See task_run1 for variable descriptions

%%

function task_run2()

clear; clc;

rand('state',sum(100*clock));

% Subject ID (e.g. 001)
id = input('Input subject ID (e.g. 001): ','s');

global windowH windowW window

%% THESE SEQUENCES ARE USED

strName = [14,95]; % if id is odd number then 34, 45 when even
iRun = 2;

%%
strUsed = -rem(str2num(id),2)+2;

% Load the predefined sequence
dgp = load(['../raw_sequence/data_mri_main_20150205_' num2str(strName(strUsed)) '.mat']);
nTotal = dgp.nTotal; % Number of trials in one session (=24)
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

load(['../../data/beh/cdm' id '/seqOrder.mat']);

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
    [ keyIsDown, ~, keyCode ] = KbCheck([-1]);
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
