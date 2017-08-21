%% TEMPORAL CHANGE DETECTION TASK FOR BEHAVIORAL EXPERIMENT PRACTICE (3/19/2015)

% In this version, the stimuli is Bernoulli and changepoints are Bernoulli

%%

clear; clc;

rand('state',sum(100*clock));

global windowH windowW window

% Input the subject trial information
id = input('Input subject ID (e.g. 001): ','s'); % subject ID (e.g. 001);

% Create folders to save 
mkdir(['../../data/beh/cdm' id ]);

% Load the predefined sequence
dgp = load(['../raw_sequence/data_mri_prac_20150205.mat']);
nRun = dgp.nRun; nTotal = dgp.nTotal;
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

% Randomize the predefined sequence
seqOrder = transpose(reshape(1:nRun*nTotal,nTotal,nRun));
for i=1:nTotal
    seqOrder(:,i) = seqOrder(randperm(nRun),i);
end
for i=1:nRun
    seqOrder(i,:) = seqOrder(i,randperm(nTotal));
end
iRun = 1;

%% Practice start

jitter = 2; % no need to randomize in behavioral exp

% Ready screen
Screen(window,'TextSize',30);
disp_text = 'Press any key to begin practice.';
DrawFormattedText(window, disp_text, 'center', 'center', whiteRGB);
clear disp_text;
Screen('Flip',window);

% Task initiated by any key pressing
KbWait([-1],2);

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
    
    switch cnd(iTrial,1)
        
        case {10} % control condition
            
            % Fixation
            fixation(jitter, greenRGB);
            
            % Image presentation
            [runStart(iTrial), runEnd(iTrial), detectTime{iTrial}, stimStart{iTrial}, stimEnd{iTrial}] = ...
                oddball_temporal_mri(stimBinary{iTrial}, oddballLocation{iTrial});
            
        case {1,2} % main condition
            
            % Fixation
            fixation(jitter, whiteRGB);
            
            % Image presentation
            [runStart(iTrial), runEnd(iTrial), detectTime{iTrial}, stimStart{iTrial}, stimEnd{iTrial}] = ...
                temporal_change_mri(stimBinary{iTrial});
            
    end
    
end

% End of run screen
Screen(window,'TextSize',30);
disp_text = 'End of practice. Press any key to progress.';
DrawFormattedText(window, disp_text, 'center', 'center', whiteRGB);
clear disp_text;
Screen('Flip',window);

KbWait([-1],2);

save (['../../data/beh/cdm' id '/practice.mat'],'nSeq','cnd','stimBinary','changeLocation','pValue','oddballLocation','detectTime','runStart','runEnd','stimStart','stimEnd');

load('gong.mat'); sound(y)

%% Experiment complete screen
Screen(window,'TextSize',30);
disp_text = 'Waiting for anatomical scan to complete...';
DrawFormattedText(window, disp_text, 'center', 'center', whiteRGB);
clear disp_text;
Screen('Flip',window);

%% Exit by any key pressing
KbWait;
while KbCheck; end;
ListenChar(0);
ShowCursor;
Screen(window,'Close');
fprintf('\nDone.\n');
