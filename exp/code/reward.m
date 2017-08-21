%% REWARD PAYMENT FOR THE TEMPORAL CHANGE DETECTION FMRI TASK (2/21/2015)

% chosenTrial: whichi trial is chosen for each run
% reward: how much is earned in each run
% rewardTotal: how much in total is earned (sum(reward))

%%

function reward()

clear; clc;

rand('state',sum(100*clock));

% subject ID (e.g. 001);
id = input('Input subject ID (e.g. 001): ','s');

nRun = 4; % number of runs;
nTrial = 24; % number of trials in one run
mainDuration = 5; % a BP within this time duration is regarded as a correct change detection
oddDuration = 1; % a BP within this time duration is regarded as a correct odd ball detection
sec = 1/8; % the duration of the image presentation
nTrialUsed = 6;

whiteRGB = [255 255 255];

% Open full screen
whichScreen = 0;
window = Screen(whichScreen, 'OpenWindow');

% Black background
black = BlackIndex(window); % pixelvalue for black
Screen('FillRect', window, black);

% Show the result
Screen(window,'TextSize',60);
disp_text = ['Press any key to start.'];
DrawFormattedText(window, disp_text, 'center', 'center', whiteRGB);
clear disp_text;
Screen('Flip',window);

KbWait([-1],2);

reward = zeros(nRun,1);
for iRun = 1:nRun
    
    data = load(['../../data/beh/cdm' id '/cdm' id '_run' num2str(iRun)]);
    selectedTrial = datasample([1:1:nTrial],nTrialUsed,'Replace',false);
    
    for iTrialUsed = 1:nTrialUsed
        
        iTrial = selectedTrial(iTrialUsed);
        bpressTime = data.detectTime{iTrial} - data.runStart(iTrial);
        
        TP = 0; FP = 0; FN = 0;
        
        switch data.cnd(iTrial)
            
            case {1,2} % main condition
                
                nChange = length(data.changeLocation{iTrial});
                
                if ~(nChange == 0)
                    for j=1:nChange
                        tmp = bpressTime - (sec + 2*sec*(data.changeLocation{iTrial}(j)-1));
                        if sum(tmp>=0)>0
                            TP = TP + (min(tmp(tmp>=0)) <= mainDuration);
                        end
                    end
                end
                
                FN = nChange - TP;
                FP = length(bpressTime) - TP;
                
            case {10}
                
                nOddball = length(data.oddballLocation{iTrial});
                
                if ~(nOddball == 0)
                    for j=1:nOddball
                        tmp = bpressTime - (sec + 2*sec*(data.oddballLocation{iTrial}(j)-1));
                        if sum(tmp>=0)>0
                            TP = TP + (min(tmp(tmp>=0)) <= oddDuration);
                        end
                    end
                end
                
                FN = length(data.oddballLocation{iTrial}) - TP;
                FP = length(bpressTime) - TP;
                
        end
        
        reward(iRun) = reward(iRun) + 2*(TP-FP);
        
    end
    
end

rewardTotal = min(max(5,sum(reward)),40);

% Show the result
Screen(window,'TextSize',60);
disp_text = ['The computer is calculating the payout...'];
DrawFormattedText(window, disp_text, 'center', 'center', whiteRGB);
clear disp_text;
Screen('Flip',window);

wait(3);

% Show the result
Screen(window,'TextSize',60);
disp_text = ['You earned $' num2str(rewardTotal+10) ' for your performance.'];
DrawFormattedText(window, disp_text, 'center', 'center', whiteRGB);
clear disp_text;
Screen('Flip',window);

wait(3);

% Show the result
Screen(window,'TextSize',60);
disp_text = ['You get $' num2str(rewardTotal+10+35) ' in total.'];
DrawFormattedText(window, disp_text, 'center', 'center', whiteRGB);
clear disp_text;
Screen('Flip',window);

save(['../payment/cdm' id], 'selectedTrial', 'reward', 'rewardTotal');

% Exit by any key pressing
KbWait([-1],2);

ListenChar(0);
ShowCursor;
Screen(window,'Close');

    function wait(T)
        t = GetSecs;
        while GetSecs - t < T
        end
    end

end