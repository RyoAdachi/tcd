% Odd ball trials for temporal change detection
% N is the condition number
% T is the total duration of one trial

function [runStart, runEnd, detect, stimStart, stimEnd] = oddball_temporal_mri(S, oddballTime)

global windowH windowW window

hSIZE = 70; % horizontal length of the image shown
vSIZE = hSIZE * 7/18; % vertical length of the image shown
sec = 1/8; % how long one image should stay on the screen

% Specify the rectangle where the image is drawn
[im1,~,alpha1] = imread('../../../image/sushi1.png');
im1(:,:,4) = alpha1(:,:);
t1Ptr = Screen('MakeTexture', window, im1);
t1Rect = Screen('Rect',t1Ptr);

[im2,~,alpha2] = imread(['../../../image/sushi1_54.png']); % oddball image
im2(:,:,4) = alpha2(:,:);
t2Ptr = Screen('MakeTexture', window, im2);
t2Rect = Screen('Rect',t2Ptr);
clear alpha1 alpha2 map1 map2;

runStart = GetSecs; 

n = 0; prekey = 0; detect = []; s_flag = zeros(1,sum(S)); 
stimStart = zeros(1,sum(S)); stimEnd = zeros(1,sum(S)); 

% Presentation of a sequence of pictures
while GetSecs - runStart <= 30
    
    [keyIsDown,timeSecs,keyCode] = KbCheck([-1]);
    if keyIsDown && ~prekey
        if strcmp(KbName(keyCode),'2@')
            n = n+1;
            detect(n) = timeSecs;
        end
        %while KbCheck; end;
    end
    prekey = keyIsDown;
    
    indexStimOn = find(S == 1);
    [dist, I] = min(abs(sec+2*sec*(indexStimOn-1) - (GetSecs - runStart)));
    if dist < sec/2
        if s_flag(I) == 0
            s_flag(I) = 1; stimStart(I) = GetSecs;
        end
        d1Rect = [windowW/2-hSIZE windowH/2-vSIZE windowW/2+hSIZE windowH/2+vSIZE];
        if indexStimOn(I) == oddballTime
            Screen('DrawTexture', window, t2Ptr, t2Rect, d1Rect);
        else
            Screen('DrawTexture', window, t1Ptr, t1Rect, d1Rect);
        end
        if s_flag(I) == 1
            stimEnd(I) = GetSecs;
        end
    end
    Screen('Flip',window);
    
end

runEnd = GetSecs;

end