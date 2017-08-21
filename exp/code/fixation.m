% i_Run: current trial index 
% tN: total number of the trial in one session

function fixation(T, RGB)

global window windowW windowH

Screen('glTranslate', window, windowW/2, windowH/2);
Screen('DrawLine', window, RGB, 0, -12.5, 0, 12.5, 3);
Screen('DrawLine', window, RGB, -12.5, 0, 12.5, 0, 3);
Screen('glTranslate', window, -windowW/2, -windowH/2);
Screen('Flip',window);

start = GetSecs;
while GetSecs-start<T; end

end
