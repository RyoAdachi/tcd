%% 
% THIS SCRIPT CHECKS THE HEAD MOVEMENT USING THE RP.TXT CREATED IN REALIGNMENT
% Input is the subject id (e.g. '002')
%%

function motion_check(id)

if nargin == 0
    error('Use argument of subject ID (e.g. 002)');
end

nRun = 4; % 4 runs in total

figure;
for iRun = 1:nRun
    
    movFileLoc = dir(['../data/scan/' id '/rsBOLD_MB_1_000' num2str(iRun+5) '/rp*.txt']);
    movFile = spm_load(['../data/scan/' id '/rsBOLD_MB_1_000' num2str(iRun+5) ...
        '/' movFileLoc(1).name]); % Select the rp*.txt file
    
    subplot(4,2,2*iRun-1);plot(movFile(:,1:3));
    set(gca,'xlim',[0 size(movFile,1)+1]);
    subplot(4,2,2*iRun);plot(movFile(:,4:6).*(180/pi)); % Convert from radian to degree
    set(gca,'xlim',[0 size(movFile,1)+1]);
    
    % Display maximum translation and rotation
    max(max(abs(movFile(:,1:3))))
    max(max(abs(movFile(:,4:6))))*(180/pi)
    
end

end