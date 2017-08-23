%%
% Given the parameter set, calculate the response of the heuristic1 model.
% The inputs are the followings:
% (1) parameter_space: parameter vector
% (2) fit_name: e.g. main7
% (3) method_name: method name
% (4) subject_name: 'group' or subject name (e.g. '002')

% Output is the following.
% (1) 2*1 or 1*1 cell each contains 72*120 matrix of binary response for
% the sequence tested.
%%

function response = resp_heuristic1(parameter_space, fit_name, method_name, subject_name)

rp = '/home/radachi/research/CD/';

X = parameter_space(1); % Duration to consider (-X~0 seconds)
Y = parameter_space(2); % Duration to consider (-X-Y~-X seconds)
if ismember(fit_name,{'main1','main9','main13','main14','main15','main18','main20','main21'})
    theta_low_to_high_fix = parameter_space(3); % threshold to elicit button press (L to H)
    theta_high_to_low_fix = parameter_space(4); % threshold to elicit button press (H to L)
elseif ismember(fit_name,{'main2','main7','main8','main10','main16','main17','main19','main22'})
    theta_low_to_high_fix = parameter_space(3); % threshold to elicit button press (L to H)
    theta_high_to_low_fix = -parameter_space(3); % threshold to elicit button press (H to L)
end

% k
if ismember(fit_name,{'main7','main10','main17','main19'})
    k = parameter_space(4);
elseif ismember(fit_name,{'main9','main14','main18','main21'})
    k = parameter_space(5);
else
    k = 0;
end

% alpha, beta, delta
if ismember(fit_name,{'main8','main22'})
    alpha = parameter_space(4);
    beta = parameter_space(5);
    delta = parameter_space(6);
elseif ismember(fit_name,{'main10','main13','main19','main20'})
    alpha = parameter_space(5);
    beta = parameter_space(6);
    delta = parameter_space(7);
elseif ismember(fit_name,{'main14','main21'})
    alpha = parameter_space(6);
    beta = parameter_space(7);
    delta = parameter_space(8);
end

% T
if strcmp(fit_name,{'main16'})
    T = parameter_space(4);
elseif ismember(fit_name,{'main15','main17'})
    T = parameter_space(5);
elseif strcmp(fit_name,{'main18'})
    T = parameter_space(6);
end

number_of_bins = 120; % Number of bins
sequence_name = [14,95]; % Raw sequence name used in the experiment
main_sequence_index = [1:18,25:42,49:66,73:90]; % Sequence number used in main
number_of_trials = length(main_sequence_index);

% If subject ID is supplied in subject_name then response is a one by one cell
% and return result only for the sequence used by the particular subject
if ~strcmp(subject_name,'group')
    sequence_name = sequence_name(-rem(str2num(subject_name),2)+2);
    response = cell(1,1);
else
    response = cell(2,1);
end

%%

for i = 1:length(sequence_name)
    
    % Load the raw data sequence
    raw_sequence = load([rp 'exp/raw_sequence/data_mri_main_20150205_' ...
        num2str(sequence_name(i)) '.mat']);
    
    % Define the output size
    response{i} = zeros(number_of_trials,number_of_bins);
    
    % There are 72 sequences used in the change detection condition
    for j = 1:number_of_trials
        
        % Existence or nonexistence (1/0) of image in each bin (1*120 vector)
        image_on_off = raw_sequence.S{main_sequence_index(j)};
        change_point = 1; button_press = 1;
        
        theta_low_to_high = zeros(number_of_bins,1);
        theta_high_to_low = zeros(number_of_bins,1);
        
        for t = 1:number_of_bins
            
            if strcmp(fit_name,'main8') || strcmp(fit_name,'main10') || ...
                    strcmp(fit_name,'main13') || strcmp(fit_name,'main14')
                theta_low_to_high(t) = theta_low_to_high_fix - ...
                    (1-exp(-((t-change_point)/alpha)^beta))*theta_low_to_high_fix*(1-delta);
                theta_high_to_low(t) = theta_high_to_low_fix - ...
                    (1-exp(-((t-change_point)/alpha)^beta))*theta_high_to_low_fix*(1-delta);
            elseif strcmp(fit_name,'main19') || strcmp(fit_name,'main20') || ...
                    strcmp(fit_name,'main21') || strcmp(fit_name,'main22')
                if t >= floor(X/0.25)+change_point
                    theta_low_to_high(t) = theta_low_to_high_fix - ...
                        (1-exp(-((t-button_press)/alpha)^beta))*theta_low_to_high_fix*(1-delta);
                    theta_high_to_low(t) = theta_high_to_low_fix - ...
                        (1-exp(-((t-button_press)/alpha)^beta))*theta_high_to_low_fix*(1-delta);
                end
            else
                theta_low_to_high(t) = theta_low_to_high_fix;
                theta_high_to_low(t) = theta_high_to_low_fix;
            end
            
            switch method_name
                
                case 1
                    
                    % The computation of difference starts X seconds after
                    % a perceived change location. Secure data points for X first
                    % and then increase the number of data points in Y
                    
                    if t >= floor(X/0.25)+change_point
                        % Calculate the unit time mean image frequency difference
                        bins_to_includeX = t-floor(X/0.25)+1:t;
                        muX = 4*sum(image_on_off(bins_to_includeX))/length(bins_to_includeX);
                        bins_to_includeY = max(change_point,t-floor((X+Y)/0.25)+1):(t-floor(X/0.25));
                        muY = 4*sum(image_on_off(bins_to_includeY))/length(bins_to_includeY);
                        diffXY = muX - muY;
                        if diffXY >= max(theta_low_to_high(t)-k*min(muX,muY),0.001) || ...
                                diffXY <= min(theta_high_to_low(t)+k*min(muX,muY),-0.001)
                            response{i}(j,t) = 1;
                            change_point = t-floor(X/0.25)+1;
                            button_press = t;
                        end
                    end
                    
                case 2
                    
                    % The computation of difference starts X seconds after
                    % a perceived change location. Secure data points for Y first
                    % and then increase the number of data points in X
                    
                    if (change_point==1 && t >= floor(Y/0.25)+1) || ...
                            (change_point~=1 && (floor(X/0.25) >= floor(Y/0.25)) && (t >= floor(X/0.25)+change_point)) || ...
                            (change_point~=1 && (floor(X/0.25) < floor(Y/0.25)) &&  (t >= floor(Y/0.25)+change_point))
                        % Calculate the unit time mean image frequency difference
                        binY_start = max(change_point,t-floor((X+Y)/0.25)+1);
                        bins_to_includeY = binY_start:(binY_start+floor(Y/0.25)-1);
                        muY = 4*sum(image_on_off(bins_to_includeY))/length(bins_to_includeY);
                        bins_to_includeX = (binY_start+floor(Y/0.25)):t;
                        muX = 4*sum(image_on_off(bins_to_includeX))/length(bins_to_includeX);
                        diffXY = muX - muY;
                        if diffXY >= max(theta_low_to_high(t)-k*min(muX,muY),0.001) || ...
                                diffXY <= min(theta_high_to_low(t)+k*min(muX,muY),-0.001)
                            response{i}(j,t) = 1;
                            change_point = t-length(bins_to_includeX)+1;
                            button_press = t;
                        end
                    end
                    
                case 3
                    
                    % The computation of difference starts X+Y seconds after
                    % a perceived change location.
                    
                    if t >= floor((X+Y)/0.25)+change_point-1
                        % Calculate the unit time mean image frequency difference
                        muX = sum(image_on_off(t-floor(X/0.25)+1:t),2)/(floor(X/0.25)*0.25);
                        muY = sum(image_on_off((t-floor((X+Y)/0.25)+1):(t-floor(X/0.25))))/((floor((X+Y)/0.25)-floor(X/0.25))*0.25);
                        diffXY = muX - muY;
                        if diffXY >= max(theta_low_to_high(t)-k*min(muX,muY),0.001) || ...
                                diffXY <= min(theta_high_to_low(t)+k*min(muX,muY),-0.001)
                            response{i}(j,t) = 1;
                            change_point = t-floor(X/0.25)+1;
                        end
                    end
                    
                case 4
                    
                    % The computation of difference starts after T seconds
                    % from a button press. Secure data points for X first
                    % and then increase the number of data points in Y
                    
                    if (change_point == 1 && t >= max(floor(X/0.25),floor(T/0.25))+1) || ...
                            (change_point ~= 1 && t >= floor(T/0.25)+button_press)
                            % Calculate the unit time mean image frequency difference
                            bins_to_includeX = t-floor(X/0.25)+1:t;
                            muX = 4*sum(image_on_off(bins_to_includeX))/length(bins_to_includeX);
                            bins_to_includeY = max(change_point,t-floor((X+Y)/0.25)+1):(t-floor(X/0.25));
                            muY = 4*sum(image_on_off(bins_to_includeY))/length(bins_to_includeY);
                            diffXY = muX - muY;
                            if diffXY >= max(theta_low_to_high(t)-k*min(muX,muY),0.001) || ...
                                    diffXY <= min(theta_high_to_low(t)+k*min(muX,muY),-0.001)
                                response{i}(j,t) = 1;
                                change_point = t-floor(X/0.25)+1;
                                button_press = t;
                            end
                    end
                    
                case 5
                    
                    % The computation of difference starts after T seconds
                    % from a button press. Secure data points for Y first
                    % and then increase the number of data points in X
                    
                    if (change_point ~= 1 && (floor(X/0.25) + floor(T/0.25) >= floor(Y/0.25)) && (t >= floor(T/0.25)+button_press)) || ...
                            (change_point == 1 && (floor(X/0.25) + floor(T/0.25) >= floor(Y/0.25)) && (t >= max(floor(Y/0.25),floor(T/0.25))+1)) || ...
                            ((floor(X/0.25) + floor(T/0.25) < floor(Y/0.25)) &&  (t >= floor(Y/0.25)+change_point))
                        
                        % Calculate the unit time mean image frequency difference
                        binY_start = max(change_point,t-floor((X+Y)/0.25)+1);
                        bins_to_includeY = binY_start:(binY_start+floor(Y/0.25)-1);
                        muY = 4*sum(image_on_off(bins_to_includeY))/length(bins_to_includeY);
                        bins_to_includeX = (binY_start+floor(Y/0.25)):t;
                        muX = 4*sum(image_on_off(bins_to_includeX))/length(bins_to_includeX);
                        diffXY = muX - muY;
                        if diffXY >= max(theta_low_to_high(t)-k*min(muX,muY),0.001) || ...
                                diffXY <= min(theta_high_to_low(t)+k*min(muX,muY),-0.001)
                            response{i}(j,t) = 1;
                            change_point = t-length(bins_to_includeX)+1;
                            button_press = t;
                        end
                    end
                    
            end
            
        end
        
    end
    
end

end
