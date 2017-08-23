%% TRUE PRIOR

% This function implement Yu Cohen algorithm to return the binary
% response in each bin.

% Input is the following.
% (1) parameter vestor [alpha, beta, theta, H]
% (2) method to be used to reset run length after detecting changes
% (1,2,3 for 'set to zero', 'posterior mean' and 'posterior mode')
% (3) fit to group or individual data or optimal agent
% ('group', 'subject ID' or 'opt')

% Output is the following.
% (1) 2*1 or 1*1 cell each contains 72*120 matrix of binary response for
% either either sequence tested.
%%
function response = resp_yucohen4(parameter_space, fit_name, method_name, fit_or_optimal)

rp = '/home/radachi/research/CD/';

sequence_name = [14,95]; % Raw sequence name used in the experiment
main_sequence_index = [1:18,25:42,49:66,73:90]; % Sequence number used in main
number_of_trials = length(main_sequence_index);
number_of_bins = 120; % Number of bins
epsilon = [0:0.01:1]'; % How sparse do we calculate the posterior

% If subject ID is supplied in fit_or_optimal then response is a one by one cell
% and return result only for the sequence used by the particular subject
if ~ strcmp(fit_or_optimal, 'opt') && ~strcmp(fit_or_optimal,'group')
    sequence_name = sequence_name(-rem(str2num(fit_or_optimal),2)+2);
    response = cell(1,1);
else
    response = cell(2,1);
end

% theta
if strcmp(fit_name,'main23')
    theta_low_to_high_fix = parameter_space(1); % threshold to elicit button press (L to H)
    theta_high_to_low_fix = parameter_space(2); % threshold to elicit button press (H to L)
else
    theta_low_to_high_fix = parameter_space(1); % threshold to elicit button press (L to H)
    theta_high_to_low_fix = -parameter_space(1); % threshold to elicit button press (H to L)
end

% T
if ismember(fit_name, {'main3','main5','main20','main22'})
    T = parameter_space(3);
else
    T = 5;
end

% k
if ismember(fit_name, {'main4','main21'})
    k = parameter_space(3);
elseif ismember(fit_name, {'main5','main22','main23'})
    k = parameter_space(4);
else
    k = 0;
end

% alpha, beta, delta
if ismember(fit_name, {'main19'})
    alpha = parameter_space(3);
    beta = parameter_space(4);
    delta = parameter_space(5);
elseif ismember(fit_name, {'main20','main21'})
    alpha = parameter_space(4);
    beta = parameter_space(5);
    delta = parameter_space(6);
elseif ismember(fit_name, {'main22','main23'})
    alpha = parameter_space(5);
    beta = parameter_space(6);
    delta = parameter_space(7);
end

%% HAZARD RATE SPECIFICATION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if strcmp(fit_or_optimal, 'opt')
    H = (1/60).*ones(number_of_bins,1);
    H(:,100:120) = 0; H(1:19,:) = 0;
elseif strcmp(fit_name,'main23')
    H = exp(parameter_space(3))*ones(number_of_bins,1);
else
    H = exp(parameter_space(2))*ones(number_of_bins,1);
end

%%
for i = 1:length(sequence_name)
    
    % Load the sequence
    raw_sequence = load([rp 'exp/raw_sequence/data_mri_main_20150205_' ...
        num2str(sequence_name(i)) '.mat']);
    
    % Define the output size
    response{i} = zeros(number_of_trials,number_of_bins);
    
    for j = 1:length(main_sequence_index)
        
        % Data points in each bin (1/0)
        image_on_off = raw_sequence.S{main_sequence_index(j)}';
        
        % Define the size of posterior mean
        posterior = zeros(number_of_bins+1,1);
        
        % For each point in a discretized space of [0,1], calculate the
        % posterior. Here define the size of the posterior.
        P = zeros(length(epsilon),number_of_bins+1);
        
        theta_low_to_high = zeros(number_of_bins,1);
        theta_high_to_low = zeros(number_of_bins,1);
        
        % Track when the change occured
        change_point = 0;
        
        % True prior
        P(:,1) = [zeros(10,1);1/4;zeros(4,1);1/4;zeros(29,1);1/4;zeros(4,1);1/4;zeros(50,1)];
        
        posterior(1) = epsilon'*P(:,1);
        
        for t = 1:number_of_bins
            
            if strcmp(fit_name,'main19') || strcmp(fit_name,'main20') ...
                    || strcmp(fit_name,'main21') || strcmp(fit_name,'main22') ...
                    || strcmp(fit_name,'main23')
                theta_low_to_high(t) = theta_low_to_high_fix - (1-exp(-((t-change_point)/alpha)^beta))*theta_low_to_high_fix*(1-delta);
                theta_high_to_low(t) = theta_high_to_low_fix - (1-exp(-((t-change_point)/alpha)^beta))*theta_high_to_low_fix*(1-delta);
            else
                theta_low_to_high(t) = theta_low_to_high_fix;
                theta_high_to_low(t) = theta_high_to_low_fix;
            end
            
            P(:,t+1) = (epsilon.^(image_on_off(t))).*((1-epsilon).^(1-image_on_off(t)))...
                .*((1-H(t)*ones(length(epsilon),1)).*P(:,t) ...
                + H(t)*ones(length(epsilon),1).*P(:,1));
            P(:,t+1) = P(:,t+1)./sum(P(:,t+1));
            
            % Compute the posterior mean/mode
            switch method_name
                case {1,3} % Posterior mean
                    posterior(t+1) = epsilon'*P(:,t+1);
                case {2,4} % Posterior mode
                    [~,maxI] = max(P(:,t+1));
                    posterior(t+1) = epsilon(maxI);
            end
            
            % Calculate the maximum disrepancy in the posterior mean/model.
            % If it is greater than threshold theta, then report a change.
            % If the change reported, running average is computed from this point.
            
            if sum(method_name==[1,2]) || ...
                    (sum(method_name==[3,4]) && (t-change_point) >= floor(T/0.25))
                tmp = max(change_point+1,t-floor(T/0.25)+1):t;
                posterior_tmp = posterior(tmp);              
                [~,ind] = max(abs(posterior(t+1) - posterior_tmp));
                val = posterior(t+1) - posterior_tmp(ind);
                if val >= max(theta_low_to_high(t)-k*min(posterior_tmp(ind),posterior(t+1)),0.001) || ...
                        val <= min(theta_high_to_low(t)+k*min(posterior_tmp(ind),posterior(t+1)),-0.001)
                    response{i}(j,t) = 1;
                    change_point = t;
                end
            end
            
        end
        
    end
    
end

end