%%
% This function implement Adams MacKay algorithm and returns the binary
% response(72*120) in each bin as a cell.

% Input is the following.
% (1) parameter_space: parameter vector
% (2) fit_name: e.g. main7
% (3) method_name: method name
% (4) fit to group or individual data or optimal agent
% ('group', 'subject ID' or 'opt')

% Output is the following.
% (1) 2*1 or 1*1 cell each contains 72*120 matrix of binary response for
% either either sequence tested.
%%
function response = resp_adamsmackay3(parameter_space, fit_name, method_name, fit_or_optimal)

rp = '/home/radachi/research/CD/';

sequence_name = [14,95]; % Raw sequence name used in the experiment
main_sequence_index = [1:18,25:42,49:66,73:90]; % Sequence number used in main
number_of_trials = length(main_sequence_index);
prune_threshold = 0.0001; % Threshold for pruning nodes for faster computations
number_of_bins = 120; % Number of bins
epsilon = [0.1;0.15;0.45;0.5]; % the occurrence probability used

% If subject ID is supplied in fit_or_optimal then response is a one by one cell
% and return result only for the sequence used by the particular subject
if ~strcmp(fit_or_optimal,'opt') && ~strcmp(fit_or_optimal,'group')
    sequence_name = sequence_name(-rem(str2num(fit_or_optimal),2)+2);
    response = cell(1,1);
else
    response = cell(2,1);
end

% theta
if strcmp(fit_or_optimal, 'opt')
    theta_fix = 0.5;
else
    theta_fix = parameter_space(1);
end

% T
if ismember(fit_name, {'main32','main34','main36','main38'})
    T = parameter_space(3);
else
    T = 5;
end

% alpha,beta,delta
if ismember(fit_name, {'main36','main38'})
    alpha = parameter_space(4);
    beta = parameter_space(5);
    delta = parameter_space(6);
elseif ismember(fit_name, {'main35','main37'})
    alpha = parameter_space(3);
    beta = parameter_space(4);
    delta = parameter_space(5);
end

%% HAZARD RATE SPECIFICATION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if strcmp(fit_or_optimal, 'opt')
%     H = (1/60).*triu(ones(number_of_bins,number_of_bins),0);
%     H(:,100:120) = 0; H(1:19,:) = 0;
    H = (1/80).*triu(ones(number_of_bins,number_of_bins),0);
else
    H = exp(parameter_space(2))*ones(number_of_bins,number_of_bins);
end

%%

numerator_component_positive = repmat(epsilon,1,number_of_bins);
numerator_component_negative = 1-numerator_component_positive;

for i = 1:length(sequence_name)
    
    % Load the raw data sequence
    raw_sequence = load([rp 'exp/raw_sequence/data_mri_main_20150205_' ...
        num2str(sequence_name(i)) '.mat']);
    
    % Define the output size
    response{i} = zeros(number_of_trials,number_of_bins);
    
    for j = 1:number_of_trials
        
        % Existence or nonexistence of image in each bin (1/0)
        image_on_off = raw_sequence.S{main_sequence_index(j)};
        
        % Define the matrix to store the run length distribution
        prior_run_length = zeros(number_of_bins+1,1); 
        % Initialize the run length (e.g. P(r_0=1)=1)
        prior_run_length(1,1) = 1;
        
        aT = ones(length(epsilon),1); bT = ones(length(epsilon),1);
        change_probability = zeros(1,number_of_bins);
        no_change_run_length = 1; change_point = 1; button_press = 1;
        
        theta = zeros(number_of_bins,1);
        
        for t = 1:number_of_bins
            
            if strcmp(fit_name,'main35') || strcmp(fit_name,'main36')
                theta(t) = theta_fix - (1-exp(-((t-change_point)/alpha)^beta))*theta_fix*(1-delta);
            elseif strcmp(fit_name,'main37') || strcmp(fit_name,'main38')
                theta(t) = theta_fix - (1-exp(-((t-button_press)/alpha)^beta))*theta_fix*(1-delta);               
            else
                theta(t) = theta_fix;
            end
            
            posterior_run_length = zeros(number_of_bins+1,1);
            
            % Update the no change run length
            no_change_run_length = no_change_run_length + 1;
            
            % Evaluate the predictive distribution for the new datum under
            % each run length.
            if image_on_off(t) == 1
                numerator_component = numerator_component_positive(:,1:t);
            elseif image_on_off(t) == 0
                numerator_component = numerator_component_negative(:,1:t);
            end
            denom_component = aT.*bT;
            predictive_distribution = ones(1,length(epsilon))*(numerator_component.*denom_component)...
                ./(ones(1,length(epsilon))*denom_component);            
            
            tmp = prior_run_length(1:t,1).*predictive_distribution';
            % Evaluate the no change run length probabilities.
            % Shift the probabilities up and to the right, scaled by the
            % transition probability and the predictive distribution.
            posterior_run_length(2:t+1,1) = tmp.*(1-H(1:t,t));
           
            % Evaluate the zero run length (change) probability by
            % accumulating the mass back down
            posterior_run_length(1,1) = H(1:t,t)'*tmp;
            
            % Node pruning
            posterior_run_length((posterior_run_length(:,1) < prune_threshold),1) = 0;
            
            % Normalize the run length distribution
            posterior_run_length(:,1) = posterior_run_length(:,1)./sum(posterior_run_length(:,1));
            
            % Evaluate the change point probability (probability of a change)
            % within the last 5 seconds. (21=correctT/0.25+1)
            ind = min(floor(T/0.25)+1,no_change_run_length-1);
            change_probability(t) = sum(posterior_run_length(1:ind,1));
            
            % Check for the button press, if yes, reset the run length
            if (sum(method_name==[1,2,3]) || (sum(method_name==[4,5,6]) && ...
                    no_change_run_length > floor(T/0.25)+1)) && ...
                    (change_probability(t) >= theta(t))
                response{i}(j,t) = 1;
                button_press = t;
                switch method_name
                    case {1,4} % Set to 0
                        posterior_run_length(:,1) = 0;
                        posterior_run_length(1,1) = 1;
                        no_change_run_length = 1;
                        change_point = t;
                    case {2,5} % Posterior mean
                        posterior_mean = round([0:1:number_of_bins]*posterior_run_length(:,1))+1;
                        posterior_run_length(:,1) = 0;
                        posterior_run_length(posterior_mean,1) = 1;
                        no_change_run_length = posterior_mean;
                        change_point = t-posterior_mean+1;
                    case {3,6} % Posterior mode
                        [~, I] = max(posterior_run_length(:,1));
                        posterior_run_length(:,1) = 0; 
                        posterior_run_length(I,1) = 1;
                        no_change_run_length = I;
                        change_point = t-I+1;
                end
            end
            
            % Update the parameter sets for each possible run length.
            if image_on_off(t) == 1
                aT = [ones(length(epsilon),1), aT.*numerator_component];
                bT = [ones(length(epsilon),1), bT];
            else
                aT = [ones(length(epsilon),1), aT];
                bT = [ones(length(epsilon),1), bT.*numerator_component];
            end
            
            % Posterior becomes the prior
            prior_run_length = posterior_run_length;
            
        end
        
    end
    
end

end