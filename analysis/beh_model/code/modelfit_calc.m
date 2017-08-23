%%
% Returns the resulting fmeasure from modelfit.

% The inputs are the followings.
% modelName: string (e.g. 'adamsmackay3')
% fitName: string (e.g. 'main17')
% methodName: integer (e.g. 4)
% distanceTime: integer (e.g. 5)
% parameterValue:
% subjectName: string (e.g. '002')
% cond: use all trials ('') or hold-out validation ('_pred')
%%

function out_val = modelfit_calc(modelName, fitName, methodName, distanceTime, ...
    parameterValue, subjectName, cond)

% Get the model response for each sequence and for each methods.
% Output 'response' is 2*1 cell: column is sequence number
eval(['response = resp_', modelName, '(parameterValue, fitName, methodName, subjectName);']);
% fIndividualCrossValidationH is performance for [odd, even] trials
[fIndividual,~,~,~,~,~,~,~,~,~,~,~,~,~,~,~,~,~,fIndividualCrossValidationH] = ...
    model_subject_response_evaluator(response, distanceTime, subjectName);

if strcmp(cond,'')
    out_val = -fIndividual;
elseif strcmp(cond,'_pred')
    out_val = -fIndividualCrossValidationH(1);
end

end