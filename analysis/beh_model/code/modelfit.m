%%
% Saves the resulting fmeasure from modelfit.
% modelNumber: integer (e.g. 1)
% subjectName: string (e.g. '002')
%%

function modelfit(modelNumber, subjectName)

rp = '/home/radachi/research/CD/';

% Load the parameter structure
strc = load([rp 'analysis/modelfit_param_structure.mat']);
modelName = strc.param(modelNumber).modelName;
fitName = strc.param(modelNumber).fitName;
methodName = strc.param(modelNumber).methodName;
distanceTime = strc.param(modelNumber).distanceTime;
cond = strc.param(modelNumber).cond;

%% Create directories to save the result
if ~exist([rp 'analysis/beh_model/' modelName '/result/fmin'], 'dir')
    mkdir([rp 'analysis/beh_model/' modelName '/result/'], 'fmin');
end

if ~exist([rp 'analysis/beh_model/' modelName '/result/fmin/' subjectName], 'dir')
    mkdir([rp 'analysis/beh_model/' modelName '/result/fmin/'], subjectName);
end

%% Specify the parameters and options for fmin search

options = optimset('Display','iter','TolX',1e-4, 'TolFun',1e-5);
x0 = strc.param(modelNumber).x0;
lb = strc.param(modelNumber).lb;
ub = strc.param(modelNumber).ub;

%% Run fmin search

tStart = tic;
for iRep = 1:size(x0,1)
    [x(iRep,:),fval(iRep),exitflag(iRep),output(iRep)]= ...
        fminsearchbnd(@(x)modelfit_calc(modelName, fitName, ...
        methodName, distanceTime, x, subjectName, cond),x0(iRep,:),lb,ub,options);
end
tEnd = toc(tStart)/3600; % hours

% Save the result
save([rp 'analysis/beh_model/' modelName '/result/fmin/' subjectName ...
    '/fit_fmin_' fitName '_method' num2str(methodName) '_distancetime' ...
    num2str(2*distanceTime-1) cond '.mat' ],'x','fval','exitflag','output','tEnd');

end
