%% 
% THIS SCRIPT CREATES FIGURES IN MANUSCRIPT RELATED TO BEHAVIORAL SUMMARY

% FIG1; Figure 2A
% FIG2; Figure 2B
% FIG3; Figure S2B
% FIG4; Figure S2C
% FIG5; Figure S2C
%%

function beh_figure(odd)

close all

cd(fileparts(mfilename('fullpath')));

figOn = [1,1,1,1,1];

LGFS = 20; AXFS = 28; LBFS = 28; mSize = 10;

% Load the subject IDs used for the analysis
load('../../../subjectList.mat');

for iSub = 1:length(subjects)
    
    basic = load(['../result/cdm' subjects{iSub} '.mat']);
    fDifficulty(iSub,:) = [basic.fMainDifficulty, basic.fOdd, basic.fMain];
    precisionDifficulty(iSub,:) = [basic.precisionMainDifficulty, basic.precisionOdd];
    sensitiveDifficulty(iSub,:) = [basic.sensitiveMainDifficulty, basic.sensitiveOdd];
    bPressDifficulty(iSub,:) = [basic.mainBpressDifficulty, basic.oddBpress];
    fMainOdd(iSub,:) = [basic.fMain, basic.fOdd];
    precisionMainOdd(iSub,:) = [basic.precisionMain, basic.precisionOdd];
    sensitiveMainOdd(iSub,:) = [basic.sensitiveMain, basic.sensitiveOdd];
    bPressMainOdd(iSub,:) = [basic.mainBpress, basic.oddBpress];
    
    clear basic;
    
end

iFig=0; x=[1,1.7];
if figOn(1)
    iFig = iFig + 1;
    figure(iFig);
    for i=1:barN
        h=bar(x(i),mean(bPressDifficulty(:,i)),0.4,'LineWidth',2,'FaceColor','none','EdgeColor',col{i});
        hold on
        xdata = get(h,'XData');
        plot(xdata,bPressDifficulty(:,i),'o','Color',col{i},'MarkerSize',mSize);
        h1 = errorbar(xdata,mean(bPressDifficulty(:,i)),...
            std(bPressDifficulty(:,i),0,1)./sqrt(length(subjects)),'Color',col{i},'LineWidth',2);
        removeErrorBarEnds(h1);
    end
    [h,p,~,stats] = ttest(bPressDifficulty(:,1),bPressDifficulty(:,2))
    [h,p,~,stats] = ttest(bPressDifficulty(:,1)-1.0278)
    [h,p,~,stats] = ttest(bPressDifficulty(:,2)-1.0278)
    
    set(h1,'linestyle','none');
    plot([0,2.5],[1.0278,1.0278],'k--','LineWidth',2)
    set(gca,'Box','off','TickDir','Out');
    set(gca,'XLim',[0.6 2.1]);
    set(gca,'XTick',[]);
    set(gca,'YLim',[0 2.55]);
    set(gca,'YTick',[0:0.5:2.5]);
    ylabel('Button press per trial','FontSize',LBFS);
    set(gca,'FontSize',AXFS);
end

if figOn(2)
    rng('default');
    for iSub = 1:length(subjects)
        data = load(['../result/cdm' subjects{iSub} '.mat']);
        % Button press per bin for each subject
        nBP(iSub) = data.mainBpress/120;
    end
    
    nSim = 5000; distanceTime = 5;
    
    % Define the output size
    fIndividual = zeros(length(subjects),nSim);
    
    for iSim = 1:nSim       
        for iSub = 1:length(subjects)            
            % For individual
            % Get the response
            response = {rand(72,120) <= nBP(iSub)};
            % Evaluate the model output according to f measure
            f(iSub,iSim) = subject_actual_change_evaluator(response, distanceTime, subjects{iSub});
        end        
    end
    
    fIndividual = mean(f,2);
    meanfIndividual = mean(fIndividual);
    semIndividual = std(fIndividual,0,1)/sqrt(length(subjects));
    
    iFig = iFig + 1;
    figure(iFig);
    for i = 1:barN
        h=bar(x(i),mean(fDifficulty(:,i)),0.6,'LineWidth',2,'FaceColor','none','EdgeColor',col{i});
        hold on
        xdata = get(h,'XData');
        plot(xdata,fDifficulty(:,i),'o','Color',col{i},'MarkerSize',mSize);
        h1 = errorbar(xdata,mean(fDifficulty(:,i)),...
            std(fDifficulty(:,i),0,1)./sqrt(length(subjects)),'Color',col{i},'LineWidth',2);
        removeErrorBarEnds(h1);
        set(h1,'linestyle','none');
    end
    h=bar(3.1,meanfIndividual,0.6,'LineWidth',2,'FaceColor','none','EdgeColor',grey);
    xdata = get(h,'XData');
    plot(xdata,fIndividual,'o','Color',grey,'MarkerSize',mSize);  
    h1 = errorbar(xdata,meanfIndividual,semIndividual,'Color',grey,'LineWidth',2);
    removeErrorBarEnds(h1);
    set(gca,'Box','off','TickDir','Out');
    set(gca,'XTick',[]);
    set(gca,'XLim',[0.5 3.6]);
    set(gca,'YLim',[0 1.05]);
    set(gca,'YTick',[0:0.2:1.0]);
    ylabel('Performance f measure','FontSize',LBFS);
    set(gca,'FontSize',AXFS);
    [h p,~,stats] = ttest(fDifficulty(:,1),fIndividual)
    [h p,~,stats] = ttest(fDifficulty(:,2),fIndividual) 
end

x=[1,1.6]; y=[2.5,3.1];
if figOn(3)
    iFig = iFig + 1;
    figure(iFig);
    for i=1:barN
        h=bar(x(i),mean(precisionDifficulty(:,i)),0.4,'LineWidth',2,'FaceColor','none','EdgeColor',col{i});
        hold on
        xdata = get(h,'XData');
        plot(xdata,precisionDifficulty(:,i),'o','Color',col{i},'MarkerSize',mSize);
        h1 = errorbar(xdata,mean(precisionDifficulty(:,i)),...
            std(precisionDifficulty(:,i),0,1)./sqrt(length(subjects)),'Color',col{i},'LineWidth',2);
        removeErrorBarEnds(h1);
    end
    [h,p] = ttest(precisionDifficulty(:,1),precisionDifficulty(:,2))
    for i=1:barN
        h=bar(y(i),mean(sensitiveDifficulty(:,i)),0.4,'LineWidth',2,'FaceColor','none','EdgeColor',col{i});
        hold on
        xdata = get(h,'XData');
        plot(xdata,sensitiveDifficulty(:,i),'o','Color',col{i},'MarkerSize',mSize);
        h1 = errorbar(xdata,mean(sensitiveDifficulty(:,i)),...
            std(sensitiveDifficulty(:,i),0,1)./sqrt(length(subjects)),'Color',col{i},'LineWidth',2);
        removeErrorBarEnds(h1);
    end
    [h,p] = ttest(sensitiveDifficulty(:,1),sensitiveDifficulty(:,2))
    set(h1,'linestyle','none');
    set(gca,'Box','off','TickDir','Out');
    set(gca,'XLim',[0.5 3.6]);
    set(gca,'XTick',[]);
    set(gca,'YLim',[0 1]);
    set(gca,'YTick',[0:0.2:1.05]);
    ylabel('Performance','FontSize',LBFS);
    set(gca,'FontSize',AXFS);
end

x1=[1,1.6]; x2=[2.5,3.1]; x3=[4,4.6];
col = {[0,0,0],[1,0.4,0]};
if figOn(4)
    iFig = iFig + 1;
    figure(iFig);
    for i=1:barN
        h=bar(x1(i),mean(fMainOdd(:,i)),0.4,'LineWidth',2,'FaceColor','none','EdgeColor',col{i});
        hold on
        xdata = get(h,'XData');
        plot(xdata,fMainOdd(:,i),'o','Color',col{i},'MarkerSize',mSize);
        h1 = errorbar(xdata,mean(fMainOdd(:,i)),...
            std(fMainOdd(:,i),0,1)./sqrt(length(subjects)),'Color',col{i},'LineWidth',2);
        removeErrorBarEnds(h1);
    end
    [h,p,~,stats] = ttest(fMainOdd(:,1),fMainOdd(:,2))
    for i=1:barN
        h=bar(x2(i),mean(precisionMainOdd(:,i)),0.4,'LineWidth',2,'FaceColor','none','EdgeColor',col{i});
        hold on
        xdata = get(h,'XData');
        plot(xdata,precisionMainOdd(:,i),'o','Color',col{i},'MarkerSize',mSize);
        h1 = errorbar(xdata,mean(precisionMainOdd(:,i)),...
            std(precisionMainOdd(:,i),0,1)./sqrt(length(subjects)),'Color',col{i},'LineWidth',2);
        removeErrorBarEnds(h1);
    end
    [h,p,~,stats] = ttest(precisionMainOdd(:,1),precisionMainOdd(:,2))
    for i=1:barN
        h=bar(x3(i),mean(sensitiveMainOdd(:,i)),0.4,'LineWidth',2,'FaceColor','none','EdgeColor',col{i});
        hold on
        xdata = get(h,'XData');
        plot(xdata,sensitiveMainOdd(:,i),'o','Color',col{i},'MarkerSize',mSize);
        h1 = errorbar(xdata,mean(sensitiveMainOdd(:,i)),...
            std(sensitiveMainOdd(:,i),0,1)./sqrt(length(subjects)),'Color',col{i},'LineWidth',2);
        removeErrorBarEnds(h1);
    end
    pbaspect([3 2 1]);
    [h,p,~,stats] = ttest(sensitiveMainOdd(:,1),sensitiveMainOdd(:,2))
    set(h1,'linestyle','none');
    set(gca,'Box','off','TickDir','Out');
    set(gca,'XLim',[0.5 5.1]);
    set(gca,'XTick',[]);
    set(gca,'YLim',[0 1.1]);
    set(gca,'YTick',[0:0.5:1]);
    ylabel('Performance','FontSize',LBFS);
    set(gca,'FontSize',AXFS);
end

x1=[1,1.6]; 
col = {[0,0,0],[1,0.4,0]};
if figOn(5)
    iFig = iFig + 1;
    figure(iFig);
    for i=1:barN
        h=bar(x1(i),mean(bPressMainOdd(:,i)),0.4,'LineWidth',2,'FaceColor','none','EdgeColor',col{i});
        hold on
        xdata = get(h,'XData');
        plot(xdata,bPressMainOdd(:,i),'o','Color',col{i},'MarkerSize',mSize);
        h1 = errorbar(xdata,mean(bPressMainOdd(:,i)),...
            std(bPressMainOdd(:,i),0,1)./sqrt(length(subjects)),'Color',col{i},'LineWidth',2);
        removeErrorBarEnds(h1);
    end
    pbaspect([1 1.5 1]);
    [h,p,~,stats] = ttest(bPressMainOdd(:,1),bPressMainOdd(:,2))
    set(h1,'linestyle','none');
    set(gca,'Box','off','TickDir','Out');
    set(gca,'XLim',[0.5 2.1]);
    set(gca,'XTick',[]);
    set(gca,'YLim',[0 2.1]);
    set(gca,'YTick',[0:1:2]);
    ylabel('Button press per trial','FontSize',LBFS);
    set(gca,'FontSize',AXFS);
end

end