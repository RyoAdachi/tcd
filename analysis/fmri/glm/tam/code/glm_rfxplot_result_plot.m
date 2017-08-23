%% Plot the result of rfxplot analysis

function glm_rfxplot_result_plot(glmName,contrastName,roiNameSize,method)

close all

%AXFS = 28;
AXFS = 20;
LBFS = 28; 
LW = 2;

for i = 1:length(roiNameSize)
    
    figure(i)
    
%     res = load(['/home/radachi/research/CD/analysis/fmri/glm/glm_rfxplot/result/',...
%         glmName '_' contrastName '_' roiNameSize{i} '.mat']);

    res = load(['/home/radachi/research/CD/analysis/fmri/glm/glm_rfxplot/result/',...
        glmName '_' contrastName '_' roiNameSize{i} method '.mat']);
    
    mean_beta = nanmean(res.beta,1);
    sem_beta = nanstd(res.beta,0,1)./sqrt(21);
    
    figure;
    h = bar(1:1:length(mean_beta),mean_beta, 0.6,'LineWidth',LW,'FaceColor','none','EdgeColor','k');
    hold on
    xdata = get(h,'XData');
    h1 = errorbar(xdata,mean_beta,sem_beta,'black','LineWidth',LW);
    removeErrorBarEnds(h1);
    set(h1,'linestyle','none');
    set(gca,'Box','off','TickDir','Out');
    pbaspect([1,1.7,1]);
    %set(gca,'XLim',[0.5,length(mean_beta)+0.5]);
    set(gca,'XLim',[0.25,length(mean_beta)+0.75]);
    set(gca,'XTick',[]);
    %set(gca,'YLim',[1, 5.5]);
    %set(gca,'YTick',1:2:5);
    set(gca,'YLim',[-0.048, 0.03]);
    set(gca,'YTick',-0.04:0.02:0.02);
    ylabel('Effect size (a.u.)','FontSize',LBFS);
    set(gca,'FontSize',AXFS); set(gca,'FontName','Arial');
    saveas(gcf, ['/home/radachi/research/CD/analysis/fmri/glm/glm_rfxplot/figure/' ...
       glmName '_' roiNameSize{i} method '.eps'],'psc2');
    
    [~,p] = ttest(res.beta(:,1),res.beta(:,2))
    % [~,p] = ttest(res.beta(:,2),res.beta(:,3))
    % p = signrank(res.beta(:,1),res.beta(:,2))
    % p = signrank(res.beta(:,2),res.beta(:,3))
    
end

end