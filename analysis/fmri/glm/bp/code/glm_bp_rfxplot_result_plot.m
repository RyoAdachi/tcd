%% Plot the result of rfxplot analysis

function glm_bp_rfxplot_result_plot(roiNameSize)

close all

rp = '/home/radachi/research/CD/analysis/fmri/glm';
contrastName = {'mainBP', 'controlBP'};

AXFS = 28; LBFS = 28; LW = 3;

for i = 1:length(roiNameSize)
    
    figure(i);
    
    for iCon = 1:2
        
        res = load([rp '/glm_rfxplot/result/glm31_' contrastName{iCon} '_' roiNameSize{i} '.mat']);
        
        beta(:,iCon) = res.beta;
        mean_beta = nanmean(res.beta,1);
        sem_beta = nanstd(res.beta,0,1)./sqrt(21);
        
        h = bar(iCon,mean_beta,0.8,...
            'EdgeColor','black','FaceColor','white', 'LineWidth',LW);
        hold on
        xdata = get(h,'XData');
        h1 = errorbar(xdata,mean_beta,sem_beta,'black','LineWidth',LW);
        removeErrorBarEnds(h1);
        set(h1,'linestyle','none');
        
    end
    
    set(gca,'Box','off','TickDir','Out');
    set(gca,'XLim',[0.5,2.5]);
    set(gca,'XTick',[]);
    %set(gca,'YLim',[-0.025,0.025]);
    %set(gca,'YTick',-0.02:0.02:0.02);
    ylabel('Effect size (a.u.)','FontSize',LBFS);
    set(gca,'FontSize',AXFS); set(gca,'FontName','Arial');
    %saveas(gcf, ['/home/radachi/research/CD/analysis/fmri/glm/glm_rfxplot/figure/' ...
    %    glmName '_' roiNameSize '.eps'],'psc2');
    
    [h,p] = ttest(beta(:,1),beta(:,2))
    
end

end