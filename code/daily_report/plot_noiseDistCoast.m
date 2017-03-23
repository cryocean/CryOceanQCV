function plot_noiseDistCoast(fn,x,y1,y2,col1,col2,ylab,pPos,pos,dmax,Ylim1, ...
    Ylim2,dYTick,l)
if nargin < 14
    l = '-o';
end

if ~isempty(y2)
    pos(2) = pos(2)*1.8;
    pPos(2) = pPos(2)*2;
end

figure(fn);
set(gcf,'PaperUnits', 'centimeters','PaperPosition', [0.5 0.5 pPos]);
set(gca,'Units','centimeters','Position',[1.2 1 pos])


if ~isempty(y2)
    hs1 = subplot(2,1,1);
    hp1 = plot(x,y2,l,'color',col2,'LineWidth',0.7,'MarkerSize',4, ...
        'MarkerFaceColor',col2);
    
    hs2 = subplot(2,1,2);
    hp2 = plot(x,y1,l,'color',col1,'LineWidth',0.7,'MarkerSize',4, ...
    'MarkerFaceColor',col1);

    leg = {'PLRM','LRM'};
    hl1 = legend(hp1,leg{1},'Location','NorthEast');
    hl2 = legend(hp2,leg{2},'Location','NorthEast');
    set(hl1, 'Box', 'off')
    set(hl1, 'Color', 'none')
    set(hl2, 'Box', 'off')
    set(hl2, 'Color', 'none')
    set(hs1,'LineWidth',1,'FontSize',8,'XTick',0:5:dmax, ...
        'YTick',Ylim2(1):dYTick:Ylim2(2))
    set(hs2,'LineWidth',1,'FontSize',8,'XTick',0:5:dmax, ...
        'YTick',Ylim1(1):dYTick:Ylim1(2))
    ylabel(hs1,ylab,'FontSize',8)
    ylabel(hs2,ylab,'FontSize',8)
    xlim(hs1,[0 dmax+1])
    xlim(hs2,[0 dmax+1])
    ylim(hs1,Ylim2)
    ylim(hs2,Ylim1)
else
    plot(x,y1,l,'color',col1,'LineWidth',0.7,'MarkerSize',4, ...
        'MarkerFaceColor',col1);
    ylabel(ylab,'FontSize',8)
    set(gca,'LineWidth',1,'FontSize',8,'XTick',0:5:50, ...
        'YTick',Ylim1(1):dYTick:Ylim1(2))
    xlim([0 dmax+1])
    ylim(Ylim1)
end
box on
