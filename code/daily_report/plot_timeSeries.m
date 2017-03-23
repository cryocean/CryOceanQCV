function [hdot1,hdot2] = plot_timeSeries(fn,x,y1,y2,col1,col2,ylab,xticl,leg,pPos,pos,l)
if nargin < 12
    l = '-o';
end
figure(fn);
set(gcf,'PaperUnits', 'centimeters','PaperPosition', [0.5 0.5 pPos]);
set(gca,'Units','centimeters','Position',[1.2 1 pos])
hold on
if ~isempty(leg)
    hp1 = plot(x,y1,l,'color',col1,'LineWidth',0.7,'MarkerSize',4, ...
        'MarkerFaceColor',col1);
    hdot1 = plot(x(end),y1(end),'o','color','black','MarkerSize',4, ...
        'MarkerFaceColor',col1,'Linewidth',1.5);
else
    hp1 = plot(x,y1,'-','color',col1,'LineWidth',1);
    plot(x(end),y1(end),'o','color','black','MarkerSize',4, ...
        'MarkerFaceColor',col1,'Linewidth',1.5)
end

if ~isempty(y2)
    hp2 = plot(x,y2,l,'color',col2,'LineWidth',0.7,'MarkerSize',4, ...
        'MarkerFaceColor',col2);
    hdot2 = plot(x(end),y2(end),'o','color','black','MarkerSize',4, ...
        'MarkerFaceColor',col2,'Linewidth',1.5);
else
    hp2 = [];
end
if ~isempty(leg)
    hl = legend([hp1 hp2],leg,'Location','SouthWest');
%     hline = findobj(hl, 'type', 'line');
%     posh = get(hline(1),'XData');
%     set(hline(1),'XData',posh.*0.5);
%     htext = findobj(hl,'type','text');
%     pos = get(htext,'Position');
%     set(hl,'Position',0.7.*pos)
    set(hl, 'Box', 'off')
    set(hl, 'Color', 'none')
end
ylabel(ylab,'FontSize',8)
dx = ceil(length(y1)/7);
kx = [1 (1+dx):dx:(length(y1)-dx/2) length(y1)];
set(gca,'LineWidth',1,'FontSize',8,'XTick',x(kx), ...
    'XTickLabel',xticl(kx),'YTick',0:10:100)
if x(end) <= 100
    dx = 1;
else
    dx = 10;
end
xlim([min(x)-dx max(x)+dx])
ymin = min([y1(:);y2(:)]);
ymin = min(round(ymin-28),60);
ylim([ymin 105])
box on
