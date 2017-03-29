function [hdot1,hdot2] = plot_stdCrossover(fn,x,y1,y2,col1,col2,ylab,xticl,leg,pPos,pos,l)
if nargin < 12
    l = '-o';
end
y1(y1<=0) = NaN;
y2(y2<=0) = NaN;
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
    hp1 = plot(x,y1,'-','color',col2,'LineWidth',1);
    hdot1 = plot(x(end),y1(end),'o','color','black','MarkerSize',4, ...
        'MarkerFaceColor',col1,'Linewidth',1.5);
end

if ~isempty(y2)
    hp2 = plot(x,y2,l,'color',col2,'LineWidth',0.7,'MarkerSize',4, ...
        'MarkerFaceColor',col2);
    hdot2 = plot(x(end),y2(end),'o','color','black','MarkerSize',4, ...
        'MarkerFaceColor',col2,'Linewidth',1.5);
else
    hp2 = [];
    hdot2 = [];
end
if ~isempty(leg)
    hl = legend([hp1 hp2],leg,'Location','northoutside','Orientation','horizontal');
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
ymin = min([y1(:);y2(:)])-0.5;
ymax = max([y1(:);y2(:)])+0.5;
dtmp = ymax - ymin;
dy = round(dtmp/5);
if dy == 0
    dy = round(dtmp*10/5)/10;
end
dx = ceil(length(y1)/7);
kx = [1 (1+dx):dx:(length(y1)-dx/2) length(y1)];
set(gca,'LineWidth',1,'FontSize',8,'XTick',x(kx), ...
    'XTickLabel',xticl(kx),'YTick',0:dy:15)
if x(end) <= 100
    dx = 1;
else
    dx = 10;
end
xlim([min(x)-dx max(x)+dx])
ylim([ymin ymax])
box on
