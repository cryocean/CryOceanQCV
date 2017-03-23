function [hdot1] = plot_percData(fn,x,y1,col1,ylab,xticl,leg,pPos,pos)
figure(fn);
set(gcf,'PaperUnits', 'centimeters','PaperPosition', [0.5 0.5 pPos]);
set(gca,'Units','centimeters','Position',[1.2 1 pos])

hold on
plot(x,y1,'-o','color',col1,'LineWidth',0.7,'MarkerSize',4, ...
    'MarkerFaceColor',col1);
hdot1 = plot(x(end),y1(end),'o','color','black','MarkerSize',4, ...
    'MarkerFaceColor',col1,'Linewidth',1.5);

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
ymin = min(y1(:));
ymin = min(round(ymin-28),60);
ylim([ymin 105])
text(1.3924,ymin+8,leg,'FontSize',8);
box on