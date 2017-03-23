function boxWhiskers(z,g,units,lab,c1,c2,p)

gu = unique(g);
ndays = length(gu);
boxplot(z,g,'labels',repmat({''},1,length(z)),'widths',0.6)
labu = [lab(1:5:end);lab(end)];
[~,ix] = unique(labu);
labu = labu(sort(ix));

% set labels
ylabel(units,'FontSize',10)
set(gca,'LineWidth',1,'FontSize',10,'XTick',unique([1:5:ndays ndays]), ...
    'XTickLabel',labu)
xlim([0 length(gu)+1])

% set box color and plot mean
set(findobj(gcf,'tag','Upper Whisker'),'LineStyle','-','LineWidth',0.8)
set(findobj(gcf,'tag','Lower Whisker'),'LineStyle','-','LineWidth',0.8)
hbox = flipud(findobj(gcf,'tag','Box'));
hmed = flipud(findobj(gcf,'tag','Median'));
hold on
for j=1:length(hbox)
    patch(get(hbox(j),'XData'),get(hbox(j),'YData'),c1,'EdgeColor','black', ...
        'LineWidth',0.8);
    plot(get(hmed(j),'XData'),get(hmed(j),'YData'),'color',c2,'LineWidth',1.5)
end

% set whiskers end to the percentiles given by p
prc = zeros(2,length(gu));
for j=1:length(gu)
    kj = g == gu(j);
    prc(:,j) = prctile(z(kj),p);
end
    
hw = flipud(findobj(gca,'tag','Upper Whisker'));
for j=1:length(hw);
    ydata = get(hw(j),'YData');
    ydata(2) = prc(2,j);
    set(hw(j),'YData',ydata);
end

ha = flipud(findobj(gca,'tag','Upper Adjacent Value'));
for j=1:length(ha);
    ydata = get(ha(j),'YData');
    ydata(:) = prc(2,j);
    set(ha(j),'YData',ydata);
end

hw = flipud(findobj(gca,'tag','Lower Whisker'));
for j=1:length(hw);
    ydata = get(hw(j),'YData');
    ydata(1) = prc(1,j);
    set(hw(j),'YData',ydata);
end

ha = flipud(findobj(gca,'tag','Lower Adjacent Value'));
for j=1:length(ha);
    ydata = get(ha(j),'YData');
    ydata(:) = prc(1,j);
    set(ha(j),'YData',ydata);
end





