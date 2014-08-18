

xtick = get(gca,'xtick');
xticklabel = get(gca,'xticklabel');

h = get(gca,'xlabel');
xlabelstring = get(h,'string');
xlabelposition = get(h,'position');

yposition = xlabelposition(2);
yposition = repmat(yposition,length(xticklabel),1);

set(gca,'xtick',[]);

hnew = text(xtick, yposition, xticklabel);
set(hnew,'rotation',90,'horizontalalignment','right');
