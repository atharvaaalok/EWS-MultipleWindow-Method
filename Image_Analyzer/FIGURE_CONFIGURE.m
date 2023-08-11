x = rand(20,1) ;
y = rand(20,1) ;
plot(x,y)
% axis off
a = gca;
% a.Box = "off";
a.FontSize = 0.1;
% a.LabelFontSizeMultiplier = 1;

a.XAxis.Visible = 'off';
a.YAxis.Visible = 'off';

set(gca,'LooseInset',get(gca,'TightInset'));

saveas(gcf,'junk.jpg')