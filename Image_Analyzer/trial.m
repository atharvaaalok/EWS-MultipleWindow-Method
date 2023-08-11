openfig("Prediction_Map_WindowSize_TT20_SW300000.00_OR0.99_SL0.05000_RPV0.0001.fig");

% axis off
a = gca;

a.FontSize = 0.01;

a.XAxis.Visible = 'off';
a.YAxis.Visible = 'off';
set(gca,'XAxisLocation','top')
set(gca,'LooseInset',get(gca,'TightInset'));

% The limits can be set according to what is the minimum window size.
xlim([80, Inf]);
ylim([300000, Inf]);

savefig("My_fig.fig");
saveas(gca, "I_am.png")