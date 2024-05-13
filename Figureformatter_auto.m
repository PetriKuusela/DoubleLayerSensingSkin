minval = 1e10;
maxval = -1e10;
sminval = 1e10;
smaxval = -1e10;

for ii = 1:length(files)
    filen = files{ii};
    load(filen);
    if minval > min(reco.coupl)
        minval = min(reco.coupl);
    end
    if maxval < max(reco.coupl)
        maxval = max(reco.coupl);
    end
    if sminval > min(reco.sigma)
        sminval = min(reco.sigma);
    end
    if smaxval < max(reco.sigma)
        smaxval = max(reco.sigma);
    end
    close all;
end

figure();
colormap jet;
clim([minval maxval]);
ccc = colorbar('southoutside');
ccc.FontSize = 15;
saveas(gcf, cbname);
tempim = imread(cbname);
tempim = tempim(round(0.79*size(tempim,1)):round(1*size(tempim,1)),round(0.1*size(tempim,2)):round(0.95*size(tempim,2)),:);
imwrite(tempim, cbname);

figure();
colormap jet;
clim([sminval smaxval]);
ccc = colorbar('southoutside');
ccc.FontSize = 15;
saveas(gcf, scbname);
tempim = imread(scbname);
tempim = tempim(round(0.79*size(tempim,1)):round(1*size(tempim,1)),round(0.1*size(tempim,2)):round(0.95*size(tempim,2)),:);
imwrite(tempim, scbname);


for ii = 1:length(files)

    pl = Plotter_estsig(simesh.g, simesh.H, ginv, Hinv);

    filen = files{ii};
    load(filen);
    pl.plot(reco);
    set(0, 'CurrentFigure', pl.fig2);
    axis off;
    ccc = colorbar;
    ccc.FontSize = 16;
    colorbar off;
    clim([minval maxval]);

    saveas(gcf, [filen(1:end-4) '_mod.png']);


    set(0, 'CurrentFigure', pl.fig);
    axis off;
    ccc = colorbar;
    ccc.FontSize = 16;
    colorbar off;
    clim([sminval smaxval]);

    saveas(gcf, [filen(1:end-4) '_mod_s.png']);
    close all;

end