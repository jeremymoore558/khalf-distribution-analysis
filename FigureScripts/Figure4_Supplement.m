%% Figure 4 supplement
% In this figure, we try to plot the mean meAsp K1/2 as a function of
% background L-asp. Ideally, and according to the model, this relationship
% should be linear. However, non-linearities in this relationship could
% explain some of the bizzarre features of figure 4, such as the inability
% to use the same ki values when Lasp or meAsp are in the background

close all; clear;
addpath('./Functions')

figure()

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% ====Extract K1/2 values===
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
CVFiles = ["./Outputs/MeAspDR_OneConcentrationSets/", ...
    "./Outputs/0_01LAsp_meAspDR/",...
     "./Outputs/1LAsp_meAspDR/",...
     "./Outputs/10LAspBackMeAspDose/",...
    ];

% metadata
backgroundLasp = [0, 0.01, 1, 10];

[CV, means, stds, CVerr, CVlowerbound, CVupperbound, meanerr,...
    meanlowbound, meanupperbound] = extract_cv_from_files(CVFiles);

plot(backgroundLasp, means, 'o-')
hold on
errorbar(backgroundLasp, means, means-meanlowbound, meanupperbound-means)
xlabel("[L-asp]")
ylabel("\langle meAsp K_{1/2}\rangle")


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Plot Khalf CDF w/ w/o Ser
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Load data and declare plot parameters
% Files1 = {["./Outputs/0BackLAspDose/", "./Outputs/10serBackLAspDose/"]};
Files1 = {["./Outputs/0BackSerDose/", "./Outputs/10serBackLAspDose/"]};

% Color Palette To Color Each plot. Rows organized the same way
colrs = {["#4E4E4E", "#008080", "#008080", "#e02459",  "#e02459"]};

lnstyles = {["-", "--", "-", "--", "-"]};

CDFFunc = @(p, L) logncdf(L, p(1), p(2));
PDFFunc = @(p, L) normpdf(L, p(1), p(2));

plotInds = [1, 3];

xlims = [[10^-3, 10^3];...
       [10^-3, 10^3];...
       [10^-3, 10^3];...
       [10^-3, 10^3]];

xlabels = ["[L] \muM", "", "[meAsp] \muM", ""];

figure()
for ii = 1:length(Files1)
    FilesTemp = Files1{ii};
    colr = colrs{ii};
    xlimits = xlims(ii, :);
    lnstyl = lnstyles{ii};

    for jj = 1:length(FilesTemp)
        %Load Data
        load([convertStringsToChars(FilesTemp(jj)), 'plotData.mat'])
        Lplot = CDFPlotData.Lplot;

        %Generate sample curves for error bars
        rnd = CDFPlotData.MCMCSamples; %Dimensions: [mu, sigma, A (not used)]
        

        %Plot PDF on top
        subplot(4, 2, plotInds(ii, 1))
        hold on
        LplotPDF = log(min(CDFPlotData.Lplot)):0.01:log(max(CDFPlotData.Lplot));
        y = PDFFunc(CDFPlotData.p_opt, LplotPDF);
        plot(LplotPDF, y ./ max(y), 'color', colr(jj), 'Linewidth', 1.8,...
            'LineStyle', lnstyl(jj))
%         ylabel(['normalized', newline, 'pdf'])
        xlim([log(xlimits(1)), log(xlimits(2))])
        set(gca, 'XTick', [])
        if plotInds(ii, 2) /2 == round(plotInds(ii, 2)/2)
            set(gca, 'YTick', [])
        end
        if plotInds(ii, 2) /2 ~= round(plotInds(ii, 2)/2)
            ylabel('pdf')
        end
        if(jj == length(FilesTemp))
            set(gca, 'Position', get(gca, 'Position') + [0, -0.02, 0, -0.1])
        end
        

        %Plot CDF on bottom
        curveSamples = zeros(size(rnd,1), length(Lplot));
        for j = 1:size(rnd, 1)
            curveSamples(j, :) = CDFFunc(rnd(j, :), Lplot);
        end
        curvePosErr = prctile(curveSamples, 97.5);
        curveNegErr = prctile(curveSamples, 2.5);

        subplot(4, 2, plotInds(ii, 2))
        hold on
        shadeLineError(CDFPlotData.Lplot, curvePosErr, curveNegErr, colr(jj))
        errorbar(CDFPlotData.concLevels, CDFPlotData.CDFPoints, CDFPlotData.CDFPointNegErr,...
            CDFPlotData.CDFPointPosErr, 'o', 'color', colr(jj));
        C(jj) = plot(CDFPlotData.Lplot, CDFFunc(CDFPlotData.p_opt, CDFPlotData.Lplot), ...
            'color', colr(jj), 'Linewidth', 1.8, 'LineStyle', lnstyl(jj));
        xlabel(xlabels(ii))
        xlim(xlimits)
        
        set(gca, 'xscale', 'log', 'XTick', [10^-2, 10^-1, 10^0, 10^1, 10^2])
        if plotInds(ii, 2) /2 ~= round(plotInds(ii, 2)/2)
            ylabel('P(K_{1/2} < [L])')
        else
            set(gca, 'YTick', [])
        end
        
    end
end

set(gcf, 'Position', [384,0,844,898])
