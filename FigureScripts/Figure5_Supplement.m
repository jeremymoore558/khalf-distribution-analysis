%This script analyzes data where meAsp and Serine are alternating.
close all; 
%     clearvars -except dataparms;
addpath('./Functions')
addpath('./Data')
addpath('./Outputs')
parms = parameters(); %Struct containing system parameters


%% Datasets: Load 
Files = ["./Data/211204/FOV1.mat", "./Data/211204/FOV2.mat", "./Data/211214/FOV1.mat", "./Data/211214/FOV2.mat"];
OutFolder = 'MeAspSerAlternating/'; 
concLevels = [1, 0.02; 1, 0.02; 1, 0.02; 1, 0.02];
dataparms.concLevels = concLevels;

%% Make output destination
OutputDest = ['./Outputs/', OutFolder];
mkdir(OutputDest)

%% Initial data pre-processing
EfretData = combineFretFiles(Files, concLevels); %Combine files into one dataset
n = length(unique(concLevels)); %Number of stimulus levels
m = 15; %Number of responses to include per stimulus level

%Reformat data into a struct with n matrices with m measurements each
doseData = reorganizeData(EfretData, n, m, OutputDest);

%Calculate response-amplitude for each individual stimulus
doseData = calcResponseAmp2(doseData, n, parms);

%Calculate population-average and single-cell dose-response data
[popAvgData, SCAvgData] = calcPopAvgDoseResponse2(m, doseData, parms, dataparms);

%% Plot response to Tar Vs. Response to Tsr
tarResponses = squeeze(popAvgData.allRData_norm(:, 1, :));
tsrResponses = squeeze(popAvgData.allRData_norm(:, 2, :));
avgTarResponses = mean(tarResponses, 1);
avgTsrResponses = mean(tsrResponses, 1);
    
figure()
plot(avgTarResponses, avgTsrResponses, 'o')
xlabel("R_{meAsp}")
ylabel("R_{ser}")
ylim([-1,.1])
xlim([-1,0.1])
set(gcf, 'Position', [598,421,451,342])
JMAxes

corrcoef(avgTarResponses, avgTsrResponses)
saveas(gcf, './Outputs/FigureS4a.svg')

%% Plot K1/2 distribution going to n distribution

Files1 = ["./Outputs/MeAspDR_OneConcentrationSets/"];

legend1 = ["MeAsp", "Glu", "L-Asp", "Theoretical"];

KiValues = [18, 700, 0.9];

ii = 1;
% Calculate n distribution
load([convertStringsToChars(Files1(ii)), 'plotData.mat'])
Lplot = CDFPlotData.Lplot;

%CDF Function used to fit the K1/2 data
CDFFunc = @(p, L) logncdf(L, p(1), p(2));
PDFFunc = @(p, L) lognpdf(L, p(1), p(2));

%Get fitted K1/2 parameters
p_opt = CDFPlotData.p_opt; %[mu, sigma]

%Anonymous function to convert K1/2 to n
Ki = KiValues(ii); %uM, Shimizu 2010
pn = @(n) PDFFunc(p_opt, log(2).*Ki./n).*(log(2).*Ki)./(n.^2);

%Fit a log-normal distribution to pn
p0N = [1, 1];
n_range = 0.1:0.01:100;
nSamples = pn(n_range);
fitFuncN = @(p, n) lognpdf(n, p(1), p(2));
inv_log_postN = @(p) (sum((nSamples - fitFuncN(p, n_range)).^2));
p_optN = fminunc(inv_log_postN, p0N);

%% Plot n distribution and K_1/2 distribution
figure();
%K1/2 distribution
subplot(1, 2, 1)
hold on
yK = PDFFunc(p_opt, Lplot);
plot(CDFPlotData.Lplot, yK, 'k', 'Linewidth', 3)
xlabel("K_{1/2}")
ylabel("pdf")
xlim([10^-3,10^2])
set(gca, 'xscale', 'log')
camroll(-90)


%%N distribution from my data
subplot(1, 2, 2)
hold on
n_range = 0.1:0.01:100;
yN = pn(n_range);
plot(n_range, yN, 'color', [.2, .2, .2], 'Linewidth', 3)
xlabel('n')
ylabel('pdf')
xticks([0, 20, 40, 60, 80, 100])
xticklabels({'0', '20', '40', '60', '80', '100'})

set(gcf, 'Position', [533,516,498,388])
camroll(-90)
saveas(gcf, './Outputs/FigureS4b.svg')