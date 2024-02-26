%% Control script for running sensitivity distribution analysis
%This version of the code is designed to be more friendly when analyzing
%multiple datasets at the same time. The analyses are kept independent
%until it comes time to plot dose-response and cumulative distribution
%functions

clear all;
close all;
dataparms = datasets();

for ii = 1:length(dataparms)
    analyzeSensitivityDistribution(dataparms(ii));
end



