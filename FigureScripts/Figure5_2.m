% 230907
% This script is meant to integrate the equation rho(x|n)=rho_0*(1+N(0|sigma))^eta*n
% which is an important quantity for determining how variation in n affects
% localization of populations around gaussian signals.
%
% To solve this equation, we will use MCMC sampling. Then, to compute the
% performance we'll count all the samples that lie between -a and a. One
% problem I forsee is that the bounds may be difficult to set. For that
% reason, I plan to start by bounding the samples to +- 4*sigma and 
% adjusting as necessary.

%% Defining relevant parameters
close all
clear

rho_0 = 1;
L0_list = [1, 6, 100];
Ki = 18;
sigma = 1;
eta = 1;
n_range = 10.^[0:0.1:2];


%% Perform sampling
initial = 0.1;
nsamples = 10000;
bounds = 5*sigma;
stored_samples = zeros(nsamples, length(n_range), length(L0_list));

% for jj = 1:length(L0_list)
%     for ii = 1:length(n_range)
%         L0 = L0_list(jj);
%         n = n_range(ii);
%         pdf = @(x) unifpdf(x, -bounds, bounds)*rho_0*(1 + (L0/Ki)*exp(-x^2/(2*sigma^2)))^(eta*n);
%         stored_samples(:, ii, jj) = slicesample(initial, nsamples, "pdf", pdf);
%     end
% end
% save('PerformanceData')
load("./Outputs/Figure6PerformanceMap.mat")

%% Plot distribution for example
figure()
subplot(2, 1, 1)
hold on
histogram(stored_samples(:, 1, 1), 'Normalization','pdf')
histogram(stored_samples(:, end, 1), 'Normalization', 'pdf')
xlabel("x")
ylabel("pdf")
title("\rho(x|n) when L_0 = 1")
legend({['n = ', num2str(n_range(1))], ['n = ', num2str(n_range(end))]})
legend boxoff

%% Calculate performance
thresh = sigma;
success_matrix = stored_samples > -thresh & stored_samples < thresh;
performance_1 = sum(success_matrix(:, :, 1), 1)./nsamples;
performance_2 = sum(success_matrix(:, :, 2), 1)./nsamples;
performance_3 = sum(success_matrix(:, :, 3), 1)./nsamples;

subplot(2, 1, 2)
hold on
lw = 3;
plot(n_range, performance_1, 'Linewidth', lw, 'color', [0.8, 0.8, 0.8])
plot(n_range, performance_2, 'Linewidth', lw, 'color', [0.5, 0.5, 0.5])
plot(n_range, performance_3, 'Linewidth', lw, 'color', [0.3, 0.3, 0.3])
legend({['L_0 = ', num2str(L0_list(1))], ['L_0 = ', num2str(L0_list(2))],...
    ['L_0 = ', num2str(L0_list(3))]}, 'Position', [0.717906188293507,0.127097823560857,0.26859956544129,0.101506026336946])
legend boxoff
xlabel("n")
ylabel("performance")
title("Q(n) depends on L_0")
set(gca, 'xscale', 'log')

% set(gcf, 'Position', [675,151,547,804], 'Renderer', 'painters')
set(gcf, 'Position', [488, 98, 366, 664], 'Renderer', 'painters')
saveas(gcf, './Outputs/Figure6_part2.png')
saveas(gcf, './Outputs/Figure6_part2.svg')
