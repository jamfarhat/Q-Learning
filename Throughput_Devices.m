%%
close all;
clear;
clc;

%% Parameters
rand('state',0);
randn('state',0);

%Devices
Devices=linspace(100,1500,15);
%Devices=linspace(100,1500,29);
P_dBm=10;
P = 10^((P_dBm-30)/10);

FdB = 6;
F = 10^(FdB/10);
N0dB = -204;
N0 = 10^(N0dB/10);
B = 100e3;
N = N0*B*F;

% Relays
%Pr_dBm=15;
%Pr = 10^((Pr_dBm-30)/10);

Relays=4;
Channels_Relays=[1 2 3 4];

r=3;

%Satellites
altitude=780e3;

%Stochastic Geometry Parameters
cell_radius=5e3;

%Simulation Parameters
runs=100;
slots=100;
frames=50;

alpha=0.1;
gamma=0.5;

%Nakagami-m
m=2;

%% Stochastic Geometry
Distance=StochasticGeometry(max(Devices), max(Relays), cell_radius, runs);
h_Nakagami = abs((sqrt(gamrnd(m/2, 1, max(Devices), max(Relays), runs)) + j*sqrt(gamrnd(m/2, 1, max(Devices), max(Relays), runs)))/sqrt(m));

%p1 = 0.5;
%p2 = 1;
%x = p1 + p2.*randn(max(Relays), runs);
%y = p1 + p2.*randn(max(Relays), runs);
%g_j_i = sqrt(x.^2 + y.^2);

for(d=1:length(Devices))
    disp(['Number of Devices:' num2str(Devices(d)) '/' num2str(max(Devices))]);
    
    alpha_k_j=10.^(-(128.1+36.7.*log10(Distance))/10);
    SNR=(P/N.*alpha_k_j.*(h_Nakagami).^2);

    %epsilon_j=10.^(-(103.4+24.2.*log10(altitude))/10);
    %SNR_Satellite=(Pr/N.*epsilon_j.*(g_j_i).^2);

    % SA-NOMA
    [NormThroughput_SA_MultipleChannels_1(d), ndist_SA_MultipleChannels_1(d), ntotal_SA_MultipleChannels_1(d)] = SlottedAloha_MultipleChannels(Devices(d), Relays, Channels_Relays(1), runs, frames, slots, SNR, N, r);
    [NormThroughput_SA_MultipleChannels_2(d), ndist_SA_MultipleChannels_2(d), ntotal_SA_MultipleChannels_2(d)] = SlottedAloha_MultipleChannels(Devices(d), Relays, Channels_Relays(2), runs, frames, slots, SNR, N, r);
    [NormThroughput_SA_MultipleChannels_3(d), ndist_SA_MultipleChannels_3(d), ntotal_SA_MultipleChannels_3(d)] = SlottedAloha_MultipleChannels(Devices(d), Relays, Channels_Relays(3), runs, frames, slots, SNR, N, r);
    [NormThroughput_SA_MultipleChannels_4(d), ndist_SA_MultipleChannels_4(d), ntotal_SA_MultipleChannels_4(d)] = SlottedAloha_MultipleChannels(Devices(d), Relays, Channels_Relays(4), runs, frames, slots, SNR, N, r);

    % SA - without NOMA
    [NormThroughput_SA_MultipleChannels_1_NoNOMA(d), ndist_SA_MultipleChannels_1_NoNOMA(d), ntotal_SA_MultipleChannels_1_NoNOMA(d)] = SlottedAloha_MultipleChannels_NoNOMA(Devices(d), Relays, Channels_Relays(1), runs, frames, slots, SNR, N, r);
    [NormThroughput_SA_MultipleChannels_2_NoNOMA(d), ndist_SA_MultipleChannels_2_NoNOMA(d), ntotal_SA_MultipleChannels_2_NoNOMA(d)] = SlottedAloha_MultipleChannels_NoNOMA(Devices(d), Relays, Channels_Relays(2), runs, frames, slots, SNR, N, r);
    [NormThroughput_SA_MultipleChannels_3_NoNOMA(d), ndist_SA_MultipleChannels_3_NoNOMA(d), ntotal_SA_MultipleChannels_3_NoNOMA(d)] = SlottedAloha_MultipleChannels_NoNOMA(Devices(d), Relays, Channels_Relays(3), runs, frames, slots, SNR, N, r);
    [NormThroughput_SA_MultipleChannels_4_NoNOMA(d), ndist_SA_MultipleChannels_4_NoNOMA(d), ntotal_SA_MultipleChannels_4_NoNOMA(d)] = SlottedAloha_MultipleChannels_NoNOMA(Devices(d), Relays, Channels_Relays(4), runs, frames, slots, SNR, N, r);
        
    % Q-Learning
    QTable=InitializeQTable(Devices(d), Channels_Relays(1), slots, runs, true);
    [NormThroughput_Qlearning_MultipleChannels_1(d), ndist_Qlearning_MultipleChannels_1(d), ntotal_Qlearning_MultipleChannels_1(d)] = Qlearning_MultipleChannels(Devices(d), Relays, Channels_Relays(1), runs, frames, slots, SNR, N, r, QTable, alpha, gamma);
    QTable=InitializeQTable(Devices(d), Channels_Relays(2), slots, runs, true);
    [NormThroughput_Qlearning_MultipleChannels_2(d), ndist_Qlearning_MultipleChannels_2(d), ntotal_Qlearning_MultipleChannels_2(d)] = Qlearning_MultipleChannels(Devices(d), Relays, Channels_Relays(2), runs, frames, slots, SNR, N, r, QTable, alpha, gamma);
    QTable=InitializeQTable(Devices(d), Channels_Relays(3), slots, runs, true);
    [NormThroughput_Qlearning_MultipleChannels_3(d), ndist_Qlearning_MultipleChannels_3(d), ntotal_Qlearning_MultipleChannels_3(d)] = Qlearning_MultipleChannels(Devices(d), Relays, Channels_Relays(3), runs, frames, slots, SNR, N, r, QTable, alpha, gamma);
    QTable=InitializeQTable(Devices(d), Channels_Relays(4), slots, runs, true);
    [NormThroughput_Qlearning_MultipleChannels_4(d), ndist_Qlearning_MultipleChannels_4(d), ntotal_Qlearning_MultipleChannels_4(d)] = Qlearning_MultipleChannels(Devices(d), Relays, Channels_Relays(4), runs, frames, slots, SNR, N, r, QTable, alpha, gamma);
        
    % Q-Learning without NOMA
    QTable=InitializeQTable(Devices(d), Channels_Relays(1), slots, runs, true);
    [NormThroughput_Qlearning_MultipleChannels_1_NoNOMA(d), ndist_Qlearning_MultipleChannels_1_NoNOMA(d), ntotal_Qlearning_MultipleChannels_1_NoNOMA(d)] = Qlearning_MultipleChannels_NoNOMA(Devices(d), Relays, Channels_Relays(1), runs, frames, slots, SNR, N, r, QTable, alpha, gamma);
    QTable=InitializeQTable(Devices(d), Channels_Relays(2), slots, runs, true);
    [NormThroughput_Qlearning_MultipleChannels_2_NoNOMA(d), ndist_Qlearning_MultipleChannels_2_NoNOMA(d), ntotal_Qlearning_MultipleChannels_2_NoNOMA(d)] = Qlearning_MultipleChannels_NoNOMA(Devices(d), Relays, Channels_Relays(2), runs, frames, slots, SNR, N, r, QTable, alpha, gamma);
    QTable=InitializeQTable(Devices(d), Channels_Relays(3), slots, runs, true);
    [NormThroughput_Qlearning_MultipleChannels_3_NoNOMA(d), ndist_Qlearning_MultipleChannels_3_NoNOMA(d), ntotal_Qlearning_MultipleChannels_3_NoNOMA(d)] = Qlearning_MultipleChannels_NoNOMA(Devices(d), Relays, Channels_Relays(3), runs, frames, slots, SNR, N, r, QTable, alpha, gamma);
    QTable=InitializeQTable(Devices(d), Channels_Relays(4), slots, runs, true);
    [NormThroughput_Qlearning_MultipleChannels_4_NoNOMA(d), ndist_Qlearning_MultipleChannels_4_NoNOMA(d), ntotal_Qlearning_MultipleChannels_4_NoNOMA(d)] = Qlearning_MultipleChannels_NoNOMA(Devices(d), Relays, Channels_Relays(4), runs, frames, slots, SNR, N, r, QTable, alpha, gamma);

    % Q-Learning [4] - DQL-based JRSAC
    %QTable=InitializeQTable(Devices(d), Channels_Relays(1), slots, runs, true);
    %[NormThroughput_Qlearning_UniqueChannel_1(d), ndist_Qlearning_UniqueChannel_1(d), ntotal_Qlearning_UniqueChannel_1(d)] = Qlearning_UniqueChannel(Devices(d), Relays, Channels_Relays(1), runs, frames, slots, SNR, N, r, QTable, alpha, gamma);
    %QTable=InitializeQTable(Devices(d), Channels_Relays(2), slots, runs, true);
    %[NormThroughput_Qlearning_UniqueChannel_2(d), ndist_Qlearning_UniqueChannel_2(d), ntotal_Qlearning_UniqueChannel_2(d)] = Qlearning_UniqueChannel(Devices(d), Relays, Channels_Relays(2), runs, frames, slots, SNR, N, r, QTable, alpha, gamma);
    %QTable=InitializeQTable(Devices(d), Channels_Relays(3), slots, runs, true);
    %[NormThroughput_Qlearning_UniqueChannel_3(d), ndist_Qlearning_UniqueChannel_3(d), ntotal_Qlearning_UniqueChannel_3(d)] = Qlearning_UniqueChannel(Devices(d), Relays, Channels_Relays(3), runs, frames, slots, SNR, N, r, QTable, alpha, gamma);
    %QTable=InitializeQTable(Devices(d), Channels_Relays(4), slots, runs, true);
    %[NormThroughput_Qlearning_UniqueChannel_4(d), ndist_Qlearning_UniqueChannel_4(d), ntotal_Qlearning_UniqueChannel_4(d)] = Qlearning_UniqueChannel(Devices(d), Relays, Channels_Relays(4), runs, frames, slots, SNR, N, r, QTable, alpha, gamma);
end

redundantRate_SA_MultipleChannels_1 = 100*(1 - ndist_SA_MultipleChannels_1./ntotal_SA_MultipleChannels_1);
redundantRate_SA_MultipleChannels_2 = 100*(1 - ndist_SA_MultipleChannels_2./ntotal_SA_MultipleChannels_2);
redundantRate_SA_MultipleChannels_3 = 100*(1 - ndist_SA_MultipleChannels_3./ntotal_SA_MultipleChannels_3);
redundantRate_SA_MultipleChannels_4 = 100*(1 - ndist_SA_MultipleChannels_4./ntotal_SA_MultipleChannels_4);
redundantRate_Qlearning_MultipleChannels_1 = 100*(1 - ndist_Qlearning_MultipleChannels_1./ntotal_Qlearning_MultipleChannels_1);
redundantRate_Qlearning_MultipleChannels_2 = 100*(1 - ndist_Qlearning_MultipleChannels_2./ntotal_Qlearning_MultipleChannels_2);
redundantRate_Qlearning_MultipleChannels_3 = 100*(1 - ndist_Qlearning_MultipleChannels_3./ntotal_Qlearning_MultipleChannels_3);
redundantRate_Qlearning_MultipleChannels_4 = 100*(1 - ndist_Qlearning_MultipleChannels_4./ntotal_Qlearning_MultipleChannels_4);

%redundantRate_Qlearning_UniqueChannel_1 = 100*(1 - ndist_Qlearning_UniqueChannel_1./ntotal_Qlearning_UniqueChannel_1);
%redundantRate_Qlearning_UniqueChannel_2 = 100*(1 - ndist_Qlearning_UniqueChannel_2./ntotal_Qlearning_UniqueChannel_2);
%redundantRate_Qlearning_UniqueChannel_3 = 100*(1 - ndist_Qlearning_UniqueChannel_3./ntotal_Qlearning_UniqueChannel_3);
%redundantRate_Qlearning_UniqueChannel_4 = 100*(1 - ndist_Qlearning_UniqueChannel_4./ntotal_Qlearning_UniqueChannel_4);

redundantRate_SA_MultipleChannels_1_NoNOMA = 100*(1 - ndist_SA_MultipleChannels_1_NoNOMA./ntotal_SA_MultipleChannels_1_NoNOMA);
redundantRate_SA_MultipleChannels_2_NoNOMA = 100*(1 - ndist_SA_MultipleChannels_2_NoNOMA./ntotal_SA_MultipleChannels_2_NoNOMA);
redundantRate_SA_MultipleChannels_3_NoNOMA = 100*(1 - ndist_SA_MultipleChannels_3_NoNOMA./ntotal_SA_MultipleChannels_3_NoNOMA);
redundantRate_SA_MultipleChannels_4_NoNOMA = 100*(1 - ndist_SA_MultipleChannels_4_NoNOMA./ntotal_SA_MultipleChannels_4_NoNOMA);
redundantRate_Qlearning_MultipleChannels_1_NoNOMA = 100*(1 - ndist_Qlearning_MultipleChannels_1_NoNOMA./ntotal_Qlearning_MultipleChannels_1_NoNOMA);
redundantRate_Qlearning_MultipleChannels_2_NoNOMA = 100*(1 - ndist_Qlearning_MultipleChannels_2_NoNOMA./ntotal_Qlearning_MultipleChannels_2_NoNOMA);
redundantRate_Qlearning_MultipleChannels_3_NoNOMA = 100*(1 - ndist_Qlearning_MultipleChannels_3_NoNOMA./ntotal_Qlearning_MultipleChannels_3_NoNOMA);
redundantRate_Qlearning_MultipleChannels_4_NoNOMA = 100*(1 - ndist_Qlearning_MultipleChannels_4_NoNOMA./ntotal_Qlearning_MultipleChannels_4_NoNOMA);


%% Plot
%devices_plot = 1:2:length(Devices);
devices_plot = 1:length(Devices);

figure(1)
plot(Devices, (NormThroughput_Qlearning_MultipleChannels_1),'o-','LineWidth',1.5,'MarkerIndices',devices_plot,'MarkerSize',10,'MarkerFaceColor','w');
hold on;
plot(Devices, (NormThroughput_SA_MultipleChannels_1),'*-','LineWidth',1.5,'MarkerIndices',devices_plot,'MarkerSize',10);
set(gca,'ColorOrderIndex',1)
plot(Devices, (NormThroughput_Qlearning_MultipleChannels_1_NoNOMA),'o:','LineWidth',1.5,'MarkerIndices',devices_plot,'MarkerSize',10,'MarkerFaceColor','w');
plot(Devices, (NormThroughput_SA_MultipleChannels_1_NoNOMA),'*:','LineWidth',1.5,'MarkerIndices',devices_plot,'MarkerSize',10);
hold off;
grid on
xlabel('Number of Devices (D)','fontsize',14);
ylabel('Normalized Throughput (\tau) [bps/Hz]', 'fontsize', 14);
legend('QL-NOMA', 'SA-NOMA', 'QL without NOMA', 'SA without NOMA', 'fontsize', 12)

figure(2)
plot(Devices, redundantRate_Qlearning_MultipleChannels_1,'o-','LineWidth',1.5,'MarkerIndices',devices_plot,'MarkerSize',10,'MarkerFaceColor','w');
hold on;
plot(Devices, redundantRate_SA_MultipleChannels_1,'*-','LineWidth',1.5,'MarkerIndices',devices_plot,'MarkerSize',10);
set(gca,'ColorOrderIndex',1)
plot(Devices, redundantRate_Qlearning_MultipleChannels_1_NoNOMA,'o:','LineWidth',1.5,'MarkerIndices',devices_plot,'MarkerSize',10,'MarkerFaceColor','w');
plot(Devices, redundantRate_SA_MultipleChannels_1_NoNOMA,'*:','LineWidth',1.5,'MarkerIndices',devices_plot,'MarkerSize',10);
hold off;
grid on
xlabel('Number of Devices (D)','fontsize',14);
ylabel('Rate of Redundant Messages (\rho) [%]', 'fontsize', 14);
legend('QL-NOMA', 'SA-NOMA', 'QL without NOMA', 'SA without NOMA', 'fontsize', 12)

% figure(3)
% yyaxis left
% plot(inf, inf, 'ko-', 'LineWidth', 1, 'MarkerSize', 8);
% hold on;
% plot(inf, inf, 'ko--', 'LineWidth', 1, 'MarkerSize', 8);
% plot(inf, inf, 'k*-', 'LineWidth', 1, 'MarkerSize', 8);
% plot(inf, inf, 'ko:', 'LineWidth', 1, 'MarkerSize', 8);
% plot(inf, inf, 'k*:', 'LineWidth', 1, 'MarkerSize', 8);
% plot(Devices, smooth(NormThroughput_Qlearning_MultipleChannels_3),'o-','LineWidth',1.5,'MarkerIndices',devices_plot,'MarkerSize',10,'MarkerFaceColor','w');
% plot(Devices, smooth(NormThroughput_Qlearning_UniqueChannel_3),'o--','LineWidth',1.5,'MarkerIndices',devices_plot,'MarkerSize',10,'MarkerFaceColor','w');
% plot(Devices, (NormThroughput_SA_MultipleChannels_3),'*-','LineWidth',1.5,'MarkerIndices',devices_plot,'MarkerSize',10);
% plot(Devices, (NormThroughput_Qlearning_MultipleChannels_3_NoNOMA),'o:','LineWidth',1.5,'MarkerIndices',devices_plot,'MarkerSize',10,'MarkerFaceColor','w');
% plot(Devices, (NormThroughput_SA_MultipleChannels_3_NoNOMA),'*:','LineWidth',1.5,'MarkerIndices',devices_plot,'MarkerSize',10);
% hold off;
% grid on
% xlabel('Number of Devices (D)','fontsize',14);
% ylabel('Normalized Throughput (\tau) [bps/Hz]', 'fontsize', 14);
% 
% yyaxis right
% plot(Devices, redundantRate_Qlearning_MultipleChannels_3,'o-','LineWidth',1.5,'MarkerIndices',devices_plot,'MarkerSize',10,'MarkerFaceColor','w');
% hold on;
% plot(Devices, redundantRate_Qlearning_UniqueChannel_3,'o--','LineWidth',1.5,'MarkerIndices',devices_plot,'MarkerSize',10,'MarkerFaceColor','w');
% plot(Devices, redundantRate_SA_MultipleChannels_3,'*-','LineWidth',1.5,'MarkerIndices',devices_plot,'MarkerSize',10);
% plot(Devices, redundantRate_Qlearning_MultipleChannels_3_NoNOMA,'o:','LineWidth',1.5,'MarkerIndices',devices_plot,'MarkerSize',10,'MarkerFaceColor','w');
% plot(Devices, redundantRate_SA_MultipleChannels_3_NoNOMA,'*:','LineWidth',1.5,'MarkerIndices',devices_plot,'MarkerSize',10);
% hold off;
% ylabel('Rate of Redundant Messages (\rho) [%]', 'fontsize', 14);
% legend('QL-NOMA', 'DQL-based JRSAC', 'SA-NOMA', 'QL without NOMA', 'SA without NOMA', 'fontsize', 12)
% 
% 
% 
% figure(4)
% yyaxis left
% plot(inf, inf, 'ko-', 'LineWidth', 1, 'MarkerSize', 8);
% hold on;
% plot(inf, inf, 'ko--', 'LineWidth', 1, 'MarkerSize', 8);
% plot(inf, inf, 'k*-', 'LineWidth', 1, 'MarkerSize', 8);
% plot(inf, inf, 'ko:', 'LineWidth', 1, 'MarkerSize', 8);
% plot(inf, inf, 'k*:', 'LineWidth', 1, 'MarkerSize', 8);
% plot(Devices, smooth(NormThroughput_Qlearning_MultipleChannels_4),'o-','LineWidth',1.5,'MarkerIndices',devices_plot,'MarkerSize',10,'MarkerFaceColor','w');
% plot(Devices, smooth(NormThroughput_Qlearning_UniqueChannel_4),'o--','LineWidth',1.5,'MarkerIndices',devices_plot,'MarkerSize',10,'MarkerFaceColor','w');
% plot(Devices, (NormThroughput_SA_MultipleChannels_4),'*-','LineWidth',1.5,'MarkerIndices',devices_plot,'MarkerSize',10);
% plot(Devices, (NormThroughput_Qlearning_MultipleChannels_4_NoNOMA),'o:','LineWidth',1.5,'MarkerIndices',devices_plot,'MarkerSize',10,'MarkerFaceColor','w');
% plot(Devices, (NormThroughput_SA_MultipleChannels_4_NoNOMA),'*:','LineWidth',1.5,'MarkerIndices',devices_plot,'MarkerSize',10);
% hold off;
% grid on
% xlabel('Number of Devices (D)','fontsize',14);
% ylabel('Normalized Throughput (\tau) [bps/Hz]', 'fontsize', 14);
% 
% yyaxis right
% plot(Devices, redundantRate_Qlearning_MultipleChannels_4,'o-','LineWidth',1.5,'MarkerIndices',devices_plot,'MarkerSize',10,'MarkerFaceColor','w');
% hold on;
% plot(Devices, redundantRate_Qlearning_UniqueChannel_4,'o--','LineWidth',1.5,'MarkerIndices',devices_plot,'MarkerSize',10,'MarkerFaceColor','w');
% plot(Devices, redundantRate_SA_MultipleChannels_4,'*-','LineWidth',1.5,'MarkerIndices',devices_plot,'MarkerSize',10);
% plot(Devices, redundantRate_Qlearning_MultipleChannels_4_NoNOMA,'o:','LineWidth',1.5,'MarkerIndices',devices_plot,'MarkerSize',10,'MarkerFaceColor','w');
% plot(Devices, redundantRate_SA_MultipleChannels_4_NoNOMA,'*:','LineWidth',1.5,'MarkerIndices',devices_plot,'MarkerSize',10);
% hold off;
% ylabel('Rate of Redundant Messages (\rho) [%]', 'fontsize', 14);
% legend('QL-NOMA', 'DQL-based JRSAC', 'SA-NOMA', 'QL without NOMA', 'SA without NOMA', 'fontsize', 12)