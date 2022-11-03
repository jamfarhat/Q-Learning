%%
close all;
clear;
clc;

%% Parameters
rand('state',0);
randn('state',0);

%Devices
Devices=1500;

%P_dBm=-20:2.5:20;
P_dBm=-20:5:30;
P_range = 10.^((P_dBm-30)./10);

FdB = 6;
F = 10^(FdB/10);
N0dB = -204;
N0 = 10^(N0dB/10);
B = 100e3;
N = N0*B*F;

% Relays
%Pr_dBm=15;
%Pr = 10^((Pr_dBm-30)/10);

Relays=3;
Channels_Relays=[1 2 3];

r=3;

%Stochastic Geometry Parameters
%cell_radius=100e3;
cell_radius=5e3;

%Satellites
altitude=780e3;

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

for(d=1:length(P_range))
    disp(['Power:' num2str(P_dBm(d)) '/' num2str(max(P_dBm))]);
    P=P_range(d);
    
    alpha_k_j=10.^(-(128.1+36.7.*log10(Distance))/10);
    SNR=(P/N.*alpha_k_j.*(h_Nakagami).^2);

    %epsilon_j=10.^(-(103.4+24.2.*log10(altitude))/10);
    %SNR_Satellite=(Pr/N.*epsilon_j.*(g_j_i).^2);

    [NormThroughput_SA_MultipleChannels_1(d), ndist_SA_MultipleChannels_1(d), ntotal_SA_MultipleChannels_1(d)] = SlottedAloha_MultipleChannels(Devices, Relays, Channels_Relays(1), runs, frames, slots, SNR, N, r);
    [NormThroughput_SA_MultipleChannels_2(d), ndist_SA_MultipleChannels_2(d), ntotal_SA_MultipleChannels_2(d)] = SlottedAloha_MultipleChannels(Devices, Relays, Channels_Relays(2), runs, frames, slots, SNR, N, r);
    [NormThroughput_SA_MultipleChannels_3(d), ndist_SA_MultipleChannels_3(d), ntotal_SA_MultipleChannels_3(d)] = SlottedAloha_MultipleChannels(Devices, Relays, Channels_Relays(3), runs, frames, slots, SNR, N, r);

    QTable=InitializeQTable(Devices, Channels_Relays(1), slots, runs, true);
    [NormThroughput_Qlearning_MultipleChannels_1(d), ndist_Qlearning_MultipleChannels_1(d), ntotal_Qlearning_MultipleChannels_1(d)] = Qlearning_MultipleChannels(Devices, Relays, Channels_Relays(1), runs, frames, slots, SNR, N, r, QTable, alpha, gamma);
    QTable=InitializeQTable(Devices, Channels_Relays(2), slots, runs, true);
    [NormThroughput_Qlearning_MultipleChannels_2(d), ndist_Qlearning_MultipleChannels_2(d), ntotal_Qlearning_MultipleChannels_2(d)] = Qlearning_MultipleChannels(Devices, Relays, Channels_Relays(2), runs, frames, slots, SNR, N, r, QTable, alpha, gamma);
    QTable=InitializeQTable(Devices, Channels_Relays(3), slots, runs, true);
    [NormThroughput_Qlearning_MultipleChannels_3(d), ndist_Qlearning_MultipleChannels_3(d), ntotal_Qlearning_MultipleChannels_3(d)] = Qlearning_MultipleChannels(Devices, Relays, Channels_Relays(3), runs, frames, slots, SNR, N, r, QTable, alpha, gamma);

    QTable=InitializeQTable(Devices, Channels_Relays(1), slots, runs, true);
    [NormThroughput_Qlearning_UniqueChannel_1(d), ndist_Qlearning_UniqueChannel_1(d), ntotal_Qlearning_UniqueChannel_1(d)] = Qlearning_UniqueChannel(Devices, Relays, Channels_Relays(1), runs, frames, slots, SNR, N, r, QTable, alpha, gamma);
    QTable=InitializeQTable(Devices, Channels_Relays(2), slots, runs, true);
    [NormThroughput_Qlearning_UniqueChannel_2(d), ndist_Qlearning_UniqueChannel_2(d), ntotal_Qlearning_UniqueChannel_2(d)] = Qlearning_UniqueChannel(Devices, Relays, Channels_Relays(2), runs, frames, slots, SNR, N, r, QTable, alpha, gamma);
    QTable=InitializeQTable(Devices, Channels_Relays(3), slots, runs, true);
    [NormThroughput_Qlearning_UniqueChannel_3(d), ndist_Qlearning_UniqueChannel_3(d), ntotal_Qlearning_UniqueChannel_3(d)] = Qlearning_UniqueChannel(Devices, Relays, Channels_Relays(3), runs, frames, slots, SNR, N, r, QTable, alpha, gamma);

end

redundantRate_SA_MultipleChannels_1 = 100*(1 - ndist_SA_MultipleChannels_1./ntotal_SA_MultipleChannels_1);
redundantRate_SA_MultipleChannels_2 = 100*(1 - ndist_SA_MultipleChannels_2./ntotal_SA_MultipleChannels_2);
redundantRate_SA_MultipleChannels_3 = 100*(1 - ndist_SA_MultipleChannels_3./ntotal_SA_MultipleChannels_3);
redundantRate_Qlearning_MultipleChannels_1 = 100*(1 - ndist_Qlearning_MultipleChannels_1./ntotal_Qlearning_MultipleChannels_1);
redundantRate_Qlearning_MultipleChannels_2 = 100*(1 - ndist_Qlearning_MultipleChannels_2./ntotal_Qlearning_MultipleChannels_2);
redundantRate_Qlearning_MultipleChannels_3 = 100*(1 - ndist_Qlearning_MultipleChannels_3./ntotal_Qlearning_MultipleChannels_3);
redundantRate_Qlearning_UniqueChannel_1 = 100*(1 - ndist_Qlearning_UniqueChannel_1./ntotal_Qlearning_UniqueChannel_1);
redundantRate_Qlearning_UniqueChannel_2 = 100*(1 - ndist_Qlearning_UniqueChannel_2./ntotal_Qlearning_UniqueChannel_2);
redundantRate_Qlearning_UniqueChannel_3 = 100*(1 - ndist_Qlearning_UniqueChannel_3./ntotal_Qlearning_UniqueChannel_3);

%% Plot
figure(1)
plot(P_dBm, redundantRate_SA_MultipleChannels_1, 's-', 'LineWidth', 1.5, 'Color', [0, 0.4470, 0.7410]);
hold on;
plot(P_dBm, redundantRate_Qlearning_MultipleChannels_1, 's--', 'LineWidth', 1.5, 'Color', [0, 0.4470, 0.7410]);
plot(P_dBm, redundantRate_SA_MultipleChannels_2, 'd-', 'LineWidth', 1.5, 'Color', [0.8500, 0.3250, 0.0980]);
plot(P_dBm, redundantRate_Qlearning_MultipleChannels_2, 'd--', 'LineWidth', 1.5, 'Color', [0.8500, 0.3250, 0.0980]);
plot(P_dBm, redundantRate_SA_MultipleChannels_3, 'o-', 'LineWidth', 1.5, 'Color', 	[0.9290, 0.6940, 0.1250]);
plot(P_dBm, redundantRate_Qlearning_MultipleChannels_3, 'o--', 'LineWidth', 1.5, 'Color', 	[0.9290, 0.6940, 0.1250]);
grid on
legend('SA-NOMA - C=1', 'QL-NOMA - C=1', 'SA-NOMA - C=2', 'QL-NOMA - C=2', 'SA-NOMA - C=3', 'QL-NOMA - C=3')
xlabel('Power (dBm)','Fontsize',14);
ylabel('Rate of redundant messages','Fontsize',14);


%P_plot = 1:2:length(P_dBm);
P_plot = 1:1:length(P_dBm);

figure(2)
plot(P_dBm, (NormThroughput_Qlearning_MultipleChannels_1), 's-', 'LineWidth', 1.5, 'MarkerIndices', P_plot, 'MarkerFaceColor', 'w', 'MarkerSize', 10);
hold on;
%plot(P_dBm, (NormThroughput_Qlearning_UniqueChannel_1), 's--', 'LineWidth', 1.5, 'MarkerIndices', P_plot, 'MarkerFaceColor', 'w', 'MarkerSize', 10);
plot(P_dBm, (NormThroughput_SA_MultipleChannels_1), 's:', 'LineWidth', 1.5, 'MarkerIndices', P_plot, 'MarkerFaceColor', 'w', 'MarkerSize', 10);
set(gca,'ColorOrderIndex',1)
plot(P_dBm, (NormThroughput_Qlearning_MultipleChannels_3), 'o-', 'LineWidth', 1.5, 'MarkerIndices', P_plot, 'MarkerFaceColor', 'w', 'MarkerSize', 10);
plot(P_dBm, (NormThroughput_SA_MultipleChannels_3), 'o:', 'LineWidth', 1.5, 'MarkerIndices', P_plot, 'MarkerFaceColor', 'w', 'MarkerSize', 10);
plot(P_dBm, (NormThroughput_Qlearning_UniqueChannel_3), '*-', 'LineWidth', 1.5, 'MarkerIndices', P_plot, 'MarkerFaceColor', 'w', 'MarkerSize', 10);
hold off
grid on
%legend('QL-NOMA - C=1', 'QL-NOMA - C=2', 'QL-NOMA - C=3', 'SA-NOMA - C=1', 'SA-NOMA - C=2', 'SA-NOMA - C=3', 'fontsize',12)
legend('QL-NOMA (C=1)', 'SA-NOMA (C=1)', 'QL-NOMA (C=3)', 'SA-NOMA (C=3)', 'DQL-based JRSAC [4] (C=3)', 'fontsize',12)
xlabel('Transmit Power (P) [dBm]','Fontsize',14);
ylabel('Normalized Throughput (\tau) [bps/Hz]','Fontsize',14);