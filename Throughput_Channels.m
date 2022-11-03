%%
close all;
clear;
clc;

%% Parameters
rand('state',0);
randn('state',0);

%Devices
Devices=1500;

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
Channels_Relays=1:1:10;

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
    
for(d=1:length(Channels_Relays))
    disp(['Channels:' num2str(Channels_Relays(d)) '/' num2str(max(Channels_Relays))]);
    
    alpha_k_j=10.^(-(128.1+36.7.*log10(Distance))/10);
    SNR=(P/N.*alpha_k_j.*(h_Nakagami).^2);

    %epsilon_j=10.^(-(103.4+24.2.*log10(altitude))/10);
    %SNR_Satellite=(Pr/N.*epsilon_j.*(g_j_i).^2);

    [NormThroughput_SA_MultipleChannels(d), ndist_SA_MultipleChannels(d), ntotal_SA_MultipleChannels(d)] = SlottedAloha_MultipleChannels(Devices, Relays, Channels_Relays(d), runs, frames, slots, SNR, N, r);
    
    QTable=InitializeQTable(Devices, Channels_Relays(d), slots, runs, true);
    [NormThroughput_Qlearning_MultipleChannels(d), ndist_Qlearning_MultipleChannels(d), ntotal_Qlearning_MultipleChannels(d)] = Qlearning_MultipleChannels(Devices, Relays, Channels_Relays(d), runs, frames, slots, SNR, N, r, QTable, alpha, gamma);
end

redundantRate_SA_MultipleChannels = 100*(1 - ndist_SA_MultipleChannels./ntotal_SA_MultipleChannels);
redundantRate_Qlearning_MultipleChannels = 100*(1 - ndist_Qlearning_MultipleChannels./ntotal_Qlearning_MultipleChannels);

%% Plot
figure(1);
yyaxis left
plot(inf, inf, 'ko-', 'LineWidth', 1, 'MarkerSize', 8);
hold on;
plot(inf, inf, 'k*-', 'LineWidth', 1, 'MarkerSize', 8);
plot(Channels_Relays, smooth(NormThroughput_Qlearning_MultipleChannels), 'o-', 'LineWidth', 1.5, 'MarkerSize', 10, 'MarkerFaceColor', 'w')
plot(Channels_Relays, smooth(NormThroughput_SA_MultipleChannels), '*-', 'LineWidth', 1.5, 'MarkerSize', 10, 'MarkerFaceColor', 'w')
hold off;
grid on
xlabel('Number of Channels (C)','fontsize',14);
ylabel('Normalized Throughput (\tau) [bps/Hz]', 'fontsize', 14);

yyaxis right
plot(Channels_Relays, (redundantRate_Qlearning_MultipleChannels), 'o:', 'LineWidth', 1.5, 'MarkerSize', 10, 'MarkerFaceColor', 'w')
hold on
plot(Channels_Relays, (redundantRate_SA_MultipleChannels), '*:', 'LineWidth', 1.5, 'MarkerSize', 10, 'MarkerFaceColor', 'w')
hold off;
ylabel('Rate of Redundant Messages (\rho) [%]', 'fontsize', 14);
%legend('\tau: QL-NOMA', '\tau: SA-NOMA', '\rho: QL-NOMA', '\rho: SA-NOMA','fontsize',12);
legend('QL-NOMA', 'SA-NOMA','fontsize',12);