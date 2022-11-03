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

Relays_Range=1:15;
Channels_Relays_Range=Relays_Range;

r=3;

% Relays
%Pr_dBm=15;
%Pr = 10^((Pr_dBm-30)/10);

%Satellites
altitude=780e3;

%Stochastic Geometry Parameters
%cell_radius=100e3;
cell_radius=5e3;

%Simulation Parameters
runs=100;
slots=100;
frames=10;

alpha=0.1;
gamma=0.5;

%Nakagami-m
m=2;

%% Stochastic Geometry
Distance=StochasticGeometry(max(Devices), max(Relays_Range), cell_radius, runs);
h_Nakagami = abs((sqrt(gamrnd(m/2, 1, max(Devices), max(Relays_Range), runs)) + j*sqrt(gamrnd(m/2, 1, max(Devices), max(Relays_Range), runs)))/sqrt(m));

for(d=1:length(Relays_Range))
    disp(['Number of Relays:' num2str(Relays_Range(d)) '/' num2str(max(Relays_Range))]);
    
    Relays=Relays_Range(d);
    Channels_Relays=Channels_Relays_Range(d);

    %p1 = 0.5;
    %p2 = 1;
    %x = p1 + p2.*randn(max(Relays), runs);
    %y = p1 + p2.*randn(max(Relays), runs);
    %g_j_i = sqrt(x.^2 + y.^2);
    
    alpha_k_j=10.^(-(128.1+36.7.*log10(Distance))/10);
    SNR=(P/N.*alpha_k_j.*(h_Nakagami).^2);

    %epsilon_j=10.^(-(103.4+24.2.*log10(altitude))/10);
    %SNR_Satellite=(Pr/N.*epsilon_j.*(g_j_i).^2);

    [NormThroughput_SA_MultipleChannels(d), ndist_SA_MultipleChannels(d), ntotal_SA_MultipleChannels(d)] = SlottedAloha_MultipleChannels(Devices, Relays, Channels_Relays, runs, frames, slots, SNR, N, r);
    
    QTable=InitializeQTable(Devices, Channels_Relays, slots, runs, true);
    [NormThroughput_Qlearning_MultipleChannels(d), ndist_Qlearning_MultipleChannels(d), ntotal_Qlearning_MultipleChannels(d)] = Qlearning_MultipleChannels(Devices, Relays, Channels_Relays, runs, frames, slots, SNR, N, r, QTable, alpha, gamma);
    
    QTable=InitializeQTable(Devices, Channels_Relays, slots, runs, true);
    [NormThroughput_Qlearning_UniqueChannel(d), ndist_Qlearning_UniqueChannel(d), ntotal_Qlearning_UniqueChannel(d)] = Qlearning_UniqueChannel(Devices, Relays, Channels_Relays, runs, frames, slots, SNR, N, r, QTable, alpha, gamma);
end

redundantRate_SA_MultipleChannels = 100*(1 - ndist_SA_MultipleChannels./ntotal_SA_MultipleChannels);
redundantRate_Qlearning_MultipleChannels = 100*(1 - ndist_Qlearning_MultipleChannels./ntotal_Qlearning_MultipleChannels);
redundantRate_Qlearning_UniqueChannel = 100*(1 - ndist_Qlearning_UniqueChannel./ntotal_Qlearning_UniqueChannel);

%% Plot
figure(1)
bar_handle=bar(Relays_Range,[redundantRate_SA_MultipleChannels; redundantRate_Qlearning_MultipleChannels; redundantRate_Qlearning_UniqueChannel]', 'group');
hold on;
bar_handle(1).FaceColor = [0, 0.4470, 0.7410];
bar_handle(2).FaceColor = [0.8500, 0.3250, 0.0980];
bar_handle(3).FaceColor = [0.9290, 0.6940, 0.1250];
grid on
xlabel('Number of Relays = Number of Channels','fontsize',13);
ylabel('Rate of Redundant Messages (%)', 'fontsize', 13);
legend('SA-NOMA', 'QL-NOMA', 'QL Orthogonal Relays/Channels');
xlim([min(Relays_Range)-1.5 max(Relays_Range)+1.5]);



figure(2)
%yyaxis left
plot(Relays_Range, (NormThroughput_Qlearning_MultipleChannels), 'o-', 'LineWidth', 1.5, 'MarkerSize', 10, 'MarkerFaceColor', 'w')
hold on
plot(Relays_Range, (NormThroughput_Qlearning_UniqueChannel), 's-', 'LineWidth', 1.5, 'MarkerSize', 10, 'MarkerFaceColor', 'w')
plot(Relays_Range, (NormThroughput_SA_MultipleChannels), '*-', 'LineWidth', 1.5, 'MarkerSize', 10, 'MarkerFaceColor', 'w')
hold off;
grid on
xlabel('Number of Relays = Number of Channels (R = C)','fontsize',14);
ylabel('Normalized Throughput (\tau) [bps/Hz]', 'fontsize', 14);
legend('QL-NOMA', 'DQL-based JRSAC [4]', 'SA-NOMA','fontsize',12);

%yyaxis right
%plot(Relays_Range, smooth(redundantRate_Qlearning_MultipleChannels), 'o--', 'LineWidth', 1.5, 'MarkerSize', 10, 'MarkerFaceColor', 'w')
%hold on
%plot(Relays_Range, smooth(redundantRate_Qlearning_UniqueChannel), 's--', 'LineWidth', 1.5, 'MarkerSize', 10, 'MarkerFaceColor', 'w')
%plot(Relays_Range, smooth(redundantRate_SA_MultipleChannels), '*--', 'LineWidth', 1.5, 'MarkerSize', 10, 'MarkerFaceColor', 'w')
%hold off;
%ylabel('Rate of Redundant Messages (\tau) [%]', 'fontsize', 14);
%legend('\rho: QL-NOMA', '\rho: QL-NOMA Orthogonal Relays/Channels', '\rho: SA-NOMA', '\tau: QL-NOMA', '\tau: QL-NOMA Orthogonal Relays/Channels', '\tau: SA-NOMA','fontsize',12);