%This Matlab script generates Figure 6 in the paper:
%
%Emil Bjornson, Luca Sanguinetti, "Making Cell-Free Massive MIMO
%Competitive With MMSE Processing and Centralized Implementation,"
%IEEE Transactions on Wireless Communications, To appear.
%
%Download article: https://arxiv.org/abs/1903.10611
%
%This is version 1.0 (Last edited: 2019-03-19)
%
%License: This code is licensed under the GPLv2 license. If you in any way
%use this code for research that results in publications, please cite our
%paper as described above.


close all;
clear;

%Number of pilots
tau_p = 10;

%Range of length of coherence blocks
tau_c = tau_p:200;

%Number of antennas per AP
N = 4;

%Number of UEs
K = 40;


%Compute number of scalars to transmit per coherence block and per AP
level4 = tau_c*N;
level23 = (tau_c - tau_p)*K;


%% Plot simulation results
figure;
hold on; box on;
plot(tau_c,level4./tau_c,'r-','LineWidth',2);
plot(tau_c,level23./tau_c,'k-.','LineWidth',2);
xlabel('Length of coherence block ($\tau_c$)','Interpreter','Latex');
ylabel('Number of complex scalars','Interpreter','Latex');
legend({'Level 4','Level 2 or 3'},'Interpreter','Latex','Location','NorthWest');
