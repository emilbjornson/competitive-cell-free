%This Matlab script generates Figure 5 in the paper:
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


%% Define simulation setup

%Number of setups with random UE locations
nbrOfSetups = 200;

%Number of channel realizations per setup
nbrOfRealizations = 1000;

%Number of APs in the cell-free network
L = 400;

%Number of UEs in the network
K = 40;

%Number of antennas per AP
N = 1;

%Length of the coherence block
tau_c = 200;

%Number of pilots per coherence block
tau_p = 10;

%Uplink transmit power per UE (mW)
p = 100;


%Prepare to save simulation results
sumSE_AP_MR_tot = zeros(nbrOfSetups,4);
sumSE_AP_MMSE_tot = zeros(nbrOfSetups,4);
sumSE_AP_SIC_tot = zeros(nbrOfSetups,1);


%% Go through all setups
for n = 1:nbrOfSetups
    
    %Display simulation progress
    disp(['Setup ' num2str(n) ' out of ' num2str(nbrOfSetups)]);
    
    %Generate one setup with UEs at random locations
    [R_AP,~,pilotIndex] = generateSetup(L,K,N,1,1);
    
    %Generate channel realizations, channel estimates, and estimation
    %error correlation matrices for all UEs to the APs 
    [Hhat_AP,H_AP,B_AP] = functionChannelEstimates(R_AP,nbrOfRealizations,L,K,N,tau_p,pilotIndex,p);
    
    %Compute SE for the Cell-free mMIMO system with Monte Carlo simulations
    [SE_AP_MR,SE_AP_MMSE,sumSE_SIC] = functionComputeSE_AP_uplink(Hhat_AP,H_AP,R_AP,B_AP,tau_c,tau_p,nbrOfRealizations,N,K,L,p);


    %Save SE values
    sumSE_AP_MR_tot(n,:) = sum(SE_AP_MR,1);
    sumSE_AP_MMSE_tot(n,:) = sum(SE_AP_MMSE,1);
    sumSE_AP_SIC_tot(n) = sumSE_SIC;
    
    %Remove large matrices at the end of analyzing this setup
    clear B_AP D_AP H_AP Hhat_AP R_AP;
    
end


%% Plot simulation results
figure;
hold on; box on;
plot(sort(sumSE_AP_SIC_tot),linspace(0,1,nbrOfSetups),'r-.','LineWidth',4);
plot(sort(sumSE_AP_MMSE_tot(:,4)),linspace(0,1,nbrOfSetups),'b-','LineWidth',2);
plot(sort(sumSE_AP_MMSE_tot(:,3)),linspace(0,1,nbrOfSetups),'b--','LineWidth',2);
plot(sort(sumSE_AP_MMSE_tot(:,2)),linspace(0,1,nbrOfSetups),'b-.','LineWidth',2);
plot(sort(sumSE_AP_MMSE_tot(:,1)),linspace(0,1,nbrOfSetups),'b:','LineWidth',2);
xlabel('Sum spectral efficiency [bit/s/Hz]','Interpreter','Latex');
ylabel('CDF','Interpreter','Latex');
legend({'L4 (MMSE-SIC)','L4 (MMSE)','L3 (L-MMSE)','L2 (L-MMSE)','L1 (Small cells)'},'Interpreter','Latex','Location','NorthWest');
