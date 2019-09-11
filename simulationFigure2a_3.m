%This Matlab script generates Figure 2(a) and Figure 3 in the paper:
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

%Number of BSs in the cellular network (cannot be changed)
nbrBSs = 4;

%Number of antennas at the 4 BSs
M = 100;

%Number of UEs in the network
K = 40;

%Number of antennas per AP
N = 1;

%Length of the coherence block
tau_c = 200;

%Compute number of pilots per coherence block
tau_p = K/nbrBSs;

%Uplink transmit power per UE (mW)
p = 100;


%Prepare to save simulation results
SE_BS_MR_tot = zeros(K,nbrOfSetups);
SE_BS_RZF_tot = zeros(K,nbrOfSetups);
SE_BS_MMMSE_tot = zeros(K,nbrOfSetups);
SE_AP_MR_tot = zeros(K,4,nbrOfSetups);
SE_AP_MMSE_tot = zeros(K,4,nbrOfSetups);


%% Go through all setups
for n = 1:nbrOfSetups
    
    %Display simulation progress
    disp(['Setup ' num2str(n) ' out of ' num2str(nbrOfSetups)]);
    
    %Generate one setup with UEs at random locations
    [R_AP,R_BS,pilotIndex,BSassignment] = generateSetup(L,K,N,M,1);

    
    %Generate channel realizations, channel estimates, and estimation
    %error correlation matrices for all UEs to the cell-free APs
    [Hhat_AP,H_AP,B_AP] = functionChannelEstimates(R_AP,nbrOfRealizations,L,K,N,tau_p,pilotIndex,p);
    
    %Generate channel realizations, channel estimates, and estimation
    %error correlation matrices for all UEs to the cellular BSs
    [Hhat_BS,~,B_BS] = functionChannelEstimates(R_BS,nbrOfRealizations,nbrBSs,K,M,tau_p,pilotIndex,p);   
    

    %Compute SE for the Cell-free mMIMO system with Monte Carlo simulations
    [SE_AP_MR,SE_AP_MMSE] = functionComputeSE_AP_uplink(Hhat_AP,H_AP,R_AP,B_AP,tau_c,tau_p,nbrOfRealizations,N,K,L,p);

    %Compute SE for the Cellular mMIMO system with Monte Carlo simulations
    [SE_BS_MR,SE_BS_RZF,SE_BS_MMMSE] = functionComputeSE_BS_uplink(Hhat_BS,R_BS,B_BS,BSassignment,tau_c,tau_p,nbrOfRealizations,M,K,nbrBSs,p);


    %Save SE values
    SE_BS_MR_tot(:,n) = SE_BS_MR;
    SE_BS_RZF_tot(:,n) = SE_BS_RZF;
    SE_BS_MMMSE_tot(:,n) = SE_BS_MMMSE;
    SE_AP_MR_tot(:,:,n) = SE_AP_MR;
    SE_AP_MMSE_tot(:,:,n) = SE_AP_MMSE;
    
    %Remove large matrices at the end of analyzing this setup
    clear B_AP B_BS D_AP D_BS H_AP Hhat_AP Hhat_BS R_AP R_BS;
    
end


%% Plot simulation results
figure;
hold on; box on;
plot(sort(SE_BS_MMMSE_tot(:)),linspace(0,1,K*nbrOfSetups),'k-','LineWidth',2);
plot(sort(reshape(SE_AP_MMSE_tot(:,4,:),[K*nbrOfSetups 1])),linspace(0,1,K*nbrOfSetups),'b-','LineWidth',2);
plot(sort(reshape(SE_AP_MMSE_tot(:,3,:),[K*nbrOfSetups 1])),linspace(0,1,K*nbrOfSetups),'b--','LineWidth',2);
plot(sort(reshape(SE_AP_MMSE_tot(:,2,:),[K*nbrOfSetups 1])),linspace(0,1,K*nbrOfSetups),'b-.','LineWidth',2);
plot(sort(reshape(SE_AP_MMSE_tot(:,1,:),[K*nbrOfSetups 1])),linspace(0,1,K*nbrOfSetups),'b:','LineWidth',2);
xlabel('Spectral efficiency [bit/s/Hz]','Interpreter','Latex');
ylabel('CDF','Interpreter','Latex');
legend({'Cellular','L4 (MMSE)','L3 (L-MMSE)','L2 (L-MMSE)','L1 (Small cells)'},'Interpreter','Latex','Location','NorthWest');
xlim([0 10]);


figure;
hold on; box on;
plot(sort(SE_BS_MMMSE_tot(:)),linspace(0,1,K*nbrOfSetups),'k-','LineWidth',2);
plot(sort(reshape(SE_AP_MR_tot(:,4,:),[K*nbrOfSetups 1])),linspace(0,1,K*nbrOfSetups),'b-','LineWidth',2);
plot(sort(reshape(SE_AP_MR_tot(:,3,:),[K*nbrOfSetups 1])),linspace(0,1,K*nbrOfSetups),'b--','LineWidth',2);
plot(sort(reshape(SE_AP_MR_tot(:,2,:),[K*nbrOfSetups 1])),linspace(0,1,K*nbrOfSetups),'b-.','LineWidth',2);
plot(sort(reshape(SE_AP_MR_tot(:,1,:),[K*nbrOfSetups 1])),linspace(0,1,K*nbrOfSetups),'b:','LineWidth',2);
xlabel('Spectral efficiency [bit/s/Hz]','Interpreter','Latex');
ylabel('CDF','Interpreter','Latex');
legend({'Cellular','L4 (MR)','L3 (MR)','L2 (MR)','L1 (Small cells)'},'Interpreter','Latex','Location','SouthEast');
xlim([0 10]);
