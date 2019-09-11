%%This Matlab script generates Figure 4 in the paper:
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

%Number of Monte Carlo setups
nbrOfSetups = 200;

%Number of channel realizations per setup
nbrOfRealizations = 1000;

%Number of APs in the cell-free network
L = 100;

%Number of UEs
K = 40;

%Number of antennas per AP
N = 1;

%Length of the coherence block
tau_c = 200;

%Number of pilots per coherence block
tau_p = 20;

%Uplink transmit power per UE (mW)
p = 100;


%Prepare to save simulation results
SE_AP_MR_tot = zeros(K,4,nbrOfSetups);
SE_AP_MR_exact_tot = zeros(K,3,nbrOfSetups);
SE_AP_MMSE_tot = zeros(K,4,nbrOfSetups);
SE_AP2_MR_tot = zeros(K,4,nbrOfSetups);
SE_AP2_MR_exact_tot = zeros(K,3,nbrOfSetups);
SE_AP2_MMSE_tot = zeros(K,4,nbrOfSetups);
APselectionSC = zeros(nbrOfSetups,1);
APselectionSC2 = zeros(nbrOfSetups,1);


%% Go through all setups
for n = 1:nbrOfSetups
    
    %Display simulation progress
    disp(['Setup ' num2str(n) ' out of ' num2str(nbrOfSetups)]);
    
    %Generate one setup with UEs at random locations
    [gainOverNoisedB,R,pilotIndexCF,pilotIndexSC] = generateSetup_threeslope(L,K,N,tau_p,1,p);
    betaVal = db2pow(gainOverNoisedB);
    
    
    %Full transmit power case
    
    %Generate channel realizations, channel estimates, and estimation
    %error correlation matrices for all UEs to the APs 
    [Hhat_AP,H_AP,B_AP] = functionChannelEstimates(R,nbrOfRealizations,L,K,N,tau_p,pilotIndexCF,p);

    %Compute SE using Monte Carlo simulations
    [SE_AP_MR,SE_AP_MMSE] = functionComputeSE_AP_uplink(Hhat_AP,H_AP,R,B_AP,tau_c,tau_p,nbrOfRealizations,N,K,L,p);
    
    %Compute SE using closed-form expressions
    SE_CF = functionComputeSE_CF_uplink_ngo(p,p*ones(K,1),L,K,tau_p,tau_c,pilotIndexCF,betaVal);
    [SE_SC_beta,SE_SC_improved,bestAP,bestAP_improved] = functionComputeSE_SC_uplink_ngo(p,p*ones(K,1),L,K,tau_p,tau_c,pilotIndexSC,betaVal);

    
    %Compute the fraction of UEs that connect to the same AP as in [15]
    APselectionSC(n) = sum(bestAP==bestAP_improved)/K;
    
    %Save SE values
    SE_AP_MR_tot(:,:,n) = SE_AP_MR; %MR combining with Monte Carlo simulations
    SE_AP_MMSE_tot(:,:,n) = SE_AP_MMSE; %MMSE/L-MMSE combining with Monte Carlo simulations
    SE_AP_MR_exact_tot(:,1,n) = SE_SC_beta; %Small cells as in [15]
    SE_AP_MR_exact_tot(:,2,n) = SE_SC_improved; %Small cells as in [15], but with max-SE AP allocation
    SE_AP_MR_exact_tot(:,3,n) = SE_CF; %Level 2 using closed-form expression from [15]
    
    
    
    %Max-min fairness transmit power case
    
    %Extract terms in the numerator and denominator of the SINRs
    [signalCF,interferenceCF,signalSC,interferenceSC,signalSC2,interferenceSC2] = functionSINRterms_uplink_ngo(p,L,K,tau_p,pilotIndexCF,pilotIndexSC,betaVal);
    
    %Optimize the transmit powers for max-min fairness - Cell-free mMIMO
    [~,pBest_CF] = functionPowerOptimization_maxmin(signalCF,interferenceCF,p,1);
    
    %Optimize the transmit powers for small cells - largest large-scale
    %fading association
    [~,pBest_SC_beta] = functionPowerOptimization_maxmin(signalSC,interferenceSC,p,1);
    
    %Optimize the transmit powers for small cells - largest SE association
    [~,pBest_SC_improved] = functionPowerOptimization_maxmin(signalSC2,interferenceSC2,p,1);
    
    
    %Compute SE using Monte Carlo simulations
    [SE_AP_MR,SE_AP_MMSE] = functionComputeSE_AP_uplink(Hhat_AP,H_AP,R,B_AP,tau_c,tau_p,nbrOfRealizations,N,K,L,pBest_CF,pBest_SC_improved);
    
    %Compute SE using closed-form expressions
    SE_CF = functionComputeSE_CF_uplink_ngo(p,pBest_CF,L,K,tau_p,tau_c,pilotIndexCF,betaVal);
    [SE_SC_beta,SE_SC_improved,bestAP,bestAP_improved] = functionComputeSE_SC_uplink_ngo(p,pBest_SC_beta,L,K,tau_p,tau_c,pilotIndexSC,betaVal,pBest_SC_improved);

    
    %Compute the fraction of UEs that connect to the same AP as in [15]
    APselectionSC2(n) = sum(bestAP==bestAP_improved)/K;
    
    %Save SE values
    SE_AP2_MR_tot(:,:,n) = SE_AP_MR; %MR combining with Monte Carlo simulations
    SE_AP2_MMSE_tot(:,:,n) = SE_AP_MMSE; %MMSE/L-MMSE combining with Monte Carlo simulations
    SE_AP2_MR_exact_tot(:,1,n) = SE_SC_beta; %Small cells as in [15]
    SE_AP2_MR_exact_tot(:,2,n) = SE_SC_improved; %Small cells as in [15], but with max-SE AP allocation
    SE_AP2_MR_exact_tot(:,3,n) = SE_CF; %Level 2 using closed-form expression from [15]
    
    %Remove large matrices at the end of analyzing this setup
    clear B_AP D_AP H_AP Hhat_AP;
    
end


%% Plot simulation results
figure;
hold on; box on;

plot(sort(reshape(SE_AP_MMSE_tot(:,2,:),[K*nbrOfSetups 1])),linspace(0,1,K*nbrOfSetups),'k-.','LineWidth',2);
plot(sort(reshape(SE_AP_MR_tot(:,1,:),[K*nbrOfSetups 1])),linspace(0,1,K*nbrOfSetups),'r-','LineWidth',2);
plot(sort(reshape(SE_AP_MR_exact_tot(:,1,:),[K*nbrOfSetups 1])),linspace(0,1,K*nbrOfSetups),'r-','LineWidth',5);
plot(sort(reshape(SE_AP_MR_exact_tot(:,2,:),[K*nbrOfSetups 1])),linspace(0,1,K*nbrOfSetups),'r:','LineWidth',2);
plot(sort(reshape(SE_AP_MR_exact_tot(:,3,:),[K*nbrOfSetups 1])),linspace(0,1,K*nbrOfSetups),'b-.','LineWidth',5);

xlabel('Spectral efficiency [bit/s/Hz]','Interpreter','Latex');
ylabel('CDF','Interpreter','Latex');
legend({'L2 (L-MMSE)','L1 (Small cells)','Ref.~[15] (Small cells)','Ref.~[15] (Improved)','Ref.~[15] (L2 MR)'},'Interpreter','Latex','Location','NorthWest');
xlim([0 2]);


figure;
hold on; box on;

plot(sort(reshape(SE_AP_MMSE_tot(:,2,:),[K*nbrOfSetups 1])),linspace(0,1,K*nbrOfSetups),'k-.','LineWidth',2);
plot(sort(reshape(SE_AP2_MR_tot(:,1,:),[K*nbrOfSetups 1])),linspace(0,1,K*nbrOfSetups),'r-','LineWidth',2);
plot(sort(reshape(SE_AP2_MR_exact_tot(:,1,:),[K*nbrOfSetups 1])),linspace(0,1,K*nbrOfSetups),'r-','LineWidth',5);
plot(sort(reshape(SE_AP2_MR_exact_tot(:,2,:),[K*nbrOfSetups 1])),linspace(0,1,K*nbrOfSetups),'r:','LineWidth',2);
plot(sort(reshape(SE_AP2_MR_exact_tot(:,3,:),[K*nbrOfSetups 1])),linspace(0,1,K*nbrOfSetups),'b-.','LineWidth',5);

xlabel('Spectral efficiency [bit/s/Hz]','Interpreter','Latex');
ylabel('CDF','Interpreter','Latex');
legend({'L2 (L-MMSE), Full power','L1 (Small cells)','Ref.~[15] (Small cells)','Ref.~[15] (Improved)','Ref.~[15] (L2 MR)'},'Interpreter','Latex','Location','SouthEast');
xlim([0 2]);


disp(['Fraction of UEs prefering max-beta AP (correlated): ' num2str(mean(APselectionSC))]);
disp(['Fraction of UEs prefering max-beta AP (uncorrelated): ' num2str(mean(APselectionSC2))]);
