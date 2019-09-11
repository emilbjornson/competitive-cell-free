function SE_CF = functionComputeSE_CF_uplink_ngo(Pmax,p,L,K,tau_p,tau_c,pilotIndex,betaVal)
%Compute the uplink SE for Cell-free mMIMO using MR combining, following
%the formulas from the paper and from [15]:
%
%H. Q. Ngo, A. Ashikhmin, H. Yang, E. G. Larsson, and T. L. Marzetta,
%"Cell-Free Massive MIMO versus Small Cells," IEEE Trans. Wireless Commun.,
%vol. 16, no. 3, pp. 1834-1850, 2017.
%
%All references to equations in the code refers to that paper.
%
%
%This function was developed as a part of the paper:
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
%
%INPUT:
%Pmax          = Maximum transmit power
%p             = Matrix K x 1 where element k is the  uplink transmit
%                power of UE k (If it is a scalar, the same value is
%                used for all users)
%L             = Number of APs
%K             = Total number of UEs
%tau_p         = Length of pilot sequences
%tau_c         = Length of coherence block
%pilotIndex    = Matrix with dimension K x 1 containing the pilot assigned
%                to the UEs
%betaVal       = Matrix with dimension L x K where element (l,k) contains
%                the pathloss normalized by the noise power between AP l
%                and UE k
%
%OUTPUT:
%SE_CF         = K x 1 matrix where the k:th element is the uplink SE of
%                UE k achieved with MR combining



%Compute the gamma parameters for all APs and UE according to (8) in [15]
gammaVal = zeros(L,K);

for l = 1:L
    for k = 1:K
        gammaVal(l,k) = Pmax*tau_p*(betaVal(l,k))^2 ./ (Pmax*tau_p*sum(betaVal(l,pilotIndex(k) == pilotIndex)) + 1);
    end
end

%Compute the numerator of (27) in [15] for all UEs, assuming full power
signal = p.*((sum(gammaVal,1)).^2)';


%Compute the denominator of (27) in [15] for all UEs, assuming full power
interference = zeros(K,1);

%Create matrix to easier multipy with transmit powers
prep = repmat(p',[L 1]);

%Go through all UEs
for k = 1:K
    
    %Add the last two terms of the denominator
    interference(k) = interference(k) + sum(gammaVal(:,k).*sum(prep.*betaVal,2)) + sum(gammaVal(:,k));
    
    %Compute the first term with interference due to pilot contamination
    coPilot = (pilotIndex(k) == pilotIndex); %Extract UEs that use same pilot
    coPilot(k) = false; %Remove the UE of interest
    samePilot = find(coPilot);
    
    for ind = 1:length(samePilot)
        interference(k) = interference(k) + p(samePilot(ind))*sum(gammaVal(:,k).*betaVal(:,samePilot(ind))./betaVal(:,k))^2;
    end
    
end


%Compute the SE with Cell-free mMIMO and MR using (27) in [15], with a
%pre-log factor accounting for uplink pilot
SE_CF = (1-tau_p/tau_c)*log2(1+signal./(interference));
