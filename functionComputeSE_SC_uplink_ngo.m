function [SE_beta,SE_improved,bestAP,bestAP_improved] = functionComputeSE_SC_uplink_ngo(pmax,p,L,K,tau_p,tau_c,pilotIndex,betaVal,p_improved,pilotIndex_improved)
%Compute the uplink SE for small cell using MR combining, following the
%formulas from the paper and from [15], but with corrections made in
%Proposition 3.
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
%pmax           = Maximum transmit power
%p              = Matrix K x 1 where element k is the  uplink transmit
%                 power of UE k (If it is a scalar, the same value is
%                 used for all users)
%L              = Number of APs
%K              = Total number of UEs
%tau_p          = Length of pilot sequences
%tau_c          = Length of coherence block
%pilotIndex     = Matrix with dimension K x 1 containing the pilot assigned
%                 to the UEs
%betaVal        = Matrix with dimension L x K where element (l,k) contains
%                 the pathloss normalized by the noise variance between
%                 AP l and UE k
%p_imp          = (Optional) Same as p but for the case when each UE is
%                 served by the AP that gives the highest SE
%pilotIndex_imp = (Optional) Same as pilotIndex but for the case when each
%                 UE is served by the AP that gives the highest SE
%
%OUTPUT:
%SE_beta         = K x 1 matrix where the k:th element is the uplink SE of
%                  UE k achieved using small cells using the methods in
%                  "Cell-Free Massive MIMO versus Small Cells"
%SE_improved     = Same as SE_beta but when the AP that gives the
%                  highest SE will serve a UE in the small-cell network
%bestAP          = K x 1 matrix where the k:th element is the AP that UE k
%                  is served by, using the AP assignment method from [15]
%bestAP_improved = K x 1 matrix where the k:th element is the AP that UE k
%                  is served by when the AP giving the highest SE serves
%                  the UE.


%Initialize optional input parameters
if nargin<9
    p_improved = p;
end

if nargin<10
    pilotIndex_improved = pilotIndex;
end

%Prepare to save results
SE_beta = zeros(K,1);
SE_improved = zeros(K,L);

%All APs are available to serve the UEs in the beginning
availableAPs = true(L,1);

%Prepare to store which AP that serves which UE
bestAP = zeros(K,1);

%Go through the UEs in "random" order
for k = 1:K
    
    %Select the AP to serve UE k using (38) in [15]
    activeAPs = find(availableAPs);
    [~,ind] = max(betaVal(activeAPs,k));
    m_k = activeAPs(ind);
    bestAP(k) = m_k;
    
    %Set the AP as non-available
    availableAPs(m_k) = false;
    
    %Compute the SE using the closed form expression in Proposition 3
    A_kl = sum(betaVal(m_k,pilotIndex(k) == pilotIndex).^2 .* p(pilotIndex(k) == pilotIndex)' )/betaVal(m_k,k).^2/p(k) - 1;
    Psi_tkl = (pmax*tau_p*sum(betaVal(m_k,pilotIndex(k) == pilotIndex)) + 1);
    w_mkk = pmax*p(k)*tau_p*(betaVal(m_k,k))^2 ./ ((betaVal(m_k,:)*p)*Psi_tkl - pmax*tau_p*sum(betaVal(m_k,pilotIndex(k) == pilotIndex).^2 .* p(pilotIndex(k) == pilotIndex)') + Psi_tkl);
    
    exponent1 = w_mkk*(1+A_kl);
    exponent2 = w_mkk*A_kl;
    
    SE_beta(k) = (1-tau_p/tau_c)/log(2) * ( exp(1/exponent1)*expint(1/exponent1) - exp(1/exponent2)*expint(1/exponent2) );
    
    %Deal with cases when the closed form expression is unstable
    if isnan(SE_beta(k))
        SE_beta(k) = (1-tau_p/tau_c)/log(2) * ( exp(1/exponent1)*expint(1/exponent1) );
    end
    
    
    %Consider the improved AP association where each UE served by the
    %AP that gives the highest SE
    if nargout > 1
        
        %Go through all APs
        for l = 1:L
            
            %Consider if UE k is served by AP l
            m_k = l;
            
            %Compute SE using the closed form expression in Proposition 3
            A_kl = sum(betaVal(m_k,pilotIndex_improved(k) == pilotIndex_improved).^2 .* p_improved(pilotIndex_improved(k) == pilotIndex_improved)' )/betaVal(m_k,k).^2/p_improved(k) - 1;
            Psi_tkl = (pmax*tau_p*sum(betaVal(m_k,pilotIndex_improved(k) == pilotIndex_improved)) + 1);
            w_mkk = pmax*p_improved(k)*tau_p*(betaVal(m_k,k))^2 ./ ((betaVal(m_k,:)*p_improved)*Psi_tkl - pmax*tau_p*sum(betaVal(m_k,pilotIndex_improved(k) == pilotIndex_improved).^2 .* p_improved(pilotIndex_improved(k) == pilotIndex_improved)') + Psi_tkl);
            
            exponent1 = w_mkk*(1+A_kl);
            exponent2 = w_mkk*A_kl;
            
            SE_improved(k,l) = (1-tau_p/tau_c)/log(2) * ( exp(1/exponent1)*expint(1/exponent1) - exp(1/exponent2)*expint(1/exponent2) );
            
            %Deal with cases when the closed form expression is unstable
            if isnan(SE_improved(k,l))
                
                SE_improved(k,l) = (1-tau_p/tau_c)/log(2) * ( exp(1/exponent1)*expint(1/exponent1) );
                
            end
            
            if isinf(SE_improved(k,l))
                
                SE_improved(k,l) = 0;
                
            end
            
        end
        
    end
    
end


if nargout > 1
    
    %Let each UE be served by the AP giving the highest SE, which is
    %different from the AP selection in [15]
    [SE_improved,bestAP_improved] = max(SE_improved,[],2);
    
end
