function [SE_MR,SE_RZF,SE_MMMSE] = functionComputeSE_BS_uplink(Hhat,R,B,BSassignment,tau_c,tau_p,nbrOfRealizations,M,K,nbrBSs,p)
%Compute uplink SE for different receive combining schemes using 
%Theorem 4.1 in 
%
%Emil Bjornson, Jakob Hoydis and Luca Sanguinetti (2017), 
%"Massive MIMO Networks: Spectral, Energy, and Hardware Efficiency", 
%Foundations and Trends in Signal Processing: Vol. 11, No. 3-4, 
%pp. 154-655. DOI: 10.1561/2000000093.
%
%INPUT:
%Hhat              = Matrix with dimension M*nbrBSs x nbrOfRealizations x K
%                    where (:,n,k) is the estimated collective channel from
%                    all BSs to UE k at channel realization n.
%R                 = Matrix with dimension M x M x nbrBSs x K where
%                    (:,:,l,k) is the spatial correlation matrix between BS
%                    l and UE k in setup n, normalized by the noise power
%B                 = Matrix with dimension M x M x nbrBSs x K where
%                    (:,:,l,k) is the spatial correlation matrix of the
%                    estimate between BS l and UE k in setup n, normalized
%                    by the noise power
%BSassignment      = Matrix with dimension K x 1 containing the
%                    index of the BS that serves a particular UE
%tau_c             = Length of the coherence block
%tau_p             = Length of pilot sequences and number of UEs per cell
%nbrOfRealizations = Number of channel realizations
%M                 = Number of antennas per BS
%K                 = Total number of UEs
%nbrBSs            = Number of BSs and cells
%p                 = Uplink transmit power per UE (same for everyone)
%
%OUTPUT:
%SE_MR    = K x 1 matrix where the k:th element is the uplink SE of UE k
%           achieved with MR combining 
%SE_RZF   = Same as SE_MR but with RZF combining
%SE_MMMSE = Same as SE_MR but with M-MMSE combining



%Store identity matrices of different sizes
eyetaup = eye(tau_p);
eyeM = eye(M);

%Compute the pre-log factor (normalized with number of channel realizations)
%assuming only uplink transmission
prelogFactor = (tau_c-tau_p)/(tau_c*nbrOfRealizations);

%Prepare to store simulation results
SE_MR = zeros(K,1);

if nargout > 1
    SE_RZF = zeros(K,1);
end

if nargout > 2
    SE_MMMSE = zeros(K,1);
    
    %Compute sum of the estimation error correlation matrices at every BS
    C_tot = p*sum(R-B,4);
    
end



%% Go through all channel realizations
for n = 1:nbrOfRealizations
    
    %Go through all BSs
    for l = 1:nbrBSs
        
        %Extract channel estimate realizations from all UEs to BS l
        Hhatallj = reshape(Hhat(1+(l-1)*M:l*M,n,:),[M K]);
        
        %Extract which UEs are served by the BS
        servedUEs = find(BSassignment==l);
        
        
        %Compute MR combining
        V_MR = Hhatallj(:,servedUEs);
        
        if nargout > 1 %Compute RZF combining
            V_RZF = p*V_MR/(p*(V_MR'*V_MR)+eyetaup);
        end
        
        if nargout > 2 %Compute M-MMSE combining
            V_MMMSE = p*(p*(Hhatallj*Hhatallj')+C_tot(:,:,l)+eyeM)\V_MR;
        end
        
        
        
        %Go through all UEs in cell l
        for ind = 1:tau_p
            
            %Extract UE index
            k = servedUEs(ind);
            
            
            %%MR combining
            v = V_MR(:,ind); %Extract combining vector
            
            %Compute numerator and denominator of the effective SINR
            numerator = p*abs(v'*Hhatallj(:,k))^2;
            denominator = p*norm(v'*Hhatallj)^2 + v'*(C_tot(:,:,l)+eyeM)*v - numerator;
            
            %Compute instantaneous SE for one channel realization
            SE_MR(k) = SE_MR(k) + prelogFactor*real(log2(1+numerator/denominator));
            
            
            %%RZF combining
            if nargout > 1
                
                v = V_RZF(:,ind); %Extract combining vector
                
                %Compute numerator and denominator of the effective SINR
                numerator = p*abs(v'*Hhatallj(:,k))^2;
                denominator = p*norm(v'*Hhatallj)^2 + v'*(C_tot(:,:,l)+eyeM)*v - numerator;
                
                %Compute instantaneous SE for one channel realization
                SE_RZF(k) = SE_RZF(k) + prelogFactor*real(log2(1+numerator/denominator));
                
            end
            
            
            %%M-MMSE combining
            if nargout > 2
                
                v = V_MMMSE(:,ind); %Extract combining vector
                
                %Compute numerator and denominator of the effective SINR
                numerator = p*abs(v'*Hhatallj(:,k))^2;
                denominator = p*norm(v'*Hhatallj)^2 + v'*(C_tot(:,:,l)+eyeM)*v - numerator;
                
                %Compute instantaneous SE for one channel realization
                SE_MMMSE(k) = SE_MMMSE(k) + prelogFactor*real(log2(1+numerator/denominator));
                
            end
            
        end
        
    end
    
end
