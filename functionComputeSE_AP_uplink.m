function [SE_MR,SE_MMSE,sumSE_SIC] = functionComputeSE_AP_uplink(Hhat,H,R,B,tau_c,tau_p,nbrOfRealizations,N,K,L,p,p1)
%Compute uplink SE for Cell-free mMIMO for the four different receiver
%cooperation levels, using either MR or MMSE/L-MMSE combining
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
%Hhat              = Matrix with dimension N*L x nbrOfRealizations x K
%                    where (:,n,k) is the estimated collective channel from
%                    all BSs to UE k at channel realization n.
%H                 = Matrix with dimension N*L x nbrOfRealizations x K
%                    where (:,n,k) is the true collective channel from all
%                    BSs to UE k at channel realization n.
%R                 = Matrix with dimension N x N x L x K where
%                    (:,:,l,k) is the spatial correlation matrix between BS
%                    l and UE k in setup n, normalized by the noise power
%B                 = Matrix with dimension N x N x L x K where
%                    (:,:,l,k) is the spatial correlation matrix of the
%                    estimate between BS l and UE k in setup n, normalized
%                    by the noise power
%tau_c             = Length of coherence block
%tau_p             = Length of pilot sequences and number of UEs per cell
%nbrOfRealizations = Number of channel realizations
%N                 = Number of antennas per AP
%K                 = Total number of UEs
%L                 = Number of APs
%p                 = Matrix K x 1 where element k is the uplink transmit
%                    power of UE k (If it is a scalar, the same value is
%                    used for all users)
%p1                = (Optional) Same as p but only for Level 1
%
%OUTPUT:
%SE_MR     = K x 4 matrix where the (k,n):th element is the uplink SE of 
%            UE k achieved with MR combining at cooperation level n
%SE_MMMSE  = Same as SE_MR but with MMSE or L-MMSE combining
%sumSE_SIC = Scalar with the sum SE achieved using MMSE-SIC combining



%If only one transmit power is provided, use the same for all the UEs
if length(p) == 1
   p = p*ones(K,1);
end

%If no specific Level 1 transmit powers are provided, use the same as for
%the other levels
if nargin<12
    p1 = p;
end


%Store identity matrices of different sizes
eyeN = eye(N);
eyeK = eye(K);
eyeNL = eye(N*L);

%Compute the prelog factor
prelogFactor = (1-tau_p/tau_c);

%Prepare to store simulation results
SE_MR = zeros(K,4);
SE_MMSE = zeros(K,4);

%Compute sum of all estimation error correlation matrices at every BS
C_tot = zeros(N,N,L);
C1_tot = zeros(N,N,L);
for k = 1:K
    C_tot = C_tot + p(k)*(R(:,:,:,k)-B(:,:,:,k));
    C1_tot = C1_tot + p1(k)*(R(:,:,:,k)-B(:,:,:,k));
end

C_tot_blk = zeros(N*L,N*L);

for l = 1:L
    C_tot_blk(1+(l-1)*N:l*N,1+(l-1)*N:l*N) = C_tot(:,:,l);
end

%Diagonal matrix with transmit powers and its square root
Dp = diag(p);
Dp12 = diag(sqrt(p));
Dp112 = diag(sqrt(p1));


%Prepare to save simulation results
SE_MR_level1 = zeros(K,L);
SE_MMSE_level1 = zeros(K,L);

signal_MR_level23 = zeros(L,K);
scaling_MR_level23 = zeros(L,K);
Gp_MR_level23 = zeros(L,L,K);

signal_MMSE_level23 = zeros(L,K);
scaling_MMSE_level23 = zeros(L,K);
G_MMSE_level23 = zeros(L,L,K);

sumSE_SIC = 0;


%% Go through all channel realizations
for n = 1:nbrOfRealizations
    
    
    %Level 4
    
    %Extract channel estimate realizations from all UEs to all APs
    Hhatallj = reshape(Hhat(:,n,:),[N*L K]);
    
    %Compute MR combining
    V_MR = Hhatallj;
    
    %Compute MMSE combining
    V_MMSE = ((Hhatallj*Dp*Hhatallj')+C_tot_blk+eyeNL)\(Hhatallj*Dp);
    

    %Go through all UEs
    for k = 1:K
        
        %%MR combining
        v = V_MR(:,k); %Extract combining vector
        
        %Compute numerator and denominator of instantaneous SINR at Level 4
        numerator = p(k)*abs(v'*Hhatallj(:,k))^2;
        denominator = norm(v'*Hhatallj*Dp12)^2 + v'*(C_tot_blk+eyeNL)*v - numerator;
        
        %Compute instantaneous SE for one channel realization
        SE_MR(k,4) = SE_MR(k,4) + prelogFactor*real(log2(1+numerator/denominator))/nbrOfRealizations;
        
        
        %%MMSE combining
        v = V_MMSE(:,k); %Extract combining vector
        
        %Compute numerator and denominator of instantaneous SINR at Level 4
        numerator = p(k)*abs(v'*Hhatallj(:,k))^2;
        denominator = norm(v'*Hhatallj*Dp12)^2 + v'*(C_tot_blk+eyeNL)*v - numerator;
        
        %Compute instantaneous SE for one channel realization
        SE_MMSE(k,4) = SE_MMSE(k,4) + prelogFactor*real(log2(1+numerator/denominator))/nbrOfRealizations;
        
    end
    
    
    %Compute sum SE with MMSE-SIC combining at Level 4
    if nargout > 2

        sumSE_SIC = sumSE_SIC + prelogFactor*real(log2(det(eyeK + Dp*V_MR'*((C_tot_blk+eyeNL)\V_MR))))/nbrOfRealizations;
        
    end
    
    
    
    
    %Levels 1-3
    gp_MR_level23 = zeros(L,K,K);
    gp_MMSE_level23 = zeros(L,K,K);
    
    
    %Go through all APs
    for l = 1:L
        
        %Extract channel realizations from all UEs to AP l
        Hallj = reshape(H(1+(l-1)*N:l*N,n,:),[N K]);
        
        %Extract channel estimate realizations from all UEs to AP l
        Hhatallj = reshape(Hhat(1+(l-1)*N:l*N,n,:),[N K]);
        
        
        %Compute MR combining
        V_MR = Hhatallj;
        V_MMSE = ((Hhatallj*Dp*Hhatallj')+C_tot(:,:,l)+eyeN)\(V_MR*Dp);
        

        %Go through all UEs
        for k = 1:K
            
            
            %%MR combining
            v = V_MR(:,k); %Extract combining vector
            
            
            %Level 2 and Level 3
            signal_MR_level23(l,k) = signal_MR_level23(l,k) + (v'*Hallj(:,k))/nbrOfRealizations;
            gp_MR_level23(l,:,k) = gp_MR_level23(l,:,k) + (v'*Hallj)*Dp12;
            scaling_MR_level23(l,k) = scaling_MR_level23(l,k) + norm(v).^2/nbrOfRealizations;
            
            
            %Level 1
            
            %Compute numerator and denominator of instantaneous SINR
            numerator = p1(k)*abs(v'*Hhatallj(:,k))^2;
            denominator = norm((v'*Hhatallj)*Dp112)^2 + v'*(C1_tot(:,:,l)+eyeN)*v - numerator;
            
            %Compute instantaneous SE for one channel realization
            SE_MR_level1(k,l) = SE_MR_level1(k,l) + prelogFactor*real(log2(1+numerator/denominator))/nbrOfRealizations;
            
            
            
            %%MMSE combining
            v = V_MMSE(:,k); %Extract combining vector
            
            
            %Level 2 and Level 3
            signal_MMSE_level23(l,k) = signal_MMSE_level23(l,k) + (v'*Hallj(:,k))/nbrOfRealizations;
            gp_MMSE_level23(l,:,k) = gp_MMSE_level23(l,:,k) + (v'*Hallj)*Dp12;
            scaling_MMSE_level23(l,k) = scaling_MMSE_level23(l,k) + norm(v).^2/nbrOfRealizations;
            
            
            %Level 1
            
            %Compute numerator and denominator of instantaneous SINR
            numerator = p1(k)*abs(v'*Hhatallj(:,k))^2;
            denominator = norm((v'*Hhatallj)*Dp112)^2 + v'*(C1_tot(:,:,l)+eyeN)*v - numerator;
            
            %Compute instantaneous SE for one channel realization
            SE_MMSE_level1(k,l) = SE_MMSE_level1(k,l) + prelogFactor*real(log2(1+numerator/denominator))/nbrOfRealizations;
            
            
        end
        
    end
    
    %Compute averaging of terms for Level 2 and Level 3
    for k = 1:K
        
        Gp_MR_level23(:,:,k) = Gp_MR_level23(:,:,k) + gp_MR_level23(:,:,k)*gp_MR_level23(:,:,k)'/nbrOfRealizations;
        G_MMSE_level23(:,:,k) = G_MMSE_level23(:,:,k) + gp_MMSE_level23(:,:,k)*gp_MMSE_level23(:,:,k)'/nbrOfRealizations;
        
    end
    
end


%Compute SE for Level 1
SE_MR(:,1) = max(SE_MR_level1,[],2);
SE_MMSE(:,1) = max(SE_MMSE_level1,[],2);


%Compute SE for Level 2 and Level 3
for k = 1:K
    
    %With MR combining
    b = signal_MR_level23(:,k);
    A = Gp_MR_level23(:,:,k) + diag(scaling_MR_level23(:,k)) - p(k)*(b*b');
    SE_MR(k,3) = prelogFactor*real(log2(1+p(k)*b'*(A\b)));   
    SE_MR(k,2) = prelogFactor*real(log2(1+p(k)*abs(mean(b)).^2 / mean(mean(A))));   

    %With L-MMSE combining
    b = signal_MMSE_level23(:,k);
    A = G_MMSE_level23(:,:,k) + diag(scaling_MMSE_level23(:,k)) - p(k)*(b*b');
    SE_MMSE(k,3) = prelogFactor*real(log2(1+p(k)*b'*(A\b)));  
    SE_MMSE(k,2) = prelogFactor*real(log2(1+p(k)*abs(mean(b)).^2 / mean(mean(A))));   

end
