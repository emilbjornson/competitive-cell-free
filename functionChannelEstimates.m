function [Hhat,H,B] = functionChannelEstimates(R,nbrOfRealizations,L,K,N,tau_p,pilotIndex,p)
%Generate the channel realizations and estimates of these channels for all
%UEs in the entire network. The channels are modeled as correlated
%Rayleigh fading and the MMSE estimator is used.
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
%R                 = Matrix with dimension N x N x L x K where (:,:,l,k) is
%                    the spatial correlation matrix between AP l and UE k 
%                    in setup n, normalized by the noise power
%nbrOfRealizations = Number of channel realizations
%L                 = Number of APs
%K                 = Number of UEs in the network
%N                 = Number of antennas per AP
%tau_p             = Number of orthogonal pilots
%pilotIndex        = Vector containing the pilot assigned to each UE
%p                 = Uplink transmit power per UE (same for everyone)
%
%OUTPUT:
%Hhat         = Matrix with dimension L*N x nbrOfRealizations x K where
%               (:,n,k) is the estimated collective channel to UE k at
%               channel realization n.
%H            = Matrix with dimension L*N x nbrOfRealizations x K with the
%               true channel realizations. The matrix is organized in the
%               same way as Hhat_MMSE.
%B            = Matrix with dimension N x N x L x K where (:,:,l,j) is the
%               spatial correlation matrix of the estimate between AP l and
%               UE k in setup n, normalized by the noise power


%% Generate channel realizations

%Generate uncorrelated Rayleigh fading channel realizations
H = (randn(L*N,nbrOfRealizations,K)+1i*randn(L*N,nbrOfRealizations,K));


%Go through all channels and apply the spatial correlation matrices
for l = 1:L
    
    for k = 1:K
        
        %Apply correlation to the uncorrelated channel realizations
        Rsqrt = sqrtm(R(:,:,l,k));
        H((l-1)*N+1:l*N,:,k) = sqrt(0.5)*Rsqrt*H((l-1)*N+1:l*N,:,k);
        
    end
    
end


%% Perform channel estimation

%Store identity matrix of size N x N
eyeN = eye(N);

%Generate realizations of normalized noise
Np = sqrt(0.5)*(randn(N,nbrOfRealizations,L,tau_p) + 1i*randn(N,nbrOfRealizations,L,tau_p));


%Prepare to store results
Hhat = zeros(L*N,nbrOfRealizations,K);

if nargout>2
    B = zeros(size(R));
end


%Go through all APs
for l = 1:L
    
    %Go through all pilots
    for t = 1:tau_p
        
        %Compute processed pilot signal for all UEs that use pilot t
        yp = sqrt(p)*tau_p*sum(H((l-1)*N+1:l*N,:,t==pilotIndex),3) + sqrt(tau_p)*Np(:,:,l,t);
        
        %Compute the matrix that is inverted in the MMSE estimator
        PsiInv = (p*tau_p*sum(R(:,:,l,t==pilotIndex),4) + eyeN);
        
        %Go through all UEs that use pilot t
        for k = find(t==pilotIndex)'
            
            %Compute the MMSE estimate
            RPsi = R(:,:,l,k) / PsiInv;
            Hhat((l-1)*N+1:l*N,:,k) = sqrt(p)*RPsi*yp;
            
            %Compute the spatial correlation matrix of the estimate
            if nargout>2
                B(:,:,l,k) = p*tau_p*RPsi*R(:,:,l,k);
            end
            
        end
        
    end
    
end
