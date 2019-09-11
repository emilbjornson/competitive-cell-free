function [R_AP,R_BS,pilotIndex,BSassignment,gainOverNoisedB_AP,gainOverNoisedB_BS,UEpositions,APpositions,BSpositions] = generateSetup(L,K,N,M,nbrOfSetups)
%Generate realizations of the simulation setup described in Section IV.
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
%L                  = Number of APs for the Cell-free system
%K                  = Number of UEs in the network
%N                  = Number of antennas per AP
%M                  = Number of antennas per cellular BS
%nbrOfSetups        = Number of setups with random UE locations
%
%OUTPUT:
%R_AP               = Matrix with dimension N x N x L x K x nbrOfSetups
%                     where (:,:,l,k,n) is the spatial correlation matrix
%                     between AP l and UE k in setup n, normalized by noise
%                     power
%R_BS               = Matrix with dimension M x M x 4 x K x nbrOfSetups
%                     where (:,:,l,k,n) is the spatial correlation matrix
%                     between BS l and UE k in setup n, normalized by noise
%                     power
%pilotIndex         = Matrix with dimension K x nbrOfSetups containing the
%                     pilot assigned to the UEs in different setups
%BSassignment       = Matrix with dimension K x nbrOfSetups containing the
%                     index of the BS that serves a particular UE
%gainOverNoisedB_AP = Matrix with dimension L x K x nbrOfSetups where
%                     element (l,k,n) is the channel gain (normalized by
%                     the noise power) between AP l and UE k in setup n
%gainOverNoisedB_BS = Matrix with dimension 4 x K x nbrOfSetups where
%                     element (l,k,n) is the channel gain (normalized by
%                     the noise power) between BS l and UE k in setup n
%UEpositions        = Vector of length K with UE positions, where the real
%                     part is the horizontal position and the imaginary
%                     part is the vertical position
%APpositions        = Vector of length L with the AP locations, measured in
%                     the same way as UEpositions
%BSpositions        = Vector of length L with the AP locations, measured in
%                     the same way as UEpositions



%% Define simulation setup

%Size of the coverage area (as a square with wrap-around)
squareLength = 1000; %meter

%Communication bandwidth
B = 20e6;

%Noise figure (in dB)
noiseFigure = 5;

%Compute noise power
noiseVariancedBm = -174 + 10*log10(B) + noiseFigure;

%Pathloss parameters
alpha = 36.7;
constantTerm = -30.5;

%Standard deviation of the shadow fading
sigma_sf = 4;

%Decorrelation distance of the shadow fading
decorr = 9;

%Height difference between an AP and a UE
distanceVertical = 10;

%Define the antenna spacing (in number of wavelengths)
antennaSpacing = 1/2; %Half wavelength distance

%Angular standard deviation around the nominal angle (measured in degrees)
ASDdeg = 15;


%Number of cellular BSs
nbrBSs = 4;

%Number of pilots is equal to the number of UEs per cell
tau_p = K/nbrBSs;


%Number of cellular BSs per dimension on the grid
nbrBSsPerDim = sqrt(nbrBSs);

%Distance between BSs in vertical/horizontal direction
interBSDistance = squareLength/nbrBSsPerDim;

%Deploy BSs on the grid
locationsGridHorizontal = repmat(interBSDistance/2:interBSDistance:squareLength-interBSDistance/2,[nbrBSsPerDim 1]);
locationsGridVertical = locationsGridHorizontal';
BSpositions = locationsGridHorizontal(:) + 1i*locationsGridVertical(:);


%Number of APs per dimension on the grid
nbrAPsPerDim = sqrt(L);

%Distance between APs in vertical/horizontal direction
interAPDistance = squareLength/nbrAPsPerDim;

%Deploy APs on the grid
locationsGridHorizontal = repmat(interAPDistance/2:interAPDistance:squareLength-interAPDistance/2,[nbrAPsPerDim 1]);
locationsGridVertical = locationsGridHorizontal';
APpositions = locationsGridHorizontal(:) + 1i*locationsGridVertical(:);


%Prepare to save results
gainOverNoisedB_AP = zeros(L,K,nbrOfSetups);
gainOverNoisedB_BS = zeros(nbrBSs,K,nbrOfSetups);
R_BS = zeros(M,M,nbrBSs,K,nbrOfSetups);
R_AP = zeros(N,N,L,K,nbrOfSetups);
pilotIndex = zeros(K,nbrOfSetups);
BSassignment = zeros(K,nbrOfSetups);


%% Go through all setups
for n = 1:nbrOfSetups
    
    %Compute alternative BS locations by using wrap around
    wrapHorizontal = repmat([-squareLength 0 squareLength],[3 1]);
    wrapVertical = wrapHorizontal';
    wrapLocations = wrapHorizontal(:)' + 1i*wrapVertical(:)';
    BSpositionsWrapped = repmat(BSpositions,[1 length(wrapLocations)]) + repmat(wrapLocations,[nbrBSs 1]);
    
    
    %Random UE locations with uniform distribution with K/nbrBSs in the
    %coverage area of each cellular BS
    nbrOfUEsPerBS = zeros(nbrBSs,1);
    
    %Prepare to compute UE locations
    UEpositions = zeros(K,1);
    
    %Prepare to store shadowing correlation matrix
    shadowCorrMatrix = sigma_sf^2*ones(K,K);
    shadowBSrealizations = zeros(K,nbrBSs);
    
    nbrOfUEs = 0;
    
    %Add UEs until the each BS has tau_p UEs to serve
    while min(nbrOfUEsPerBS)<tau_p
        
        %Generate a random UE location in the area
        UEposition = (rand(1,1) + 1i*rand(1,1)) * squareLength;
        
        %Compute distances from the UE to each of the BSs
        [distanceBSstoUE,whichpos] = min(abs(BSpositionsWrapped - repmat(UEposition,size(BSpositionsWrapped))),[],2);
        distances = sqrt(distanceVertical^2+distanceBSstoUE.^2);
        
        
        %If this is not the first UE
        if nbrOfUEs>0
            
            %Compute distances from the new prospective UE to all other UEs
            shortestDistances = zeros(nbrOfUEs,1);
            
            for i = 1:nbrOfUEs
                shortestDistances(i) = min(abs(UEposition - UEpositions(i) + wrapLocations));
            end
            
            %Compute conditional mean and standard deviation necessary to
            %obtain the new shadow fading realizations, when the previous
            %UEs' shadow fading realization have already been generated.
            %This computation is based on Theorem 10.2 in "Fundamentals of
            %Statistical Signal Processing: Estimation Theory" by S. Kay
            newcolumn = sigma_sf^2*2.^(-shortestDistances/decorr);
            term1 = newcolumn'/shadowCorrMatrix(1:nbrOfUEs,1:nbrOfUEs);
            meanvalues = term1*shadowBSrealizations(1:nbrOfUEs,:);
            stdvalue = sqrt(sigma_sf^2 - term1*newcolumn);
            
        else %If this is the first UE
            
            %Add the UE and begin to store shadow fading correlation values
            meanvalues = 0;
            stdvalue = sigma_sf;
            newcolumn = [];
            
        end
        
        %Generate the shadow fading realizations
        shadowing = meanvalues + stdvalue*randn(1,nbrBSs);
        
        %Compute the channel gain (in dB)
        gains = constantTerm - alpha*log10(distances) + shadowing';
        
        %Find which BS the UE would like to connect to
        [~,bestBS] = max(gains);
        
        %If the BS doesn't have tau_p UEs yet
        if nbrOfUEsPerBS(bestBS)<tau_p
            
            %Add the UE to the preferred cell
            nbrOfUEsPerBS(bestBS) = nbrOfUEsPerBS(bestBS) + 1;
            
            %Compute and store the UE index
            nbrOfUEs = sum(nbrOfUEsPerBS);
            k = nbrOfUEs;
            
            %Update shadowing correlation matrix and store realizations
            shadowCorrMatrix(1:nbrOfUEs-1,nbrOfUEs) = newcolumn;
            shadowCorrMatrix(nbrOfUEs,1:nbrOfUEs-1) = newcolumn';
            shadowBSrealizations(nbrOfUEs,:) = shadowing;
            
            %Assign the UE to the BS
            BSassignment(k,n) = bestBS;
            
            %Store the UE position
            UEpositions(k) = UEposition;
            
            %Assign a pilot "randomly" to the UE (the first unused one)
            pilotIndex(k,n) = nbrOfUEsPerBS(bestBS);
            
            %Store the channel gain divided by the noise power (in dB)
            gainOverNoisedB_BS(:,k,n) = gains - noiseVariancedBm;
            
            %Go through all BSs
            for l = 1:nbrBSs
                
                %Compute nominal angle between the new UE k and AP l
                angletoUE = angle(UEpositions(k)-BSpositionsWrapped(l,whichpos(l)));
                
                %Generate normalized spatial correlation matrix using the
                %local scattering model
                R_BS(:,:,l,k,n) = db2pow(gainOverNoisedB_BS(l,k,n))*functionRlocalscattering(M,angletoUE,ASDdeg,antennaSpacing);
                
            end
            
        end
        
        
    end
    
    
    
    %Compute alternative AP locations by using wrap around
    wrapHorizontal = repmat([-squareLength 0 squareLength],[3 1]);
    wrapVertical = wrapHorizontal';
    wrapLocations = wrapHorizontal(:)' + 1i*wrapVertical(:)';
    APpositionsWrapped = repmat(APpositions,[1 length(wrapLocations)]) + repmat(wrapLocations,[L 1]);

    %Create correlated shadow fading realizations
    shadowAPrealizations = sqrtm(shadowCorrMatrix)*randn(K,L);
    
    
    %Go through all UEs
    for k = 1:K
        
        %Compute distances assuming for the UE to the APs
        [distanceAPstoUE,whichpos] = min(abs(APpositionsWrapped - repmat(UEpositions(k),size(APpositionsWrapped))),[],2);
        distances = sqrt(distanceVertical^2+distanceAPstoUE.^2);
        
        %Compute the channel gain divided by the noise power (in dB)
        gainOverNoisedB_AP(:,k,n) = constantTerm - alpha*log10(distances) + shadowAPrealizations(k,:)' - noiseVariancedBm;
        
        %Go through all APs
        for l = 1:L
            
            %Compute nominal angle between UE k and AP l
            angletoUE = angle(UEpositions(k)-APpositionsWrapped(l,whichpos(l)));
            
            %Generate normalized spatial correlation matrix using the local
            %scattering model
            R_AP(:,:,l,k,n) = db2pow(gainOverNoisedB_AP(l,k,n))*functionRlocalscattering(N,angletoUE,ASDdeg,antennaSpacing);
            
        end
        
    end
    
end
