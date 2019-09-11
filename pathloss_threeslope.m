function PL = pathloss_threeslope(dvec)
%Compute the pathloss according to the three-slope model in (52) and (53)
%of [15] using the parameters specified in that paper.
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
%dvec = Vector of arbitrarily length containing distances (in meters)
%       between UEs and APs
%
%OUTPUT:
%PL   = Pathloss without shadowing, for each of the elements of dvec


%Set distances in the three-slope model
d0 = 10; %meter
d1 = 50; %meter

%Constant term in the model from (53) in [15]
L = 140.7151;

%Compute the pathloss using the three-slope model in (52) in [15]
PL = zeros(length(dvec),1);

for ind = 1:length(dvec)
    
    d = dvec(ind);
    
    if d<=d0
        PL(ind) = -L -15*log10(d1/1000) -20*log10(d0/1000);
    elseif d<=d1
        PL(ind) = -L -15*log10(d1/1000) -20*log10(d/1000);
    else
        PL(ind) = -L -35*log10(d/1000);
    end
    
end
