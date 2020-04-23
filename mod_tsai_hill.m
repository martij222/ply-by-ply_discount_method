function [ min_mod_tsai_hill_SR, failure_ply ] = mod_tsai_hill( local_stress, sig1_T_ult, sig1_C_ult, sig2_T_ult, sig2_C_ult, tau12_ult )
%mod_tsai_hill Compute strength ratio using Tsai-Hill failure theory
%   local_stress is a vector of stresses in the local ply coordinate
%   system in the form [sig1 sig2 tau12]. The remaining arguments are the ultimate compressive, tensile, and shear
%   strengths of the lamina.

k = size(local_stress, 3);

mod_tsai_hill_SR = zeros([1 k]);
for i=1:k
    % check sig1 and assign X1 accordingly
    if local_stress(1,1,i) >= 0
        X1 = sig1_T_ult;
    else
        X1 = sig1_C_ult;
    end

    % check sig2 and assign X2 accordingly
    if local_stress(2,1,i) >= 0
        X2 = sig1_T_ult;
    else
        X2 = sig1_C_ult;
    end

    % check sig2 and assign Y accordingly
    if local_stress(2,1,i) >= 0
        Y = sig2_T_ult;
    else
        Y = sig2_C_ult;
    end

    S = tau12_ult;

    c = (local_stress(1,1,i)/X1)^2 ...
        - (local_stress(1,1,i)*local_stress(2,1,i)/ X2^2) ...
        + (local_stress(2,1,i)/Y)^2 ...
        + (local_stress(3,1,i)/S)^2;

    mod_tsai_hill_SR(i) = c^-.5;
end
temp = repelem(1:k/2,2); % matrix of ply numbers
[min_mod_tsai_hill_SR, min_i ] = min( mod_tsai_hill_SR );
failure_ply = string( temp(min_i) );

end

