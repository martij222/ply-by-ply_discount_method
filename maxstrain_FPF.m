function [ maxstrain_SR ] = maxstrain_FPF( local_strain, sig1_T_ult, sig1_C_ult, sig2_T_ult, sig2_C_ult, tau12_ult, E11, E22, G12 )
%maxstrain_FPF Compute strength ratio using maximum strain failure theory
%   local_strain(:,:,1) is a vector of stresses in the local ply coordinate
%   system in the form [sig1 sig2 tau12]. 
%   local_strain(:,:,2) and (:,:,3) are the corresponding stiffnesses and ply numbers, respectively. 
%   The remaining arguments are the ultimate compressive, tensile, and shear
%   strengths of the lamina.

eps1_T_ult = sig1_T_ult/E11; eps1_C_ult = sig1_C_ult/E11;
eps2_T_ult = sig2_T_ult/E22; eps2_C_ult = sig2_C_ult/E22;
gamma12_ult = tau12_ult/G12;

sr1_C = -eps1_C_ult/local_strain(1); sr1_T = eps1_T_ult/local_strain(1); % direction 1
sr2_C = -eps2_C_ult/local_strain(2); sr2_T = eps2_T_ult/local_strain(2); % direction 2
sr12 = gamma12_ult/abs( local_strain(3) ); % shear

temp = [sr1_C sr1_T sr2_C sr2_T sr12];
pos_ratios = temp(temp > 0); % get positive ratios only

maxstrain_SR = min(pos_ratios);

% if maxstrain_SR == sr1_C || maxstrain_SR == sr1_T
%     maxstrain_ply = string( local_strain(1,1,3) ); 
%     if maxstrain_SR == sr1_C
%         maxstrain_failure = '1C';
%     else
%         maxstrain_failure = '1T';
%     end
% elseif maxstrain_SR == sr2_C || maxstrain_SR == sr2_T
%     maxstrain_ply = string( local_strain(2,1,3) );
%     if maxstrain_SR == sr2_C
%         maxstrain_failure = '2C';
%     else
%         maxstrain_failure = '2T';
%     end
% else
%     maxstrain_ply = string( local_strain(3,1,3) );
%     maxstrain_failure = 'S';
% end

end

