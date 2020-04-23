function [T] = discountmethod(n_trials,theta,F,sd)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
failure_theories = {@maxstress_FPF @tsai_hill_FPF @tsai_wu_H_FPF @tsai_wu_MH_FPF};

NL = length(theta);
discount_method_data = cell(n_trials, 2, length(failure_theories));
ply_thickness = zeros([n_trials NL]); Theta = zeros([n_trials NL]);
E11 = zeros([n_trials NL]); E22 = zeros([n_trials NL]); G12 = zeros([n_trials NL]);
sig1_T_ult = zeros(n_trials,1); sig1_C_ult = zeros(n_trials,1);
sig2_T_ult = zeros(n_trials,1); sig2_C_ult = zeros(n_trials,1);
tau12_ult = zeros(n_trials,1);

%% discount method
for n=1:n_trials
    % generate laminate for current trial
    [ A1, B1, D1, ~, ~, ~, ~, ~, ~, Qbar1, Qbar2, Z_i, Theta_i, E11_i, E22_i, G12_i, ply_thickness_i ] = ...
        generatelaminate( theta, sd(1:6,1) );

    % material strengths
    [ sig1_T_ult_i, sig1_C_ult_i, sig2_T_ult_i, sig2_C_ult_i, tau12_ult_i ] = generaterandomstrengths( sd(7:11) );

    % perform failure analysis for each theory
    for m=1:length(failure_theories)
        f = failure_theories{m}; % select failure theory
        
        % initialize arrays
        SR = zeros([1 2*NL]); mp_strains = zeros([3 1 0]); 
        failure_plies = {}; failure_SRs = [];

        F_i = F; Qbar1_i = Qbar1; A1_i = A1; B1_i = B1; D1_i = D1; % reset all required values
        while any(Qbar1_i, 'all') % loop condition: non-zero Qbars

            % compute initial strains and stresses
            [ global_strain, mp_strain_i ] = computeglobalstrain( F_i, A1_i, B1_i, D1_i, Z_i ); 
            [local_stress, ~] = computelocalstressstrain( global_strain, Qbar1_i, Theta_i, E11_i, E22_i, G12_i );

            % compute strength ratios for the top and bottom of each ply
            for i=1:2*NL
                SR(i) = f( local_stress(:,:,i), sig1_T_ult_i, sig1_C_ult_i, sig2_T_ult_i, sig2_C_ult_i, tau12_ult_i );
            end

            temp=reshape(SR,[2 NL]); % each row represents the SR for one ply (top and bottom)
            SR = min(temp, [], 1); % get SR for each ply

            failure_idx = find(SR <= 1); % get indices of failed plies (if any)

            % check ply conditions and proceed accordingly
            if isempty(failure_idx) % no ply failure
                % get min SR
                [ SR_i, ply ] = min(SR);

                % append failure data
                failure_SRs(end+1) = SR_i; 
                failure_plies{end+1} = ply;
                mp_strains = cat(3,mp_strains,mp_strain_i(1:3));

                % multiply load by SR to get failure load
                F_i = F_i*SR_i;

                % discount damage ply and recompute stiffness matrices
                Qbar1_i = degradeQbar(Qbar1_i, ply);
                [ A1_i, B1_i, D1_i, ~, ~, ~, ~, ~, ~, Qbar1_i, Qbar2 ] = computestiffnessmatrices( Qbar1_i, Qbar2, Z_i );
            else % one or more plies failed
                % get min SR
                SR_i = min(SR);

                % append failure data
                if SR_i > 1
                    failure_SRs(end+1) = SR_i; 
                    failure_plies{end+1} = failure_idx;
                    mp_strains = cat(3,mp_strains,mp_strain_i(1:3));
                else
                    failure_plies{end} = [failure_plies{end} failure_idx];
                end

                % discount damage ply and recompute stiffness matrices
                Qbar1_i = degradeQbar(Qbar1_i, failure_idx);
                [ A1_i, B1_i, D1_i, ~, ~, ~, ~, ~, ~, Qbar1_i, Qbar2 ] = computestiffnessmatrices( Qbar1_i, Qbar2, Z_i );
            end

        end
        
        % store trial data
        discount_method_data{n,1,m} = cumprod(failure_SRs); % the cumulative product of SRs give the actual failure loads
        discount_method_data{n,2,m} = failure_plies; % each cell element is an array of failed ply/plies
    end
    
    % store laminate data before next trial
    ply_thickness(n,:) = ply_thickness_i'; Theta(n,:) = Theta_i';
    E11(n,:) = E11_i'; E22(n,:) = E22_i'; G12(n,:) = G12_i';
    sig1_T_ult(n) = sig1_T_ult_i; sig1_C_ult(n) = sig1_C_ult_i;
    sig2_T_ult(n) = sig2_T_ult_i; sig2_C_ult(n) = sig2_C_ult_i;
    tau12_ult(n) = tau12_ult_i;
end

%% create table
N = repelem(n_trials, n_trials)'; % column to indicate n_trials
SD = repelem(sd', n_trials, 1); % create table of sd values

% @maxstress_FPF @tsai_hill_FPF @tsai_wu_H_FPF @tsai_wu_MH_FPF
maxstress_failure_loads = discount_method_data(:,1,1); maxstress_failure_plies = discount_method_data(:,2,1);
tsai_hill_failure_loads = discount_method_data(:,1,2); tsai_hill_failure_plies = discount_method_data(:,2,2);
tsai_wu_H_failure_loads = discount_method_data(:,1,3); tsai_wu_H_failure_plies = discount_method_data(:,2,3);
tsai_wu_MH_failure_loads = discount_method_data(:,1,4); tsai_wu_MH_failure_plies = discount_method_data(:,2,4);
T = table(N, SD, E11, E22, G12, Theta, ply_thickness, ...
          sig1_T_ult, sig1_C_ult, sig2_T_ult, sig2_C_ult, tau12_ult, ...
          maxstress_failure_loads, maxstress_failure_plies, ...
          tsai_hill_failure_loads, tsai_hill_failure_plies, ...
          tsai_wu_H_failure_loads, tsai_wu_H_failure_plies, ...
          tsai_wu_MH_failure_loads, tsai_wu_MH_failure_plies);

end

