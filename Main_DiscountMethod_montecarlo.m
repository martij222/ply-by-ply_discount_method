clc;clear;close all;
%% initialize parameters
n_trials = 5000;
theta = [0 90 90 0 90 90 0 0]; %[0 30 -45 45 -30 0];
NL = length(theta);

% loading conditions Nx/y/xy [N/m], Mx/y/xy [Nm/m]
Nx = 1; Ny = 0; Nxy = 0; Mx = 0; My = 0; Mxy = 0;
F = [Nx; Ny; Nxy; Mx; My; Mxy];

% specify sds for each simulation for which to generate data in the following order:
% [E11 E22 G12 v12 dtheta t sig1_T_ult sig1_C_ult sig2_T_ult sig2_C_ult tau12_ult]
a=.1:.1:3;
sd = reshape( repelem(a',11),[11 length(a)] ); % this makes all material property SDs the same

%%
rng(0); % set seed for reproducibility of data
n_sims = size(sd,2);
for sim=1:n_sims
    tic
    fprintf('\nRun ' + string(sim) + ' of ' + string(n_sims) + '\n'); % update progress
    
    % save a new table for this iteration
    tablename = "T_" + string(sim);
    cmd = tablename + " = discountmethod(n_trials,theta,F,sd(:,sim));";
    eval(cmd)

    % get stacking sequence as string
    theta_str = "";
    for j=1:length(theta)
        theta_str = theta_str + string(theta(j))  + '-';
    end
    
    % create file name
    filename = "DM_n_trials_" + string(n_trials) ...
             + "_sd_" + string(sd(1,sim)) ...
             + "_theta_" + theta_str + ".mat";

    % save data in Simulation Data folder
    clearcmd = "clear " + tablename;
    save("Simulation Data\" + filename,tablename); eval(clearcmd) % clear table after saving
    toc
end
%% plot a few failure visualizations
% idx = randperm(n_trials);
% for i=1:3
%     plotfailureinfo(maxstress_failure_loads{i},maxstress_failure_plies{i})
% end