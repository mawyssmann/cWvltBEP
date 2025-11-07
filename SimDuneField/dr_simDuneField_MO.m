
% ARMA model for L (based on MO River), rest as random variables

% ----- INPUTS -----
Nsim    = 5000; % # dunes to simulate
dx      = 0.45; % BEP sample spacing in streamwise direction
% NOTE: Need to load ARIMAmdl_MORivJune_p8q6.mat to have the estimated ARMA
% model for the MO River to run this script section
% ------------------

phat2 = [5.3 3.7];          % k beta
phat3 = [8.62 0.011];       % C_H gamma
phat4 = [0 0.2627];         % (Hs-H)/H normal

Y       = simulate(EstMdl,Nsim); % simulation based on ARMA model (fitted to MO River)
L_vec   = exp(Y);
k_vec   = betarnd(phat2(1),phat2(2),Nsim,1);
a_H     = gamrnd(phat3(1),phat3(2),Nsim,1);
HR      = normrnd(phat4(1),phat4(2),Nsim,1);

H = a_H.*L_vec.^0.81;
Hs_vec = HR.*H + H; 
Hl_vec = 2*H - Hs_vec; 

% Simulate dune field
DnFld = simDuneField(dx,L_vec,k_vec,Hs_vec,Hl_vec);

% Do a quick visualization of the BEP (use zd_m, which is the detrended
% fluctuating BEP version)
figure; plot(DnFld.BEP.x_m,DnFld.BEP.zd_m,'-k');