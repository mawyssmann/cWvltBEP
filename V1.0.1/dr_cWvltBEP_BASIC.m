
% NOTES TO READER: 
% The best way to use this driver script is to run it
% section by section using the "Run Section" button while having the cursor
% within the desired section. If you have your own x and z BEP data, you
% can load it with the x_m and z_m naming convention manually or else by
% replacing the lines in Step 1 to assign them. I recommend copying this
% driver script and renaming it and then referencing the data that you
% would like to process.

% The 'amplitude-normalized' wavelet is often referenced as the L1 norm wavelet and
% the 'energy normalized' wavelet is often referenced as the L2 norm
% wavelet. 

%% Step 1: Load the data, initialize object, and do basic wavelet calculations

% ---- INPUTS ----
x_m = DnFld.BEP.x_m;    % Load vector of x locations
z_m = DnFld.BEP.zd_m;   % Load vector of detrended z bed elevation profile
% ----------------

cWv = cWvltBEP(x_m,z_m); % Initialize
cWv = getWt(cWv);

%% Step 2: Plot and visualize some of the data

% ---- INPUTS ----
x_ref = DnFld.Dunes.x_pk; % Note: make these empty if you do not want to compare with reference data... x_ref = [];
L_ref = DnFld.Dunes.L;
H_ref = DnFld.Dunes.H;
% ----------------

% Plot the L1 amplitude with P/P95 superimposed
figure; 
plotver = 3;  % 1 = energy-normalized (L2 norm) wavelet; 2 = energy-normalized power ratio relative to noise; 3 = amplitude-normalized (L1 norm) wavelet
log2C = 0; 
plot95 = 1;
plotWt(cWv,plotver,log2C,plot95)
hold on;
plot(x_ref,log2(L_ref),'ro','LineWidth',2)

% Good plot for figure
plotver = 3;
log2C = 0;
plot95 = 1;
xlims = [];
ylims = [];
clims = [];
plotWtandTS(cWv,plotver,log2C,plot95,xlims,ylims,clims)

%% Step 3: Lwx and Hwx local comparison figures

% INPUTS -----
PeriodRng   = []; % This instuitutes a manual trim of period (lambda) values contributing to wavelet-based dune scale calculations
ls          = 5*mean(DnFld.Dunes.L(:));  % this is the smoothing lengthscale, which should be about 5 times the average dune length.
% ------------

% Wavelet method coefficients (Calibrated in Wyssmann JHE paper)
PR95_th = 1; 
mpow = 3;
CLm     = 0.84; 
CL_CIs  = [0.49 1.09];
CHm     = 1.55;
CH_CIs  = [1.14 2.07];

cWv = getLwxHwx(cWv,PeriodRng,PR95_th,mpow,ls,CLm,CL_CIs,CHm,CH_CIs);

plot_EvalLwxHwx_CIs(cWv,x_ref,L_ref,H_ref,ls)