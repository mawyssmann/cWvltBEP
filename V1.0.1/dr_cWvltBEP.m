

%% pull elevation profile from MoData structure and convert all units to meters

% --- INPUTS ---
elem = 3; 
loc  = 'C';
% --------------

switch loc
    case 'C'
        x_m = MoData_int(elem).Center.x*1609.34;
        z_m = MoData_int(elem).Center.z*0.3048;
    case 'L'
        x_m = MoData_int(elem).Left.x*1609.34;
        z_m = MoData_int(elem).Left.z*0.3048;
    case 'R'
        x_m = MoData_int(elem).Right.x*1609.34;
        z_m = MoData_int(elem).Right.z*0.3048;
end



%% Initialize object and do basic wavelet calculations

cWv = cWvltBEP(x_m,z_md); % Initialize
cWv = getWt(cWv);

% Plot the L1 amplitude with P/P95 superimposed
figure; 
plotver = 3; 
log2C = 0; 
plot95 = 1;
plotWt(cWv,plotver,log2C,plot95)

%% (OPTIONAL) Overlay true data as red circles in the wavelet plot

x_true = x_peak;
L_true = L;

hold on;
plot(x_true,log2(L),'ro','LineWidth',2)

%% (Optional) Trim wavelet plot limits (x,y,c)

xlims = [1.9e5 2.1e5];
ylims = [1 512];
clims = [];

cropWt(cWv,xlims,ylims,clims)

%% (Optional) Plot wavelet with time series over top

% plotver = 3; 
% log2C = 0; 
% plot95 = 1;
% xlims = [2.025e5 2.045e5];
% ylims = [1 512];

% Good plot for figure
plotver = 3; 
log2C = 0; 
plot95 = 1;
% xlims = [2*10^5 2.04*10^5];
% xlims = [5.7*10^4 6*10^4];
% xlims = [2.05e4 2.2e4];
xlims = [];
ylims = [];
clims = [0 80];
% clims = [];

plotWtandTS(cWv,plotver,log2C,plot95,xlims,ylims,clims)

%% Get length scale

% V1 with weighted average period based on statistically significant
npow = 1; 
PeriodRng = [1 128];
Lwx = getLwx_WgtAvPeriod(cWv,npow,PeriodRng);

figure; 
plot(cWv.x,Lwx,'-'); hold on; 
plot(DnFld.Dunes.x_mid,DnFld.Dunes.L,'.');

%% DEVELOPMENT BELOW

%% Sine wave testing v4

% Inputs
lambda  = 10;
amp     = 10; 
dx      = 1;
N       = 10000;
Nv      = 12;

% Calculate BEP series
x      = 0:dx:(N*dx);
z_m  = amp*sin(2*pi*x/lambda);

% Call functions from class cWvltBEP
cWv = cWvltBEP(x,z_m);  % initialize
cWv = getWt(cWv,Nv);    % call wavelet

% % get normalized amplitude, plot coefficient and get max
max(cWv.A_L1(:))/amp;


%% Call a voice scan

% Do a voice scan
Nvmin = 2;
Nvmax = 24; 
[Nvs,varRats] = doNvScan(cWv,Nvmin,Nvmax);
    % Doesn't seem to be as important as I saw the other night. I'm
    % wondering if I had my scales messed up in that test (lambda versus
    % length of series, or something similar). 
    % Observation: Seems that you just need ENOUGH voices to capture all
    % features. Variance output stabilizes after a threshold value
    % (typically > 6-10, so the 12 value from Grinsted seems reasonable).
    % TC98 mainly argue that the tradeoff with Nv is computational time,
    % which is much better now than in 1998.


%%

figure; 
plot(z_m);
xlim([0 0.2e5])

figure;
cwt(z_m);
xlim([0 0.2e5])
