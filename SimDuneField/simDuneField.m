

% INPUTS 
% 1) dx         - dx for sampling (m)
% 2) L_vec      - dune lengths vector (m)
% 3) k_vec      - dune length ratios vector (-)
% 4) Hs_vec     - dune stoss heights vector (m)
% 5) Hl_vec     - dune lee heights vector (m)
% Notes: all units in meters (or unitless for k). Inputs 2-5 all must be
% vectors with identical lengths. 

% Simulate a dune field with one right after the other

function DnFld = simDuneField(dx,L_vec,k_vec,Hs_vec,Hl_vec)

L_vec = L_vec(:); 
k_vec = k_vec(:);
Hs_vec = Hs_vec(:);
Hl_vec = Hl_vec(:);

Ltot = sum(L_vec);      % total length of dune field

x = 0:dx:Ltot; 
z = zeros(size(x));

x_st = 0; % starting x location
z_st = 0; % starting z location
for i = 1:length(L_vec)
    x_end = x_st + L_vec(i);

    [x1,z1] = simOneDune(dx,L_vec(i),k_vec(i),Hs_vec(i),Hl_vec(i));
    z_end = z_st + (z1(end) - z1(1));

    locs = find(x >= x_st & x <= x_end);
    x_l = x(locs);
    x1 = (x1 - x1(1)) + x_st;   % normalize to have new dune starting x at x_st
    z1 = (z1 - z1(1)) + z_st;   % normalize to have new dune starting z at z_st
    z_l = interp1(x1,z1,x_l,"linear");
    z(locs) = z_l;

    % Set up the next loop
    x_st = x_end;
    z_st = z_end;
end

% get locations for dune key points
x_up = [0; cumsum(L_vec(1:end-1))];
x_dn = cumsum(L_vec);
x_pk = x_up + k_vec.*L_vec;
x_mid = 0.5*(x_up + x_dn);

% flip coordinates for all locations
xmax = max(x);

x       = -1*(x - xmax);
x_up    = -1*(x_up - xmax);
x_dn    = -1*(x_dn - xmax);
x_pk    = -1*(x_pk - xmax);
x_mid   = -1*(x_mid - xmax);

x = fliplr(x); x = x(:);
z = fliplr(z); z = z(:);

DnFld = struct; 

DnFld.Dunes.L       = flipud(L_vec);
DnFld.Dunes.k       = flipud(k_vec);
DnFld.Dunes.H       = flipud((0.5*(Hs_vec + Hl_vec)));
DnFld.Dunes.Hs      = flipud(Hs_vec);
DnFld.Dunes.Hl      = flipud(Hl_vec);
DnFld.Dunes.x_up    = flipud(x_up);
DnFld.Dunes.x_dn    = flipud(x_dn);
DnFld.Dunes.x_pk    = flipud(x_pk);
DnFld.Dunes.x_mid   = flipud(x_mid);

DnFld.BEP.x_m = x;
DnFld.BEP.z_m = z;

DnFld = getSplineDet(DnFld); % run spline detrending

end

function DnFld = getSplineDet(DnFld)

x_up    = DnFld.Dunes.x_up;
x_pk    = DnFld.Dunes.x_pk;
x_dn    = DnFld.Dunes.x_dn;
x_mid   = DnFld.Dunes.x_mid;
x_m     = DnFld.BEP.x_m; 
z_m     = DnFld.BEP.z_m; 

z_up = interp1(x_m,z_m,x_up,'spline','extrap');
z_pk = interp1(x_m,z_m,x_pk,'spline','extrap');
z_dn = interp1(x_m,z_m,x_dn,'spline','extrap');

z_baseav = (1/2)*(z_up + z_dn);
z_mean   = (1/2)*(z_baseav + z_pk);
z_s = spline(x_mid,z_mean,DnFld.BEP.x_m);

zd_m = DnFld.BEP.z_m - z_s; 

DnFld.BEP.zd_m = zd_m; 

end

function [x,z] = simOneDune(dx,L,k,Hs,Hl)

if nargin == 4
    Hl = Hs;
end

% Initialize x and z vectors
x = 0:dx:L; % first pass
x = linspace(0,L,5*length(x)); % up the sampling rate to prep for linear interpolation + use linspace to cover entire range
z = zeros(size(x));

locs_s = find(x <= k*L); % find stoss locations
locs_l = find(x >  k*L); % find lee locations

xs = x(locs_s); % pull x values for stoss

% STOSS CALCULATIONS (based on Haque Mahmood 1985)

S_s = Hs/(k*L); % this is direct calculate since I specify H & L

% Parameters that don't vary with x
A = 1/(2*pi*sin(k*pi/2));
C = log(sin(0.5*pi*(1-k)));

% Parameters that vary with x
B = log(sin((pi/L)*(xs + L/2*(1-k))));
D = xs./L;

% Calculate profile based on Haque & Mahmood 1985
z_s = S_s*L*(A*( B - C ) + D);

z(locs_s) = z_s; % assign stoss calcs into overall z vector


% LEE CALCULATIONS

% simple linear profile
x_l = x(locs_l);
z_l = Hl*(ones(size(x_l)) - (x_l-k*L)./(L*(1-k))) + (Hs-Hl);

z(locs_l) = z_l; % assign stoss calcs into overall z vector

% figure; plot(x,z);

end