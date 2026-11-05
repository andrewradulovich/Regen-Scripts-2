%% sizeRegen2 
% %% Regenerative Cooling Channel Thermal Analysis
% % Authors: Andrew Radulovich (aradulov@purdue.edu), Alex Suppiah, Jan Ayala, Kamon Blong, Zach Hodgdon
% % Based on PSP AC code Size Regen by, Andrew Radulovich, Alex Suppiah,
% Jan Ayala, Zach Hodgdon, Adam Grendy. 
% % First Created: 09/12/2026
% % Last Updated: 

   %{ 
    Description:

    Program Methodology:

    Inputs:

    Outputs: 
    - 
    
    Assumptions:
    - Steady state
    - No backside wall heat transfer
    - Equally distributed temperature inside channels
    - Wicking heat into fuel doesn't change bulk gas temperature

   %}

clc; clear; close all
addpath("539 CEA Wrapper\");
addpath("IPA Data\");
addpath("Material Data\");
addpath(("539 CEA Wrapper\CEA Calculation Files"))
%Data = LoadData();
u = convertUnits;
%% Reuse CEA? 
New_CEA = false;
fprintf("If CEA values are still stored in workspace from prev run & no changes were made to engine parameters or steps, dont need to rerun CEA. \n");
input1 = input("Run CEA? Type y / n. \n",'s');
if input1 == "y" || input1 == "Y"
    New_CEA = true;
end
%% Simulation Inputs
resolution_chamber = 100; % Number of points to evaluate before the nozzle
resolution_converging = 100; % Number of points to evaluate in the converging section of the nozzle
resolution_diverging = 100; % Number of points to evaluate in the diverging section of the nozzle

heatflux_factor = 1;
fin_efficiency_factor = 1;
coolant_heatflux_factor = 1;

latconductance = 1;
axialconductance = 1;
axialgeointegral = 0;
temperaturefieldaxial = 1;
temperaturefieldlateral = 1;

tempfieldpoints = 21; % Use odd numbers


%% Engine Inputs 

% Propellants
P_c = 250; 
OF = 1.2;
fuel = 0; % Choose fuel propellant
oxidizer = 0; % Choose Ox propellant 
fuel_weight = 100;%[75,25]; % weight of each fuel used (if only 1 fuel, set to 100)  
fuel_temp = 293.15;%[293.15,293.15]; % [K] Initial gues for fuel temperature for CEA calculation
oxidizer_temp = 90.17; % [K] Initial guess for ox temperature for CEA calculation
oxidizer_weight = 100; %

% Coolant Parameters
mdot_coolant = 3 * u.LB2KG; %lbm/s
coolant = 'methanol'; % Choose coolant fluid
%coolant = 'water';
%coolant = "isopropyl alcohol"
coolant_direction = 0;
inlet_temp = 293.15; % Initial Coolant Temperature [K]
inlet_pressure = 500; % Initial Coolant Pressure [psi]

% Chamber/Nozzle Geometry
R_t = 0.73217; % Throat radius [in]
exp_ratio = 2.9525;
con_ratio = 6.558;
conv_angle = 37.5;
conical_half_angle = 15; 
L_star = 45; 
R_fillets = [1.5,1.5,1.5];

% Channel Geometry
num_channels = 62;
% chamber, throat, nozzle end [in]
t_w = [.040 .040 .030 .040]; % Thin wall thickness 
h_c = [.080 .080 .070 .090]; % Channel height
w_c = [.080 .080 .020 .050]; % Channel width (mid arclength)
floatbool = 0;
t_w_float = .035;
h_c_float = .075;
w_c_float = .040;
AR_float = 3;

% Material & Manufacturing Method
material = 0; % 0: 6061-RAM2, 1: Inconel718, 2: GrCop42, 3: Inconel625, 4: C101/OFHC Copper, 5: AlSi10Mg

%% Material Properties
if material == 0
    mat = materialProps('al6061-ram2');
elseif material == 1 
    mat = materialProps('inco-718');
elseif material == 2
    mat = materialProps('grcop42');
elseif material == 3
    mat = materialProps('inco-625');
elseif material == 4
    mat = materialProps('c101');
elseif material == 5
    mat = materialProps('alsi10mg');
end
k_w = mat.k;
K_ad = mat.K_ad;
yts = mat.yts;
CTE = mat.CTE;
elong = mat.elong;
E = mat.E;
epsilon_emax = @(x) yts(x)*1000000 / E(x);




nu = 0.3; % poissons ratio (guess)
%e = 24 * 0.001; % surface roughness (mm) [micrometer*conversion]
roughness_table = readmatrix(pwd + "/Material Data/surface_roughness.xlsx",'Range','A12:E16');
e = [roughness_table(2,2), roughness_table(5,2)] .* 0.001; %Surface roughness (mm) [micrometer*conversion] [45, 90]



%% Initialization 

% Fuel Propellant Definition 
if fuel == 0 
    fuel_def = 'C3H8O,2propanol';
end
if oxidizer == 0
    oxidizer_def = 'O2(L)';
end

% Units
R_t_SI = R_t *u.IN2M; % [m]
%R_fillets_SI = R_fillets .* u.IN2M; % [m];
L_star_SI = L_star * u.IN2M; % [m]
t_w_SI = t_w .* u.IN2M; % [m]
h_c_SI = h_c .* u.IN2M; % [m]
w_c_SI = w_c .* u.IN2M; % [m]
t_w_float_SI = t_w_float * u.IN2M;
h_c_float_SI = h_c_float * u.IN2M;
w_c_float_SI = w_c_float * u.IN2M;





%% Generate Chamber / Nozzle Contour
% Throat channel geo
%w_c_t = w_c(2) ; 
%h_c_t = h_c(2) ; 
%t_w_t = t_w(2); 

% Throat Diameter & Radius of Curvature
D_t = 2* R_t * u.IN2M; % Throat diameter [m]
R_of_curve_up = R_fillets(2) * D_t/2;
R_of_curve_down = R_fillets(3) * D_t/2;
R_of_curve =  (R_of_curve_up + R_of_curve_down) / 2; % Avg Radius of  curvature of the throat [m] **This is what is used in the Huzel & Huang example problems**

% Generate the geometry functions
geom = engineContour2(R_t_SI, exp_ratio, con_ratio, conv_angle, conical_half_angle, L_star_SI, R_fillets);

% Define the bounds of each cell [m]
x_chamber_bounds = linspace(geom.bounds.chamber(1),geom.bounds.chamber(2),resolution_chamber+1);
x_converging_bounds = linspace(geom.bounds.converging(1),geom.bounds.converging(2),resolution_converging+1);
x_diverging_bounds = linspace(geom.bounds.diverging(1),geom.bounds.diverging(2),resolution_diverging+1);
x_bounds = [x_chamber_bounds, x_converging_bounds(2:end), x_diverging_bounds(2:end)];
r_bounds = geom.r(x_bounds); 
drdx_bounds = geom.drdx(x_bounds);
dsdx_bounds = (1 + (drdx_bounds).^2).^(1/2);
dxds_bounds = dsdx_bounds.^(-1);
kappa_f = geom.kappa(x_bounds);

% Get channel functions
geom_channels = channelCountour(geom, h_c_SI,w_c_SI,t_w_SI,floatbool,h_c_float_SI,w_c_float_SI,t_w_float_SI);


% Define the midpoints of each cell [m]
j = 1; 
x_midpoints = ones(1, length(x_bounds)-1);
ds = nan(1,length(x_bounds)-1);
ds_ch = nan(1,length(x_bounds)-1);
%r_midpoints = ones(1, length(x_bounds)-1);
%drdx_midpoints = ones(1,length(x_bounds)-1);
while j < length(x_bounds)
    x1 = x_bounds(j); % Front Boundary of Cell
    x2 = x_bounds(j+1); % Back Boundary of Cell
    ds(j) = integral(@(x) sqrt(1 + geom.drdx(x).^2), x1, x2)    ; % Arc Length of Cell
    midpointResidual = @(xm) integral(@(x) sqrt(1 + geom.drdx(x).^2),x1,xm) - ds(j)/2; % Equation whose root is the arc-length midpoint
    x_mid = fzero(midpointResidual, [x1, x2]);
    %r_mid = geom.r(x_mid);
    %drdx_mid = geom.drdx(x_mid);

    x_midpoints(j) = x_mid;
    %r_midpoints(j) = r_mid;
    %drdx_midpoints(j) = drdx_mid; % [in / in]

    % Channel-centerline arc length
    ds_ch(j) = integral(@(x) ...
        (1 + geom.kappa(x).* ...
        (geom_channels.tw(x) + geom_channels.hc(x)./2)) .* ...
        sqrt(1 + geom.drdx(x).^2), ...
        x1,x2);
    j = j+1;
end
ds_m = nan(1,length(x_midpoints)-1);
r_midpoints = geom.r(x_midpoints);
drdx_midpoints = geom.drdx(x_midpoints);
kappa = geom.kappa(x_midpoints);
dsdx_midpoints = (1 + (drdx_midpoints).^2).^(1/2);
%dxds_midpoints = dsdx_midpoints.^(-1);
dxds_midpoints = geom.dxds(x_midpoints);

% Generate Subsonic and Supersonic area ratios at each midpoint
A = pi() .* r_midpoints.^2; % Combustion area [m^2]
A_t = pi() * R_t_SI^2; % Throat area [m^2]
sub = A(resolution_chamber+1:(resolution_chamber+resolution_converging)) ./ A_t;
sup = A(resolution_chamber+resolution_converging+1:end) ./ A_t;
A_ratio = [sub sup];

% Calculate Arched channel geometry
t_w = geom_channels.tw(x_midpoints); % wall thickness at midpoints [m]
h_c = geom_channels.hc(x_midpoints); % channel heights at midpoints [m]
w_c = geom_channels.wc(x_midpoints); % channel widths at midpoints [m]
t_w_f = geom_channels.tw(x_bounds); % wall thickness at faces [m]
h_c_f = geom_channels.hc(x_bounds); % Channel heights at faces [m]
w_c_f = geom_channels.wc(x_bounds);% channel widthds at faces [m]
dt_w_fdx = geom_channels.dtwdx(x_bounds); % wall thickness derivative at faces [m/m]




r_wg = r_midpoints; % radius at hot gas wall [m]
r_wg_f = r_bounds; % radius at hot gass wall faces [m]
%r_wl = r_wg + t_w; % radius at liquid wall [in]
rho_wl = r_wg + dxds_midpoints .* t_w; % Corrected radius at liquid wall based on normal orientation [m]
%r_wl_f = r_wg_f + t_w_f; % radius at liquids wall faces [in]
rho_wl_f = r_wg_f + dxds_bounds .* t_w_f; % Corrected radius at liquid wall faces
%r_backside = r_wl + h_c; % backside radius of channel [in]
rho_backside = r_wg + dxds_midpoints .* (h_c + t_w); % backside radius of channel with normal adjustment [m];
%r_backside_f = rwl_f + h_c_f; % backsid radius of channel at faces [in]
rho_backside_f = r_wg_f + dxds_bounds .* (h_c_f + t_w_f); % normal adjusted backside radius of channel at faces [m]
%r_midchan = 0.5* (r_wl + r_backside); % mid radius of channel [in]
rho_midchan = r_wg + dxds_midpoints .* ( (h_c/2) + t_w) ; % mid radius of channel with normal adjustment [m]
%r_midchan_f = 0.5*(r_wl_f + r_backside_f); % mid radius of channel at faces [in]
rho_midchan_f = r_wg_f + dxds_bounds .* ( (h_c_f/2) + t_w_f); % mid radius of channel at faces with normal adjustment [m]
theta_channel = w_c ./rho_midchan; % clocking angle of channel width [rad]
theta_channel_f = w_c_f ./ rho_midchan_f;
rho_m = r_wg + dxds_midpoints .* (t_w/2); %Adjusted radius at midpoint of the wall
%ds_ch = (1 + kappa .* (t_w + h_c./2)) .* ds; % Arc length of channel [m]




%A_channel = theta_channel/2  * (rho_backside.^2 - rho_wl.^2); % area of channel [in^2]
A_channel = w_c .* h_c; 
%A_channel_f = theta_channel_f/2 * (rho_backside_f.^2 - rho_wl_f.^2); % area of channel at faces [in^2]
A_channel_f = w_c_f .* h_c_f; 
theta_p = 2*pi()/ num_channels; % Pitch angle of channel/land pattern [rad]
theta_land = theta_p - theta_channel; % Pitch angle of land width [rad]
theta_land_f = theta_p - theta_channel_f;
%r_m = (r_wg + r_wl)/2; % radius of the midpoint of the thin wall
w_land = rho_wl .* theta_land; % width of the base of the land
w_land_f = rho_wl_f .* theta_land_f;
AR_channel = h_c ./ w_c;
AR_fin = h_c ./ w_land;



mdot_channel = mdot_coolant / num_channels; % mdot of coolant per channel [kg/s]
%P_wet = (r_wl + r_backside).*theta_channel + 2*(r_backside - r_wl); %wetted perimeter of channel
P_wet = 2 * (w_c + h_c);
P_wet_f = 2 * (w_c_f + h_c_f);
%P_wet_f = (r_wl_f + r_backside_f).*theta_channel_f + 2*(r_backside_f - r_wl_f); %wetted perimeter of channel at faces
D_hydraulic = 4 * A_channel ./ P_wet; % hydralic diameter [m^2]
D_hydraulic_f = 4 * A_channel_f ./ P_wet_f;



% Calculate Aeff for axial conduction
nz = tempfieldpoints;
neta = tempfieldpoints;
zeta = linspace(0,1,nz);
eta = linspace(-0.5,0.5,neta);
v = 1;
dxds = dxds_midpoints;
dxds_f = dxds_bounds;
RgAC = nan(1,length(x_midpoints)-1);
RgAL = nan(1,length(x_midpoints)-1);
RgAC_midpoints = nan(1,length(x_midpoints)-1);
RgAL_midpoints = nan(1,length(x_midpoints)-1);
RgAC_integral = nan(1,length(x_midpoints)-1);
RgAL_integral = nan(1,length(x_midpoints)-1);
ds_midpoints = nan(1,length(x_midpoints)-1);
HAC = nan(length(x_midpoints)-1,nz);
HAL = nan(length(x_midpoints)-1,nz);
while v < length(x_midpoints)

    ds_midpoints(v) = (ds(v)+ds(v+1))/2; % Bad misleading name
    if abs(kappa_f(v+1)*t_w_f(v+1)) < 1e-9
        Is = r_wg_f(v+1)*t_w_f(v+1) + (dxds_f(v+1)/2)*(t_w_f(v+1)^2);
    else 
        Is = dxds_f(v+1)*t_w_f(v+1)/kappa_f(v+1) + ((r_wg_f(v+1)*kappa_f(v+1)-dxds_f(v+1))/((kappa_f(v+1)^2)))*log(1+kappa_f(v+1)*t_w_f(v+1));                  
    end
    Aeff_c_mid = Is * theta_channel_f(v+1);
    Aeff_l_mid = Is * theta_land_f(v+1);
    RgAC_midpoints(v) = ds_midpoints(v)/Aeff_c_mid;
    RgAL_midpoints(v) = ds_midpoints(v)/Aeff_l_mid;

    % Distributed axial geometry kernel for temperature-field model
    rho_z = r_wg_f(v+1) + dxds_f(v+1).*zeta.*t_w_f(v+1);
    J_z =   1 + kappa_f(v+1).*zeta.*t_w_f(v+1);
    HAC(v,:) = theta_channel_f(v+1).*t_w_f(v+1) .* rho_z ./ J_z ./ ds_midpoints(v);
    HAL(v,:) = theta_land_f(v+1).*t_w_f(v+1) .* rho_z ./ J_z ./ ds_midpoints(v);

    if axialgeointegral
        
        x1 = x_midpoints(v);
        x2 = x_midpoints(v+1);
        
        integrand_c = @(x) ((1+geom.drdx(x).^2).^(1/2))./ arrayfun(@(xx) getAeffChannel(xx,geom,geom_channels),x);
        integrand_l = @(x) ((1+geom.drdx(x).^2).^(1/2))./ arrayfun(@(xx) getAeffLand(xx,geom,geom_channels,theta_p),x);
        RgAC_integral(v) = integral(integrand_c,x1,x2);
        RgAL_integral(v) = integral(integrand_l,x1,x2);
        %RgAC(v) = RgAC_integral(v);
        %RgAL(v) = RgAL_integral(v);
        if temperaturefieldaxial
            for jz = 1:nz
                 zeta_jz = zeta(jz);
    
                 invRg_c = nan(size(eta));
                 invRg_l = nan(size(eta));
                
                 for je = 1:length(eta)
    
                     eta_j = eta(je);
    
                     Rg_c = integral(@(x) ...
                         axialRgIntegrandChannel( ...
                             x,eta_j,zeta_jz,geom,geom_channels), ...
                         x1,x2);
    
                     Rg_l = integral(@(x) ...
                         axialRgIntegrandLand( ...
                             x,eta_j,zeta_jz,geom,geom_channels,theta_p), ...
                         x1,x2);
    
                     invRg_c(je) = 1/Rg_c;
                     invRg_l(je) = 1/Rg_l;
    
                 end
    
                 HAC(v,jz) = trapz(eta,invRg_c);
                 HAL(v,jz) = trapz(eta,invRg_l);
            end
        end
    end    
    v =v+1;
end
if axialgeointegral
    RgAC = RgAC_integral;
    RgAL = RgAL_integral;
    RgAC_difference = (RgAC - RgAC_midpoints) ./ RgAC_midpoints;
    RgAL_difference = (RgAL - RgAL_midpoints) ./ RgAL_midpoints;
else 
    RgAC = RgAC_midpoints;
    RgAL = RgAL_midpoints;
end

% Construct Clat for lateral conduction with temperature field
%nz = tempfieldpoints;
%zeta = linspace(0,1,nz);

rho = nan(length(x_midpoints), nz);
Jg = nan(length(x_midpoints), nz);
Clat = nan(length(x_midpoints), nz);

for n = 1:length(x_midpoints)
    rho = r_wg(n) + dxds(n).*zeta.*t_w(n);

    Jg = 1 + kappa(n).*zeta.*t_w(n);    

    Clat(n,:) = (2*t_w(n)/(theta_p/2)) .* Jg./rho;
end

%% Get Combustion Properties
% Use CEA for this 
if New_CEA
    % Subsonic Area ratios
    [CEAdata_sub] = callCEA('fr', P_c, 'psi', 'o/f', OF, 'sub', sub, fuel_def, 'K', fuel_temp, fuel_weight, oxidizer_def, 'K', oxidizer_temp, oxidizer_weight);

    % Supersonic Area ratios 
    [CEAdata_sup] = callCEA('fr', P_c, 'psi', 'o/f', OF, 'sup', sup, fuel_def, 'K', fuel_temp, fuel_weight, oxidizer_def, 'K', oxidizer_temp, oxidizer_weight);

    % Save CEA data to .mat file
    save("RegenCEA_Data.mat","CEAdata_sub","CEAdata_sup");
else 
    RegenCEA_Data = load("RegenCEA_Data.mat");
    CEAdata_sub = RegenCEA_Data.CEAdata_sub;
    CEAdata_sup = RegenCEA_Data.CEAdata_sup;
end


% Collect CEA results 

% Hot Gas Temperature [K]
T_g_sub = squeeze(CEAdata_sub('t'));
T_g_sup = squeeze(CEAdata_sup('t'));
T_g = transpose([T_g_sub(2:end); T_g_sup(2:end)]);
T_g = [T_g(1)*ones(1,resolution_chamber), T_g];
T_g_tot = T_g_sup(1);

% hot gas specific heat pressure []
cp_g_sub = squeeze(CEAdata_sub('cp'));
cp_g_sup = squeeze(CEAdata_sup('cp'));
cp_g = transpose([cp_g_sub(2:end); cp_g_sup(2:end)]);
cp_g = [cp_g(1)*ones(1,resolution_chamber), cp_g];
cp_g_tot = cp_g_sup(1);

% hot gas thermal conductivity  []
k_g_sub = squeeze(CEAdata_sub('k'));
k_g_sup = squeeze(CEAdata_sup('k'));
k_g = transpose([k_g_sub(2:end); k_g_sup(2:end)]);
k_g = [k_g(1)*ones(1,resolution_chamber), k_g];
k_g_tot = k_g_sup(1);

% Hot Gas Mach Number []
M_sub = squeeze(CEAdata_sub('mach'));
M_sup = squeeze(CEAdata_sup('mach'));
M = transpose([M_sub(2:end); M_sup(2:end)]);
M = [M(1)*ones(1,resolution_chamber), M];
M_tot = M_sup(1);

% Hot Gas Specific heat ratio [
gamma_sub = squeeze(CEAdata_sub('gammas'));
gamma_sup = squeeze(CEAdata_sup('gammas'));
gamma = transpose([gamma_sub(2:end); gamma_sup(2:end)]);
gamma = [gamma(1)*ones(1,resolution_chamber), gamma];
gamma_tot = gamma_sup(1);

% Hot Gas Viscosity [
mu_g_sub = squeeze(CEAdata_sub('visc'));
mu_g_sup = squeeze(CEAdata_sup('visc'));
mu_g = transpose([mu_g_sub(2:end); mu_g_sup(2:end)]);
mu_g = [mu_g(1)*ones(1,resolution_chamber), mu_g];
mu_g_tot = mu_g_sup(1); 

% Hot gas prantl number 
Pr_g_sub = squeeze(CEAdata_sub('prandtl'));
Pr_g_sup = squeeze(CEAdata_sup('prandtl'));
Pr_g = transpose([Pr_g_sub(2:end); Pr_g_sup(2:end)]);
Pr_g = [Pr_g(1)*ones(1,resolution_chamber), Pr_g];
Pr_g_tot = Pr_g_sup(1);

% Hot gas pressure
P_g_sub = squeeze(CEAdata_sub('p'));
P_g_sup = squeeze(CEAdata_sup('p'));
P_g = transpose([P_g_sub(2:end); P_g_sup(2:end)]);
P_g = [P_g(1)*ones(1,resolution_chamber), P_g];
P_g_tot = P_g_sup(1);


% cstar 
cstar_sup = squeeze(CEAdata_sup('cstar'));
cstar = cstar_sup(1);

% Append chamber section area ratios to beginning of vector
A_ratio = [ones(1,resolution_chamber).* A_ratio(1), A_ratio];


% Get Combustion Data a throat
%[~, ~, ~, M_t, gamma_t, P_g_t, T_g_t, rho_g_t, mu_g_t, Pr_g_t, ~, ~, son_g_t, cp_g_t] = RunCEA(P_c, P_e, fuel, fuel_weight, fuel_temp, oxidizer, oxidizer_temp, OF_tot, 0, 1, 2, 0, 0, CEA_input_name);
% [CEAdata] = callCEA('fr', P_c, 'psi', 'o/f', OF, 'sup', 1, fuel_def, 'K', fuel_temp, fuel_weight, oxidizer_def, 'K', oxidizer_temp, oxidizer_weight);
% T_g_out = squeeze(CEAdata('t')); 
% T_g_t = T_g_out(2);
% T_g_tot = T_g_out(1);
% cp_g_out = squeeze(CEAdata('cp')); 
% cp_g_t = cp_g_out(2);
% cp_g_tot = cp_g_out(1);
% k_g_out = squeeze(CEAdata('k')); 
% k_g_t = k_g_out(2); 
% M_out = squeeze(CEAdata('mach'));
% M_t = M_out(2);
% cf_g_out = squeeze(CEAdata('cf')); 
% cf_g_t = cf_g_out(2);
% Pr_g_out = squeeze(CEAdata('prandtl')); 
% Pr_g_t = Pr_g_out(2); 
% Pr_g_tot = Pr_g_out(1);
% gamma_out = squeeze(CEAdata('gammas'));
% gamma_t = gamma_out(2);
% mu_g_out = squeeze(CEAdata('visc'));
% mu_g_t = mu_g_out(2); 
% mu_g_tot = mu_g_out(1);
% cstar_out = squeeze(CEAdata('cstar'));
% c_star = cstar_out(1);

%% Units
%ensure you are using SI units in solver
P_c_SI = P_c * u.PSI2PA;
%D_hydraulic_SI = D_hydraulic .* u.IN2M;
%A_channel_SI = A_channel .* u.IN2M^2; 
%A_channel_f_SI = A_channel_f .* u.IN2M^2;
%r_wg_SI = r_wg .* u.IN2M; 
%w_land_SI = w_land .* u.IN2M;
%w_c_SI = w_c .* u.IN2M;
%r_backside_SI = r_backside .* u.IN2M;
%h_c_SI = h_c .* u.IN2M;
%kappa = kappa .* 1 / u.IN2M; % 1 / [m]
%rho_m_SI = rho_m .* u.IN2M; 
%rho_wl_SI = rho_wl .* u.IN2M;
%rho_backside = rho_backside .* u.IN2M;
%dxds = dxds_midpoints;
%dxds_f = dxds_bounds;
%t_w_SI = t_w .* u.IN2M;
%ds_SI = ds .* u.IN2M; % Arc length for each segment [m]
%ds_ch_SI = ds_ch .* u.IN2M;



%% Initialization
steps = length(x_midpoints);
P_kpa = nan(1,steps);
% Initial Temperature guesses
Twg_c = 700;
Twl_c = 300;
Tm_c = (Twg_c + Twl_c) ./ 2;  
Twg_l = 700;
Tfinbase = 300;
Tm_l = (Twg_l + Tfinbase) ./ 2; 
ones_vector = ones(1,steps);
%T0_matrix = [ones_vector.*Twg_c; ones_vector.*Tm_c; ones_vector.*Twl_c; 
%     ones_vector.*Twg_l; ones_vector.* Tm_l; ones_vector.* Tfinbase]; %Initial Temperature Guesses
Twg_c = ones_vector.*Twg_c;
Twl_c = ones_vector.*Twl_c;
Twg_l = ones_vector.*Twg_l;
Tfinbase = ones_vector.*Tfinbase;
Tm_c = ones_vector.*Tm_c;
Tm_l = ones_vector.*Tm_l;
T0_matrix = [Tm_c; Tm_l];
T = T0_matrix(:);



% Initial Coolant States
steps_f = length(x_bounds);
ones_vector_coolant = ones(1,steps_f);
P_l = ones_vector_coolant .* (inlet_pressure * u.PSI2PA); % [Pa]
hmass_inlet = double(py.CoolProp.CoolProp.PropsSI('Hmass','T', inlet_temp, 'P', P_l(1), coolant)); % viscosity of bulk coolant [Pa-s]
hmass_l = ones_vector_coolant .* hmass_inlet; % enthalpy if bulk ccoolant [J/kg]

% Preallocations
rho_l = nan(1,steps);
T_l = nan(1,steps);
rho_l_f = nan(1,steps_f);
cp_l = nan(1,steps);
mu_l = nan(1,steps);
k_l = nan(1,steps);
%vel_coolant = nan(1,steps); % m/s
vel_coolant_f = nan(1,steps_f);
Re_l = nan(1,steps);
Re_l_f = nan(1,steps_f);
Pr_l = nan(1,steps); 
h_l = nan(1,steps);
h_g_c = nan(1,steps);
h_g_l = nan(1,steps);
sigma_l = nan(1,steps);
sigma_c = nan(1,steps);
G1 = nan(1, steps);
G2 = nan(1, steps);
G3 = nan(1, steps);
G4 = nan(1, steps);
G5 = nan(1, steps);
G6 = nan(1, steps);
G7 = nan(1, steps);
G8 = nan(1, steps);
G9 = nan(1, steps);
G10 = nan(1,steps);
G11 = nan(1,steps);
G12 = nan(1,steps);
G13 = nan(1,steps);
Ggc = nan(1,steps);
Ggl = nan(1,steps);
Glc = nan(1,steps);
Gll = nan(1,steps);
GAC = nan(1,steps-1);
GAL = nan(1,steps-1);
k3 = nan(1, steps);
k4 = nan(1, steps);
k5 = nan(1, steps);
k6 = nan(1, steps);
k7 = nan(1, steps);
k9 = nan(1, steps);
k10 = nan(1,steps-1);
k11 = nan(1,steps-1);
ed = nan(1, steps);
f =  nan(1, steps);
Nu_l = nan(1, steps);
m_fin = nan(1, steps);
eta_fin = nan(1,steps);
lambda = nan(1,steps);
Gfin = nan(1,steps);
vel_coolant = nan(1,steps);
phase = nan(1,steps);
quality_l = nan(1,steps);
phase_f = nan(1,steps_f);
quality_l_f = nan(1,steps_f);
deltaP_major = nan(1,steps);
deltaP_accel = nan(1,steps);
deltaP_minor = zeros(1,steps);
deltaP_major_cum = nan(1,steps_f);
q5_integrated = nan(1,steps);
G5_old = nan(1,steps);
QAC = nan(1,steps);
QAL = nan(1,steps);
QAC_base = nan(1,steps);
QAL_base = nan(1,steps);
QAC_corr = nan(1,steps);
QAL_corr = nan(1,steps);
QAC_recon = nan(1,steps);
QAL_recon = nan(1,steps);
qax_corr_c = nan(1,steps);
qax_corr_l = nan(1,steps);
q5_recon = nan(1,steps);
q5_corr = nan(1,steps);
q5_base = nan(1,steps);



%% Inlet Coolant calcs
% velocity 
rho_l_f(1) = double(py.CoolProp.CoolProp.PropsSI('D','Hmass',hmass_l(1), 'P', P_l(1), coolant));
vel_coolant_f(1) = mdot_channel / (rho_l_f(1) * A_channel_f(1));
phase_f(1) = py.CoolProp.CoolProp.PropsSI('PHASE','Hmass', hmass_l(1),'P', P_l(1), coolant); % phase of the coolant 
quality_l_f(1) = double(py.CoolProp.CoolProp.PropsSI('Q','Hmass',hmass_l(1),'P',P_l(1),coolant));

%% Picard Solution 
r = Pr_g .^ (1 / 3); % recovery factor (turbulent flow)
T_r = T_g .* (1 + (gamma - 1) ./ 2 .* r .* M .^ 2); % recovery temperature [K] - corrects for viscous dissipation in boundary layer (Heister EQ 6.15).

% Numerical Solution
converged = 0;
temp_tol =.1;
pressure_tol = .1;
enthalpy_tol = .1;
debug = 0;


while ~converged
    % Store Previous guesses
    T_old = T;
    % T_old_shaped = reshape(T_old,6,steps);
    % Twg_c = T_old_shaped(1,:);
    % Tm_c = T_old_shaped(2,:);
    % Twl_c = T_old_shaped(3,:);  
    % Twg_l = T_old_shaped(4,:);
    % Tm_l = T_old_shaped(5,:);
    % Tfinbase = T_old_shaped(6,:); 
    P_l_old = P_l;
    hmass_l_old = hmass_l;

    % Initialize arrays for this interation
    if debug 

    end

    for i = 1:steps

        % Conductance 1, Hot gas convection to channel hot wall
        sigma_c(i) = (.5 * Twg_c(i) / T_g_tot * (1 + (gamma(i) - 1) / 2 * M(i) ^ 2) + .5) ^ -.68 * (1 + (gamma(i) - 1) / 2 * M(i) ^ 2) ^ -.12; % film coefficient correction factor [N/A] (Huzel & Huang 86).
        h_g_c(i) = heatflux_factor * (0.026 / D_t ^ 0.2) * (mu_g_tot ^ 0.2 * cp_g_tot / Pr_g_tot ^ 0.6) * (P_c_SI / cstar) ^ 0.8 * (D_t / R_of_curve) ^ 0.1 * (1 / A_ratio(i)) ^ .9 * sigma_c(i); % gas film coefficient [W/m^2-K] - bartz equation (Huzel & Huang 86).
        %R1 = 1 / ( h_g_c * w_c_t);
        %R1 = 1 / (h_g_c * theta_channel * r_wg_t);
        G1(i) = h_g_c(i) * theta_channel(i) * r_wg(i); 

        % Conductance 2: Hot gas convection to land hot wall
        sigma_l(i) = (.5 * Twg_l(i) / T_g_tot * (1 + (gamma(i) - 1) / 2 * M(i) ^ 2) + .5) ^ -.68 * (1 + (gamma(i) - 1) / 2 * M(i) ^ 2) ^ -.12; % film coefficient correction factor [N/A] (Huzel & Huang 86).
        h_g_l(i) = heatflux_factor * (0.026 / D_t ^ 0.2) * (mu_g_tot ^ 0.2 * cp_g_tot / Pr_g_tot ^ 0.6) * (P_c_SI / cstar) ^ 0.8 * (D_t / R_of_curve) ^ 0.1 * (1 / A_ratio(i)) ^ .9 * sigma_l(i); % gas film coefficient [W/m^2-K] - bartz equation (Huzel & Huang 86). 
        %R2 = 1 / ( h_g_l * w_l_t);
        %R2 = 1 / (h_g_l * theta_land * r_wg_t);
        G2(i) = h_g_l(i) * theta_land(i) * r_wg(i);

        % Conductance 3: Conduction to middle of channel thin wall
        %k3(i) = interp1(k_w(:,1), k_w(:,2), (Twg_c(i)+Tm_c(i))/2, 'linear', 'extrap');
        k3(i) = k_w((Twg_c(i)+Tm_c(i))/2);
        %R3 = (t_w_t/2) / (k3 * w_c_t);
        %R3 = log(r_m_t / r_wg_t) / (k3 * theta_channel);
        %G3(i) = (k3(i) * theta_channel(i)) / log(r_m_SI(i) / r_wg_SI(i));
        G3(i) = (k3(i) * theta_channel(i) * (dxds(i) - kappa(i)*r_wg(i))) / log((r_wg(i)+dxds(i)*t_w(i)/2) / (r_wg(i) * (1 + kappa(i)*t_w(i)/2)));
 
        % Conductance 4: Conduction to middle of land thin wall
        %k4(i) = interp1(k_w(:,1), k_w(:,2), (Twg_l(i)+Tm_l(i))/2, 'linear', 'extrap');
        k4(i) = k_w((Twg_l(i)+Tm_l(i))/2);
        %R4 = (t_w_t/2) / (k4 * w_l_t);
        %R4 = log(r_m_t / r_wg_t) / (k4 * theta_land);
        %G4(i) = (k4(i) * theta_land(i)) / log(r_m_SI(i) / r_wg_SI(i));
        G4(i) = (k4(i) * theta_land(i) * (dxds(i) - kappa(i)*r_wg(i))) / log((r_wg(i)+dxds(i)*t_w(i)/2) / (r_wg(i) * (1 + kappa(i)*t_w(i)/2)));

        % Conductance 5: Lateral Conductance between Channel and Land
        % midpoints
        %k5(i) = interp1(k_w(:,1), k_w(:,2), (Tm_c(i)+Tm_l(i))/2, 'linear', 'extrap');
        k5(i) = k_w((Tm_c(i)+Tm_l(i))/2);
        %R4 = ((w_c_t + w_l_t)/2) / (k5 * t_w_t);
        %R5 = (theta_p/2) / (k5 * log(r_wl/r_wg_t));
        if latconductance 
            %G5(i) = (k5(i) * log(r_wl_SI(i)/r_wg_SI(i))) / (theta_p/2);
            G5(i) = (2*k5(i)/(theta_p/2)) * (kappa(i)*t_w(i)/dxds(i) + (dxds(i) - kappa(i)*r_wg(i))/(dxds(i)^2) * log((r_wg(i)+dxds(i)*t_w(i))/r_wg(i)));
            if temperaturefieldlateral
                Tc = wallTempProfile(zeta,Twg_c(i),Tm_c(i),Twl_c(i));
                Tl = wallTempProfile(zeta,Twg_l(i),Tm_l(i),Tfinbase(i));

                DeltaK = K_ad(Tc) - K_ad(Tl);

                q5_recon(i) = trapz(zeta,Clat(i,:).*DeltaK);

                q5_base(i) = G5(i)*(Tm_c(i)-Tm_l(i));
                
                q5_corr(i) = q5_recon(i)-q5_base(i);

                % q5_integrated(i) = trapz(zeta,Clat(i,:).*DeltaK);
                % 
                % if abs(Tm_c(i) - Tm_l) < 1e-8
                %     G5(i) = 0;
                % else
                %     G5(i) = q5_integrated(i) / (Tm_c(i) - Tm_l(i));
                % end
                %G5_old(i) = (2*k5(i)/(theta_p/2)) * (kappa(i)*t_w(i)/dxds(i) + (dxds(i) - kappa(i)*r_wg(i))/(dxds(i)^2) * log((r_wg(i)+dxds(i)*t_w(i))/r_wg(i)));
            else 
                q5_corr(i) = 0;
            end
        else 
            G5(i) = 0;
            q5_corr(i) = 0;
        end

        % Conductance 6: Conduction to liquid side wall from mid point of
        % channel thin wall
        %k6(i) = interp1(k_w(:,1), k_w(:,2), (Tm_c(i)+Twl_c(i))/2, 'linear', 'extrap');
        k6(i) = k_w((Tm_c(i)+Twl_c(i))/2);
        %R6 = (t_w_t/2) / (k6 * w_c_t);
        %R6 = log(r_wl / r_m_t) / (k6 * theta_channel);
        %G6(i) = (k6(i) * theta_channel(i)) / log(r_wl_SI(i) / r_m_SI(i));
        G6(i) = (k6(i) * theta_channel(i) * (dxds(i) - kappa(i)*r_wg(i))) / log(((r_wg(i) + dxds(i)*t_w(i))*(1 + kappa(i)*t_w(i)/2)) / ((r_wg(i) + dxds(i)*t_w(i)/2)*(1 + kappa(i)*t_w(i))));

        % Conductance 7:
        %k7(i) = interp1(k_w(:,1), k_w(:,2), (Tm_l(i)+Tfinbase(i))/2, 'linear', 'extrap'); 
        k7(i) = k_w((Tm_l(i)+Tfinbase(i))/2);
        %R7 = (t_w_t/2) / (k7 * w_l_t);
        %R7 = log(r_wl_SI(i) / r_m_SI(i)) / (k7(i) * theta_land(i));
        %G7(i) = (k7(i) * theta_land(i)) / log(r_wl_SI(i) / r_m_SI(i));
        G7(i) = (k7(i) * theta_land(i) * (dxds(i) - kappa(i)*r_wg(i))) / log(((r_wg(i) + dxds(i)*t_w(i))*(1 + kappa(i)*t_w(i)/2)) / ((r_wg(i) + dxds(i)*t_w(i)/2)*(1 + kappa(i)*t_w(i))));

        % Conductance 8: Cold Wall Convection 
        
        if coolant == "isopropyl alcohol"
            % P_kPa = P_l(i) / 1000;
            % T_K = T_l(i);
            % 
            % rho_l(i) = Data.F_rho(T_K, P_kPa);
            % cp_l(i)  = Data.F_cp(T_K, P_kPa) * 1000; 
            % mu_data = Data.F_mu_data(T_K);
            % k_data  = Data.F_k_data(T_K);
            % 
            % if (T_K - 273.15) < 36.5
            %     mu_fit = 4.5054 * exp(-0.031 * (T_K - 273.15)) / 1000;
            % else
            %     mu_fit = 3.3724 * exp(-0.023 * (T_K - 273.15)) / 1000;
            % end
            % 
            % mu_avg = (mu_fit + mu_data) / 2;
            % mu_l(i) = mu_avg * exp(5.4574e-9 * (P_l(i) - 101325));
            % k_l(i)  = k_data * (1 + 2.9531e-9 * (P_l(i) - 101325));
        else
            % construct midpoint state 
            hmass_l_mp = (hmass_l(i) + hmass_l(i+1)) / 2;
            P_l_mp = (P_l(i) + P_l(i+1)) / 2;

            mu_l(i) = double(py.CoolProp.CoolProp.PropsSI('V','Hmass',hmass_l_mp, 'P', P_l_mp, coolant)); % viscosity of bulk coolant [Pa-s]
            cp_l(i) = double(py.CoolProp.CoolProp.PropsSI('C' , 'Hmass', hmass_l_mp, 'P', P_l_mp, coolant)); % specific heat of coolant [J/kg-k] 
            k_l(i) = double(py.CoolProp.CoolProp.PropsSI('L', 'Hmass', hmass_l_mp, 'P', P_l_mp, coolant)); % thermal conductivity of coolant [W/m-K]
            rho_l(i) = double(py.CoolProp.CoolProp.PropsSI('D','Hmass', hmass_l_mp,'P', P_l_mp, coolant)); % density of the coolant [kg/m^3]
            T_l(i) = double(py.CoolProp.CoolProp.PropsSI('T','Hmass', hmass_l_mp,'P', P_l_mp, coolant)); % Temperature of the coolant [K]
            phase(i) = py.CoolProp.CoolProp.PropsSI('PHASE','Hmass', hmass_l_mp,'P', P_l_mp, coolant); % phase of the coolant 
            quality_l(i) = double(py.CoolProp.CoolProp.PropsSI('Q','Hmass',hmass_l_mp,'P',P_l_mp,coolant));
        end
    
        vel_coolant(i) = mdot_channel / (rho_l(i) * A_channel(i));
        Re_l(i) = (rho_l(i) * vel_coolant(i) * D_hydraulic(i)) / mu_l(i); % reynolds number for channel flow [N/A] (Huzel and Huang , pg 90) 

        % Use moody diagram to find coefficient of friction
        ed(i) = e(2)/(D_hydraulic(i)*1000); % relative roughness
        f(i) = moody(ed(i), Re_l(i)); % friction factor

        if phase(i) == 0
            Pr_l(i) = (cp_l(i) * mu_l(i)) / k_l(i); % prantl number [N/A] (Huzel and Huang, pg 90) 
            Nu_l(i) = (f(i) / 8) * (Re_l(i) - 1000) * Pr_l(i) / (1 + 12.7 * (f(i) / 8) ^ 0.5 * (Pr_l(i) ^ (2/3) - 1)); % Gnielinksy correlation nusselt number [N/A] - 0.5 < Pr < 2000, 3000 < Re < 5e6
            h_l(i) = coolant_heatflux_factor *  (Nu_l(i) * k_l(i)) / D_hydraulic(i); % liquid film coefficient [W/m^2-K] (Heister EQ 6.19)
        else 
            % handle two phase flow using Gungor and Winteron correlation
            % as described in "Terrestrail cryogenic two-Phase flow and
            % heat transfer nucleate flowb boiling heat transfer during
            % cryogenic chilldown" NTRS
            % TAKE THESE RESUTS WITH HUGE GRAIN OF SALT
            % NOT IMPLEMENTED YET JUST TAKING VALUES FROM LAST STEP SO
            % MARCH CONTINUES, ANYTHING OUTSIDE OF LIQUID PHASE IS GARBAGE
            %Re_l_sat = (rho_l(i) * vel_coolant(i) * D_hydraulic_SI(i) * (1 - quality(i))) / mu_l(i);
            Pr_l(i) = Pr_l(i-1);
            Nu_l(i) = Nu_l(i-1);
            h_l(i) = h_l(i-1);
        end

        %Nu_l(i) = 0.023 * (Re_l(i) ^ .8) * (Pr_l(i) ^ .4) * (T_wl(i) / T_l(i)) ^ -.3; % nusselt number [N/A] - applicable for Re > 10,000, .7 < Pr < 160 (Heister EQ 6.19). ****
        %Nu_l(i) = (f(i) / 8) * (Re_l(i) - 1000) * Pr_l(i) / (1 + 12.7 * (f(i) / 8) ^ 0.5 * (Pr_l(i) ^ (2/3) - 1)); % Gnielinksy correlation nusselt number [N/A] - 0.5 < Pr < 2000, 3000 < Re < 5e6
        %h_l(i) = coolant_heatflux_factor *  (Nu_l(i) * k_l(i)) / D_hydraulic_SI(i); % liquid film coefficient [W/m^2-K] (Heister EQ 6.19)

        %R8 = 1 /(h_l_t * theta_channel * r_wl);
        %G8(i) = h_l(i) * theta_channel(i) * r_wl_SI(i);
        G8(i) = h_l(i) * theta_channel(i) * (r_wg(i) + dxds(i)*t_w(i)) * (1 + kappa(i)*t_w(i));

        % Conductance 9: Fin to cooling channel Bessel function formulation
        if ~isnan(eta_fin(i))
            T_finavg = T_l(i) + eta_fin(i)*(Tfinbase(i) - T_l(i));
            %k9(i) = interp1(k_w(:,1), k_w(:,2), T_finavg, 'linear', 'extrap');
            k9(i) = k_w(T_finavg);
        else
            %k9(i) = interp1(k_w(:,1), k_w(:,2), Tfinbase(i), 'linear', 'extrap');
            k9(i) = k_w(Tfinbase(i));
        end
        %m_fin(i) = sqrt(2*h_l(i) / (k9(i)*w_land_SI(i))); % 
        %eta_fin(i) = fin_efficiency_factor * tanh(m_fin(i) * (r_backside_SI(i) - r_wl_SI(i))) / (m_fin(i)*(r_backside(i) - r_wl_SI(i)));  % Fin efficiency
        % G9(i) = 2 * eta_fin(i) * h_l(i) * (r_backside(i) - r_wl_SI(i));
        lambda(i) = sqrt(2*h_l(i) / (k9(i)*theta_land(i)*dxds(i)^2)); % Variable-area annular-sector fin parameter
        u_base = t_w(i);
        u_tip = t_w(i) + h_c(i);
        y_fin = @(u) r_wg(i) + dxds(i)*u;
        y_base = y_fin(u_base);
        y_tip = y_fin(u_tip);
        z_fin = @(u) 2*lambda(i)*sqrt(y_fin(u));
        z_base =z_fin(u_base); % Dimensionless Bessel coordinate fin base
        z_tip = z_fin(u_tip); % Dimensionless Bessel coordinate fin tip
        beta = besseli(1,z_tip) / besselk(1,z_tip);
        fin_shape = @(u) besseli(0,z_fin(u)) + beta*besselk(0,z_fin(u));
        tempIntegral = integral(@(u) fin_shape(u), u_base,u_tip);
        tempMoment = integral(@(u) u.*fin_shape(u),u_base,u_tip);
        u_eff = tempMoment/tempIntegral;
        J_eff = 1 + kappa(i)*u_eff;
        % Modified Besel functions
        I0_base = besseli(0,z_base);
        I1_base = besseli(1,z_base);
        K0_base = besselk(0,z_base);
        K1_base = besselk(1,z_base);
        I1_tip = besseli(1,z_tip);
        K1_tip = besselk(1,z_tip);
        numerator_bessel = I1_tip*K1_base - K1_tip*I1_base;
        denominator_bessel = I0_base*K1_tip + K0_base*I1_tip;
        bessel_ratio = numerator_bessel / denominator_bessel; % result of bessel solution to ODE, needs to be rederived with midarchlength
        Gfin(i) = k9(i)*theta_land(i)*dxds(i)*lambda(i)*sqrt(y_base)*bessel_ratio;
        G9(i) = Gfin(i) * fin_efficiency_factor * J_eff;
        eta_fin(i) = Gfin(i) / (2 * h_l(i) * h_c(i));
        %G9(i) = Gfin(i) * fin_efficiency_factor;    

        %Conductance 10, 11 & 12, 13: Axial Conductance Channel & land 
        if axialconductance
            if i < steps
                k10(i) = k_w((Tm_c(i)+Tm_c(i+1))/2);
                k11(i) = k_w((Tm_l(i)+Tm_l(i+1))/2);
    
                % if abs(kappa_f(i+1)*t_w_f(i+1)) < 1e-9
                %     Is = r_wg_f(i+1)*t_w_f(i+1) + (dxds_f(i+1)/2)*(t_w_f(i+1)^2);
                % else 
                %     Is = dxds_f(i+1)*t_w_f(i+1)/kappa_f(i+1) + ((r_wg_f(i+1)*kappa_f(i+1)-dxds_f(i+1))/((kappa_f(i+1)^2)))*log(1+kappa_f(i+1)*t_w_f(i+1));                  
                % end
                % GAC(i) = Is*k10(i)*theta_channel_f(i+1)/ds_midpoints(i); %% Axial conductance channel**W/K**
                % GAL(i) = Is*k11(i)*theta_land_f(i+1)/ds_midpoints(i);  %% Axial conductance land**W/K**
                if temperaturefieldaxial
                    Tc_i = wallTempProfile(zeta,Twg_c(i),Tm_c(i),Twl_c(i));
                    Tc_ip1 = wallTempProfile(zeta,Twg_c(i+1),Tm_c(i+1),Twl_c(i+1));
                    Tl_i = wallTempProfile(zeta,Twg_l(i),Tm_l(i),Tfinbase(i));
                    Tl_ip1 = wallTempProfile(zeta,Twg_l(i+1),Tm_l(i+1),Tfinbase(i+1));

                    DeltaKc = K_ad(Tc_i) - K_ad(Tc_ip1);
                    DeltaKl = K_ad(Tl_i) - K_ad(Tl_ip1);

                    % QAC(i) = trapz(zeta,HAC(i,:).*DeltaKc);
                    % QAL(i) = trapz(zeta,HAL(i,:).*DeltaKl);
                    % 
                    % DeltaTm_c = Tm_c(i) - Tm_c(i+1);
                    % DeltaTm_l = Tm_l(i) - Tm_l(i+1);
                    % if QAC(i)*DeltaTm_c < 0
                    %     fprintf('Channel axial sign mismatch at i = %d\n',i);
                    % end
                    % if QAL(i)*DeltaTm_l < 0
                    %     fprintf('Land axial sign mismatch at i = %d\n',i);
                    % end
                    % if abs(Tm_c(i) - Tm_c(i+1)) < 1e-8
                    %     GAC(i) = k10(i)/RgAC(i); %% Axial conductance channel**W/K**
                    % else
                    %     GAC(i) = QAC(i) / (Tm_c(i) - Tm_c(i+1));
                    % end
                    % if abs(Tm_l(i) - Tm_l(i+1)) < 1e-8
                    %     GAL(i) = k11(i)/RgAL(i);  %% Axial conductance land**W/K**
                    % else
                    %     GAL(i) = QAL(i) / (Tm_l(i) - Tm_l(i+1));
                    % end
                    % Positive baseline conductances
                    GAC(i) = k10(i)/RgAC(i);
                    GAL(i) = k11(i)/RgAL(i);

                    % Midpoint-model axial heat flow
                    QAC_base(i) = GAC(i) * (Tm_c(i)-Tm_c(i+1));
                    QAL_base(i) = GAL(i) * (Tm_l(i)-Tm_l(i+1));

                    % Step-3 reconstructed heat flow
                    QAC_recon(i) = trapz(zeta,HAC(i,:).*DeltaKc);
                    QAL_recon(i) = trapz(zeta,HAL(i,:).*DeltaKl);

                    % Explicit correction
                    QAC_corr(i) = QAC_recon(i) - QAC_base(i);
                    QAL_corr(i) = QAL_recon(i) - QAL_base(i);

                   
                else
                    GAC(i) = k10(i)/RgAC(i); %% Axial conductance channel**W/K**
                    GAL(i) = k11(i)/RgAL(i);  %% Axial conductance land**W/K**
                    QAC_corr(i) = 0;
                    QAL_corr(i) = 0;
                end

            end
            qax_corr_c(i) = 0;
            qax_corr_l(i) = 0;
            if i > 1
                qax_corr_c(i) = qax_corr_c(i) + QAC_corr(i-1)/ds(i);
                qax_corr_l(i) = qax_corr_l(i) + QAL_corr(i-1)/ds(i);
            end
            if i < steps
                qax_corr_c(i) = qax_corr_c(i) - QAC_corr(i)/ds(i);
                qax_corr_l(i) = qax_corr_l(i) - QAL_corr(i)/ds(i);
            end
            
            if i == 1 
                G10(i) = 0;
                G12(i) = 0;
            else 
                G10(i) = GAC(i-1)/ds(i); % [W/mK] Axial conductance from previous cell normalized by hotgas wall arclength, channel
                G12(i) = GAL(i-1)/ds(i); % [W/mK] Axial conductacnce from previous cell normalized by hotgas wall arclength, land
            end
            if i == steps 
                G11(i) = 0;
                G13(i) = 0;
            else
                G11(i) = GAC(i)/ds(i); % [W/mK] Axial conductance from next cell normalized by hotgas wall arclength, channel
                G13(i) = GAL(i)/ds(i); % [W/mK] Axial conductance from next cell normalized by hotgas wall arclength, land
            end
        else 
            G10(i) = 0;
            G11(i) = 0;
            G12(i) = 0;
            G13(i) = 0;
            qax_corr_c(i) = 0;
            qax_corr_l(i) = 0;
        end

    end

    % Calculate Condensed Conductances
    Ggc = (1./G1 + 1./G3).^(-1); % Recovery temp to mid channel temp
    Ggl = (1./G2 + 1./G4).^(-1); % Recovery temp to mid land temp
    Glc = (1./G6 + 1./G8).^(-1); % Mid channel temp to liquid temp
    Gll = (1./G7 + 1./G9).^(-1); % Mid lamd temp to liquid temp

    % Construct the conductance matrixn
    %numNodes = 6;
    numNodes = 2;
    totalNodes = numNodes * steps;
    count = 0;

    % Maximum 36 entries per 6x6 block
    % I = zeros(36*steps,1);
    % J = zeros(36*steps,1);
    % V = zeros(36*steps,1);

    % Maximum 4 local matrix entries + 2 previous neighbor entries + 2 next
    % neightbor entries: 8 per step
    I = zeros(8*steps,1);
    J = zeros(8*steps,1);
    V = zeros(8*steps,1);

    b = zeros(totalNodes,1);

    for j = 1:steps

        % Conductances at axial station i
        % g1 = G1(j);
        % g2 = G2(j);
        % g3 = G3(j);
        % g4 = G4(j);
        % g5 = G5(j);
        % g6 = G6(j);
        % g7 = G7(j);
        % g8 = G8(j);
        % g9 = G9(j);
        % 
        % Ai = [
        %         (g1+g3),        -g3,        0,        0,           0,        0;
        %             -g3, (g3+g5+g6),      -g6,        0,         -g5,        0;
        %              0,         -g6,  (g6+g8),        0,           0,        0;
        %              0,           0,        0,  (g2+g4),         -g4,        0;
        %              0,         -g5,        0,      -g4,  (g4+g5+g7),      -g7;
        %              0,           0,        0,        0,         -g7,  (g7+g9);
        %     ];
        % bi = [
        %      g1*T_r(j); 
        %              0; 
        %      g8*T_l(j); 
        %      g2*T_r(j); 
        %              0; 
        %      g9*T_l(j)
        %      ];

        % 2D Conductances at axial station i
        ggc = Ggc(j);
        ggl = Ggl(j);
        glc = Glc(j);
        gll = Gll(j);
        g5 = G5(j);

        % Axial conductances already normalized by ds(j)
        g10 = G10(j);   % channel, previous node
        g11 = G11(j);   % channel, next node
        g12 = G12(j);   % land, previous node
        g13 = G13(j);   % land, next node

        % Ai = [
        %        (ggc+glc+g5),           -g5;
        %                 -g5,  (ggl+gll+g5);     
        %      ];
        % bi = [
        %        ggc*T_r(j) + glc*T_l(j);
        %        ggl*T_r(j) + gll*T_l(j);
        %      ];
        Ai = [
               (ggc+glc+g5+g10+g11),                   -g5;
                                -g5,  (ggl+gll+g5+g12+g13);     
             ];
        bi = [
               ggc*T_r(j) + glc*T_l(j) + qax_corr_c(j) - q5_corr(j);
               ggl*T_r(j) + gll*T_l(j) + qax_corr_l(j) + q5_corr(j);
             ];
        

        % Global indices belonging to axial station i 
        %idx = (6*(j-1)+1):(2*j);
        %idx = (2*(j-1)+1):(2*j);
        idx_c = 2*j - 1;
        idx_l = 2*j;
        idx = [idx_c, idx_l];

        [ii,jj,vv] = find(Ai); 

        n = length(vv);

        I(count+1:count+n) = idx(ii);
        J(count+1:count+n) = idx(jj);
        V(count+1:count+n) = vv;

        count = count + n;

        b(idx) = bi;

        % Axial connection to previous station

        if j > 1
    
            % Channel: Tm_c(j-1)
            count = count + 1;
            I(count) = idx_c;
            J(count) = idx_c - 2;
            V(count) = -g10;
    
            % Land: Tm_l(j-1)
            count = count + 1;
            I(count) = idx_l;
            J(count) = idx_l - 2;
            V(count) = -g12;
    
        end
    
    
        % Axial connection to next station
    
        if j < steps
    
            % Channel: Tm_c(j+1)
            count = count + 1;
            I(count) = idx_c;
            J(count) = idx_c + 2;
            V(count) = -g11;
    
            % Land: Tm_l(j+1)
            count = count + 1;
            I(count) = idx_l;
            J(count) = idx_l + 2;
            V(count) = -g13;
    
        end

    end

    % Trim unused preallocation 
    I = I(1:count);
    J = J(1:count);
    V = V(1:count);

    % Consturct sparse matrix
    A = sparse(I,J,V,totalNodes,totalNodes);
    T = A \ b;

    % T_shaped = reshape(T,6,steps);
    % Twg_c = T_shaped(1,:);
    % Tm_c = T_shaped(2,:);
    % Twl_c = T_shaped(3,:);  
    % Twg_l = T_shaped(4,:);
    % Tm_l = T_shaped(5,:);
    % Tfinbase = T_shaped(6,:); 
    T_shaped = reshape(T,2,steps);
    Tm_c = T_shaped(1,:);
    Tm_l = T_shaped(2,:);

    qgc = Ggc .* (T_r - Tm_c);
    qgl = Ggl .* (T_r - Tm_l);
    qlc = Glc .* (Tm_c - T_l);
    qll = Gll .* (Tm_l - T_l);

    Twg_c = T_r - qgc./G1;
    Twg_l = T_r - qgl./G2;
    Twl_c = T_l + qlc./G8;
    Tfinbase = T_l + qll./G9;

    Qgc = qgc .* ds;
    Qgl = qgl .* ds;
    Qlc = qlc .* ds;
    Qll = qll .* ds;

    % q1 = G1 .* (T_r - Twg_c);
    % q2 = G2 .* (T_r - Twg_l);
    % q3 = G3 .* (Twg_c - Tm_c);
    % q4 = G4 .* (Twg_l - Tm_l);
    % q5 = G5 .* (Tm_c - Tm_l);
    % q6 = G6 .* (Tm_c - Twl_c);
    % q7 = G7 .* (Tm_l - Tfinbase);
    % q8 = G8 .* (Twl_c - T_l);
    % q9 = G9 .* (Tfinbase - T_l);
    % 
    % Q1 = q1 .* ds_SI;
    % Q2 = q2 .* ds_SI;
    % Q3 = q3 .* ds_SI;
    % Q4 = q4 .* ds_SI;
    % Q5 = q5 .* ds_SI;
    % Q6 = q6 .* ds_SI;
    % Q7 = q7 .* ds_SI;
    % Q8 = q8 .* ds_SI;
    % Q9 = q9 .* ds_SI;

    % Coolant March 
    for j = 2:steps_f
       cell = j - 1; % Index for the midpoint

       % velocity 
       
       rho_l_f(j) = double(py.CoolProp.CoolProp.PropsSI('D','Hmass',hmass_l(j), 'P', P_l(j), coolant));
       vel_coolant_f(j) = mdot_channel / (rho_l_f(j) * A_channel_f(j));
       deltaKE = (vel_coolant_f(j)^2 - vel_coolant_f(j-1)^2) / 2;
       deltaKE = 0;


       % enthalpy and recalculate vel_coolant
       hmass_l(j) = hmass_l(j-1) + (Qll(j-1) + Qlc(j-1))/mdot_channel - deltaKE;
       rho_l_f(j) = double(py.CoolProp.CoolProp.PropsSI('D','Hmass',hmass_l(j), 'P', P_l(j), coolant));
       vel_coolant_f(j) = mdot_channel / (rho_l_f(j) * A_channel_f(j));


       % Pressure loss
       % Use moody diagram to find coefficient of friction
       %if x_midpoints(cell) < geom.bounds.chamber(2)
       %     ed = e(2)/(D_hydraulic_SI(cell)*1000); % 90 degrees
       %else
       %     ed = e(1)/(D_hydraulic_SI(cell)*1000); % 45 degrees
       %end
       %ed = e/(hydraulic_D_x(i)*1000); % relative roughness
       %f = moody(ed, Re_l(cell)); % friction factor

       deltaP_major(j-1) = (f(cell) * ds_ch(cell) * rho_l(cell) * vel_coolant(cell)^2 ) / (2 * D_hydraulic(cell)); % Major loss using midpoint properties

       deltaP_accel(j-1) = rho_l(cell) * vel_coolant(cell) * (vel_coolant_f(j) - vel_coolant_f(j-1)); % pressure change due to fluid acceleration

       P_l(j) = P_l(j-1) - deltaP_major(j-1) - deltaP_accel(j-1); % Update pressure, minor losses should be applied outside of loop

       % Phase 
       phase_f(j) = py.CoolProp.CoolProp.PropsSI('PHASE','Hmass', hmass_l(j),'P', P_l(j), coolant); % phase of the coolant 
       quality_l_f(j) = double(py.CoolProp.CoolProp.PropsSI('Q','Hmass',hmass_l(j),'P',P_l(j),coolant));
        
    end

    % Check Convergence 

    temp_error = max(abs(T-T_old))
    pressure_error = max(abs(P_l - P_l_old))
    enthalpy_error = max(abs(hmass_l - hmass_l_old))

    % max_deltaKE = max(abs(deltaKE_vec));
    % max_deltaH_heat = max(abs(Qcell_vec ./ mdot_channel));
    % KE_ratio = max(abs(deltaKE_vec) ./ ...
    %            abs(Qcell_vec ./ mdot_channel));

    if temp_error < temp_tol && pressure_error < pressure_tol && enthalpy_error < enthalpy_tol
        converged = true;  
    end
    
end

% Calculate intermediate heat flows
q1 = G1 .* (T_r - Twg_c);
q2 = G2 .* (T_r - Twg_l);
q3 = G3 .* (Twg_c - Tm_c);
q4 = G4 .* (Twg_l - Tm_l);
q5 = G5 .* (Tm_c - Tm_l);
q6 = G6 .* (Tm_c - Twl_c);
q7 = G7 .* (Tm_l - Tfinbase);
q8 = G8 .* (Twl_c - T_l);
q9 = G9 .* (Tfinbase - T_l);
q5_old = G5_old .*(Tm_c - Tm_l);
q10 = nan(1,steps);
q11 = nan(1,steps);
q12 = nan(1,steps);
q13 = nan(1,steps);
for i = 1:steps
    if i > 1
        q10(i) = G10(i)*(Tm_c(i-1) - Tm_c(i));
        q12(i) = G12(i)*(Tm_l(i-1) - Tm_l(i));
    end
    if i < steps
        q11(i) = G11(i)*(Tm_c(i+1) - Tm_c(i));
        q13(i) = G13(i)*(Tm_l(i+1) - Tm_l(i));
    end
    if i == 1
        q10(i) = 0;
        q12(i) = 0;
    end
    if i == steps 
        q11(i) = 0;
        q13(i) = 0;
    end
end

Q1 = q1 .* ds;
Q2 = q2 .* ds;
Q3 = q3 .* ds;
Q4 = q4 .* ds;
Q5 = q5 .* ds;
Q6 = q6 .* ds;
Q7 = q7 .* ds;
Q8 = q8 .* ds;
Q9 = q9 .* ds;
Q10 = q10 .* ds;
Q11 = q11 .* ds;
Q12 = q12 .* ds;
Q13 = q13 .* ds;



% check energy conservation 
error1 = max(abs(q1-q3));
error2 = max(abs(q2-q4));
%error3 = max(abs(q3-q5-q6))
%error4 = max(abs(q4+q5-q7))
error3 = max(abs(q3-q5-q6+q10+q11));
error4 = max(abs(q4+q5-q7+q12+q13));
error5 = max(abs(q6-q8));
error6 = max(abs(q7-q9));
fprintf('max q1-q3 = %.6e W/m\n', error1)
fprintf('max q2-q4 = %.6e W/m\n', error2)
fprintf('max q6-q8 = %.6e W/m\n', error5)
fprintf('max q7-q9 = %.6e W/m\n', error6)
fprintf('max mid-channel residual = %.6e W/m\n', error3)
fprintf('max mid-land residual = %.6e W/m\n', error4)

% Show effect of temperature profile calculation
G5_change_pct = 100*(G5 - G5_old)./G5_old;
q5_change_pct = 100*(q5 - q5_old)./q5_old;

% Collect temps
% T_shaped = reshape(T,6,steps);
% Twg_c = T_shaped(1,:);
% Tm_c = T_shaped(2,:);
% Twl_c = T_shaped(3,:);  
% Twg_l = T_shaped(4,:);
% Tm_l = T_shaped(5,:);
% Tfinbase = T_shaped(6,:); 

% Calculate total pressure of coolant 
P_l_tot = P_l + .5 .* rho_l_f .* vel_coolant_f.^2;

% Calculate heat flux 
q1_flux = q1 ./ (r_wg .* theta_channel);
q2_flux = q2 ./ (r_wg .* theta_land);
q3_flux = q3 ./ (theta_channel .* (r_wg + dxds .* (t_w./4)) .*(1 + kappa.*(t_w./4) ));
q4_flux = q4 ./ (theta_land .* (r_wg + dxds .* (t_w./4)) .*(1 + kappa.*(t_w./4) ));
q5_flux = q5 ./ (2*(t_w + .5 .* kappa .* t_w.^2));
q6_flux = q6 ./ (theta_channel .* (r_wg + dxds .* (t_w.*(3/4))) .*(1 + kappa.*(t_w.*(3/4)) ));
q7_flux = q7 ./ (theta_land .* (r_wg+ dxds .* (t_w.*(3/4))) .*(1 + kappa.*(t_w.*(3/4)) ));
q8_flux = q8 ./ ( (r_wg + dxds.* t_w).*(1+kappa.*t_w).*theta_channel);
q9_flux = q9 ./ ( (2*h_c) .* (1 + kappa.*(t_w + h_c./2)));


% Calculate cumulative Major pressure loss
deltaP_major_cum(1) = 0; 
for j = 2:steps_f
    deltaP_major_cum(j) = deltaP_major_cum(j-1) + deltaP_major(j-1);
end


% organize phase change enthalpys
liquid_phasex = [];  % liquid convective heat flux
gas_phasex = [];  % gas convective heat flux
supercrit_phasex = [];
supercritgas_phasex = [];
supercritliq_phasex = [];
critpoint_phasex = [];
two_phasex = [];
hmass_twophase = [];
hmass_liq = [];
hmass_gas= [];
hmass_crit = [];
hmass_gascrit = [];
hmass_liqcrit = [];
hmass_critpoint = [];
if coolant ~= "isopropyl alcohol"
    for i = 1:steps_f
        if phase_f(i) == 0 
            liquid_phasex = [liquid_phasex, x_bounds(i)]; 
            hmass_liq = [hmass_liq, hmass_l(i)];
        elseif phase_f(i) == 1 
            supercrit_phasex = [supercrit_phasex, x_bounds(i)];
            hmass_crit = [hmass_crit, hmass_l(i)];
        elseif phase_f(i) == 2 
            supercritgas_phasex = [supercritgas_phasex, x_bounds(i)];
            hmass_gascrit = [hmass_gascrit, hmass_l(i)];
        elseif phase_f == 3
            supercritliq_phasex = [supercritliq_phasex, x_bounds(i)];
            hmass_liqcrit = [hmass_liqcrit, hmass_l(i)];
        elseif phase_f == 4
            critpoint_phasex = [critpoint_phasex, x_bounds(i)];
            hmass_critpoint = [hmass_critpoint, hmass_l(i)];
        elseif phase_f(i) == 5 
            gas_phasex = [gas_phasex, x_bounds(i)];
            hmass_gas = [hmass_gas, hmass_l(i)];
        elseif phase_f(i) == 6
            two_phasex = [two_phasex, x_bounds(i)];
            hmass_twophase = [hmass_twophase, hmass_l(i)];

        else
            warn("Unspecified Phase")
        end
    end
end

%% Structural Calculations 
% Heister Stress Check 
%deltaT1_c = Twg_c - Twl_c;
%sigma_tp = ((P_l-P_g)/2) .*((w_c_SI./t_w_SI).^2); % Tangential Pressure stress
%sigma_tp_cold = ((P_l)/2) .*((w_c_SI./t_w_SI).^2); % Tangential Pressure stress no chamber pressure







%% Throat Calc w/ Picard solution

% P_c = P_c * u.PSI2PA;
% 
% r_wg_t = D_t/2;
% r_wl = r_wg_t + t_w_t;
% r_outer = r_wl + h_c_t;
% r_midchan = 0.5* (r_wl + r_outer);
% theta_channel = w_c_t /r_midchan; % radians
% A_channel_t = theta_channel/2  * (r_outer^2 - r_wl^2);
% theta_p = 2*pi()/ num_channels; 
% theta_land = theta_p - theta_channel; 
% r_m_t = (r_wg_t + r_wl)/2;
% mdot_channel = mdot_coolant / num_channels;
% P_wet_t = (r_wl + r_outer)*theta_channel + 2*(r_outer - r_wl);
% D_hydrualic_t = 4 * A_channel_t / P_wet_t;
% 
% 
% 
% T_l_t = 290; % Coolant temp [K]
% P_l_t = 400 * u.PSI2PA; % Coolant Pressure [psi to pa]
% heatflux_factor = 1;
% fin_efficiency_factor = 1;
% coolant_heatflux_factor = 1;
% nolatconductance = 0;
% Twg_c = 700;
% Twl_c = 300;
% Tm_c = (Twg_c + Twl_c) ./ 2;  
% Twg_l = 700;
% Tfinbase = 300;
% Tm_l = (Twg_l + Tfinbase) ./ 2; 
% T = [700; 500; 300; 700; 500; 300]; %Initial Temperature Guesses
% 
% % Initial Guesses
% 
% 
% 
% 
% r = Pr_g_t ^ (1 / 3);
% T_r = T_g_t * (1 + (gamma_t - 1) / 2 * r * M_t ^ 2); % recovery temperature [K] - corrects for compressible boundry layers (Heister EQ 6.15). 
% 
% 
% % Try out the numerical solution at the throat 
% converged = 0;
% tolerance = .1;
% while ~converged
%     % Old guesses
%     T_old = T;
%     Twg_c = T_old(1);
%     Tm_c = T_old(2);
%     Twl_c = T_old(3);  
%     Twg_l = T_old(4);
%     Tm_l = T_old(5);
%     Tfinbase = T_old(6); 
% 
%     % Resistance 1, Hot gas convection to channel hot wall
%     sigma_c = (.5 * Twg_c / T_g_tot * (1 + (gamma_t - 1) / 2 * M_t ^ 2) + .5) ^ -.68 * (1 + (gamma_t - 1) / 2 * M_t ^ 2) ^ -.12; % film coefficient correction factor [N/A] (Huzel & Huang 86).
%     h_g_c = heatflux_factor * (0.026 / D_t ^ 0.2) * (mu_g_tot ^ 0.2 * cp_g_tot / Pr_g_tot ^ 0.6) * (P_c_SI / c_star) ^ 0.8 * (D_t / R_of_curve) ^ 0.1 * (1 / 1) ^ .9 * sigma_c; % gas film coefficient [W/m^2-K] - bartz equation (Huzel & Huang 86).
%     %R1 = 1 / ( h_g_c * w_c_t);
%     R1 = 1 / (h_g_c * theta_channel * r_wg_t);
% 
% 
%     % Resistance 2: Hot gas convection to land hot wall
%     sigma_l = (.5 * Twg_l / T_g_tot * (1 + (gamma_t - 1) / 2 * M_t ^ 2) + .5) ^ -.68 * (1 + (gamma_t - 1) / 2 * M_t ^ 2) ^ -.12; % film coefficient correction factor [N/A] (Huzel & Huang 86).
%     h_g_l = heatflux_factor * (0.026 / D_t ^ 0.2) * (mu_g_tot ^ 0.2 * cp_g_tot / Pr_g_tot ^ 0.6) * (P_c_SI / c_star) ^ 0.8 * (D_t / R_of_curve) ^ 0.1 * (1 / 1) ^ .9 * sigma_l; % gas film coefficient [W/m^2-K] - bartz equation (Huzel & Huang 86).
%     %R2 = 1 / ( h_g_l * w_l_t); 
%     R2 = 1 / (h_g_l * theta_land * r_wg_t);
% 
%     % Resistance 3: 
%     k3 = interp1(k_w(:,1), k_w(:,2), (Twg_c+Tm_c)/2, 'linear', 'extrap');
%     %R3 = (t_w_t/2) / (k3 * w_c_t);
%     R3 = log(r_m_t / r_wg_t) / (k3 * theta_channel);
% 
%     % Resisitance 4: 
%     k4 = interp1(k_w(:,1), k_w(:,2), (Twg_l+Tm_l)/2, 'linear', 'extrap');
%     %R4 = (t_w_t/2) / (k4 * w_l_t);
%     R4 = log(r_m_t / r_wg_t) / (k4 * theta_land);
% 
%     % Resistance 5:
%     k5 = interp1(k_w(:,1), k_w(:,2), (Tm_c+Tm_l)/2, 'linear', 'extrap');
%     %R4 = ((w_c_t + w_l_t)/2) / (k5 * t_w_t);
%     R5 = (theta_p/2) / (k5 * log(r_wl/r_wg_t));
% 
%     % Resistance 6:
%     k6 = interp1(k_w(:,1), k_w(:,2), (Tm_c+Twl_c)/2, 'linear', 'extrap');
%     %R6 = (t_w_t/2) / (k6 * w_c_t);
%     R6 = log(r_wl / r_m_t) / (k6 * theta_channel);
% 
%     % Resistance 7:
%     k7 = interp1(k_w(:,1), k_w(:,2), (Tm_l+Tfinbase)/2, 'linear', 'extrap');
%     %R7 = (t_w_t/2) / (k7 * w_l_t);
%     R7 = log(r_wl / r_m_t) / (k7 * theta_land);
% 
%     % Resistance 8: 
%     % Cold wall convection 
%     if coolant == "isopropyl alcohol"
%         P_kPa_t = P_l_t / 1000;
%         T_K_t = T_l_t;
% 
%         rho_l_t = Data.F_rho(T_K_t, P_kPa_t);
%         cp_l_t  = Data.F_cp(T_K_t, P_kPa_t) * 1000; 
%         mu_data = Data.F_mu_data(T_K_t);
%         k_data  = Data.F_k_data(T_K_t);
% 
%         if (T_K_t - 273.15) < 36.5
%             mu_fit = 4.5054 * exp(-0.031 * (T_K_t - 273.15)) / 1000;
%         else
%             mu_fit = 3.3724 * exp(-0.023 * (T_K_t - 273.15)) / 1000;
%         end
% 
%         mu_avg = (mu_fit + mu_data) / 2;
%         mu_l_t = mu_avg * exp(5.4574e-9 * (P_l(i) - 101325));
%         k_l_t  = k_data * (1 + 2.9531e-9 * (P_l(i) - 101325));
%     else
%         mu_l(i) = double(py.CoolProp.CoolProp.PropsSI('V','T', T_l_t, 'P', P_l_t, coolant)); % viscosity of bulk coolant [Pa-s]
%         cp_l(i) = double(py.CoolProp.CoolProp.PropsSI('C' , 'T', T_l_t, 'P', P_l_t, coolant)); % specific heat of coolant [J/kg-k] 
%         k_l(i) = double(py.CoolProp.CoolProp.PropsSI('L', 'T', T_l_t, 'P', P_l_t, coolant)); % thermal conductivity of coolant [W/m-K]
%         rho_l(i) = double(py.CoolProp.CoolProp.PropsSI('D','T', T_l_t,'P', P_l_t, coolant)); % density of the coolant [kg/m^3]
%     end
% 
%     vel_coolant = mdot_channel / rho_l_t / A_channel_t;
%     Re_l_t = (rho_l_t * vel_coolant * D_hydrualic_t) / mu_l_t; % reynolds number for channel flow [N/A] (Huzel and Huang , pg 90)
%     Pr_l_t = (cp_l_t * mu_l_t) / k_l_t; % prantl number [N/A] (Huzel and Huang, pg 90) 
% 
%     % Use moody diagram to find coefficient of friction
%     ed = e(2)/(D_hydrualic_t*1000); % relative roughness
%     f = moody(ed, Re_l_t); % friction factor
% 
%     %Nu_l(i) = 0.023 * (Re_l(i) ^ .8) * (Pr_l(i) ^ .4) * (T_wl(i) / T_l(i)) ^ -.3; % nusselt number [N/A] - applicable for Re > 10,000, .7 < Pr < 160 (Heister EQ 6.19). ****
%     Nu_l_t = (f / 8) * (Re_l_t - 1000) * Pr_l_t / (1 + 12.7 * (f / 8) ^ 0.5 * (Pr_l_t ^ (2/3) - 1)); % Gnielinksy correlation nusselt number [N/A] - 0.5 < Pr < 2000, 3000 < Re < 5e6
%     h_l_t = coolant_heatflux_factor *  (Nu_l_t * k_l_t) / D_hydrualic_t; % liquid film coefficient [W/m^2-K] (Heister EQ 6.19)
% 
%     R8 = 1 /(h_l_t * theta_channel * r_wl);
% 
%     % Resistance 9 
%     k9 = interp1(k_w(:,1), k_w(:,2), Tfinbase, 'linear', 'extrap');
%     r_midfin = 0.5*(r_wl + r_outer);
%     w_land_midfin = r_midfin * theta_land;
%     m_fin = sqrt(2*h_l_t / (k9*w_land_midfin)); % Constant fin width approximation use midwidth of fin
%     eta_fin = fin_efficiency_factor * tanh(m_fin * (r_outer - r_wl)) / (m_fin*(r_outer - r_wl));  % Fin efficiency
%     R9 = 1/ ( eta_fin * h_l_t * 2 * (r_outer-r_wl)); 
% 
%     % Convert resistances to Conductances
%     G1 = 1/R1; 
%     G2 = 1/R2;
%     G3 = 1/R3;
%     G4 = 1/R4; 
%     G5 = 1/R5;
%     G6 = 1/R6; 
%     G7 = 1/R7; 
%     G8 = 1/R8; 
%     G9 = 1/R9; 
%     if nolatconductance 
%         G5 = 0;
%     end
% 
%     % Construct Conductance Matrix
%     A = [-(G1+G3),          G3,        0,        0,           0,        0;
%                G3, -(G3+G5+G6),       G6,        0,          G5,        0; 
%                 0,          G6, -(G6+G8),        0,           0,        0;
%                 0,           0,        0, -(G2+G4),          G4,        0;
%                 0,          G5,        0,       G4, -(G4+G5+G7),       G7;
%                 0,           0,        0,        0,          G7, -(G7+G9) ];
%     b = [-G1*T_r; 0; -G8*T_l_t; -G2*T_r; 0; -G9*T_l_t];
% 
%     % Solve for T 
%     T = A \ b; 
% 
% 
%     % Check convergence 
%     error_max = max(abs(T-T_old))
%     if error_max < tolerance
%         converged = true;
%     end
% 
% end
% 
% h = T
% q1 = (T_r - T(1)) / R1
% q2 = (T_r - T(4)) / R2
% q3 = (T(1) - T(2)) / R3
% q4 = (T(4) - T(5)) / R4
% q5 = (T(2) - T(5)) / R5
% q6 = (T(2) - T(3)) / R6
% q7 = (T(5) - T(6)) / R7
% q8 = (T(3) - T_l_t) / R8
% q9 = (T(6) - T_l_t) / R9 
% res_Twg_c    = q1 - q3
% res_Tm_c     = q3 - q6 - q5
% res_Twl_c    = q6 - q8
% 
% res_Twg_l    = q2 - q4
% res_Tm_l     = q4 + q5 - q7
% res_Tfinbase = q7 - q9
% 
% res_global = q1 + q2 - q8 - q9

%% Structural Calcs (copied from SizeRegen1)
% i = 1; 
% 
% while i <= steps % structural calculations
%                 yield(i) = yts(Twg_c(i));
%                 E_current(i) = E(Twg_c(i));
%                 CTE_current(i) = CTE(Twg_c(i));
%                 CTE_liq_side(i) = CTE(Twl_c(i));
%                 elong(i) = elong(Twg_c(i));
%                 epsilon_emax(i) = ((yield(i)*1000000)/ E_current(i));
% 
%                 deltaT1(i) = Twg_c(i) - Twl_c(i);
%                 deltaT2(i) = ((Twg_c(i) + Twl_c(i))/2) - T_l(i); %*** look into what deltaT2 is being used for
% 
%                 sigma_tp(i) = ( ((P_l(i)-P_g(i))/2).*((w_c(i)./t_w(i)).^2) );
%                 sigma_tp_cold(i) =  ( ((P_l(i))/2).*((w_c(i)./t_w(i)).^2) );
%                 sigma_tt(i) = (E_current(i)*CTE_current(i)*qdot_g(i)*t_w_x(i))/(2*(1-nu)*k_w_current(i));
%                 sigma_t(i) = ( ((P_l(i)-P_g(i))/2).*((w_c(i)./t_w(i)).^2) ) + (E_current(i)*CTE_current(i)*qdot_g(i)*t_w_x(i))/(2*(1-nu)*k_w_current(i)); % tangential stress
%                 %sigma_l(i) = E*CTE*(T_wg(i)-T_wl(i)); % longitudinal stress (The temperatures here are wrong and I'm not sure this is applicable to rectagular channels
%                 sigma_lc(i) = E_current(i)*(CTE_liq_side(i)*(T_wl(i)-T_l(i)) + ((CTE_liq_side(i)*deltaT1(i))/(2*(1-nu))));
%                 sigma_ll(i) = E_current(i)*(CTE_current(i)*(T_wg(i)-T_l(i)) + ((CTE_current(i)*deltaT1(i))/(2*(1-nu))));
% 
% 
%                 %sigmab = ??? ; % buckling stress
%                 sigma_vc(i) = sqrt(sigma_lc(i)^2 + sigma_t(i)^2 - sigma_lc(i)*sigma_t(i));
%                 sigma_vl(i) = sqrt(sigma_ll(i)^2 + sigma_t(i)^2 - sigma_ll(i)*sigma_t(i));
% 
%                 % Calculate total Strains
%                 epsilon_lc(i) = ((CTE_liq_side(i)*deltaT1(i))/(2*(1-nu))) + CTE_liq_side(i)*(T_wl(i)-T_l(i));
%                 epsilon_ll(i) = ((CTE_current(i)*deltaT1(i))/(2*(1-nu))) + CTE_current(i)*(T_wg(i)-T_l(i));  
%                 epsilon_tp(i) = ( ((P_l(i)-P_g(i))/2).*((w_c_x(i)./t_w_x(i)).^2) )  /E_current(i);
%                 epsilon_tt(i) = ((E_current(i)*CTE_current(i)*qdot_g(i)*t_w_x(i))/(2*(1-nu)*k_w_current(i)))  /E_current(i);
%                 epsilo      n_t(i) = (( ((P_l(i)-P_g(i))/2).*((w_c_x(i)./t_w_x(i)).^2) ) + (E_current(i)*CTE_current(i)*qdot_g(i)*t_w_x(i))/(2*(1-nu)*k_w_current(i))) / E_current(i); % tangential stress
%                 epsilon_vc(i) = sqrt(epsilon_lc(i)^2 + epsilon_t(i)^2 - epsilon_lc(i)*epsilon_t(i));
%                 epsilon_vl(i) = sqrt(epsilon_ll(i)^2 + epsilon_t(i)^2 - epsilon_ll(i)*epsilon_t(i));
% 
% 
% 
%                 % New structural calcs !?
% 
%                 epsilon_tota(i) = ((CTE_current(i)*deltaT1(i))/(2*(1-nu))) + CTE_current(i) * deltaT2(i); 
%                 epsilon_tott(i) = epsilon_t(i);
%                 epsilon_toteff(i) = (2/sqrt(3)) * sqrt(((epsilon_tott(i)^2)+ epsilon_tott(i)*epsilon_tota(i) + (epsilon_tota(i))^2));
%                 sigma_a(i) = E_current(i) * epsilon_tota(i);
%                 sigma_t2(i) = E_current(i) * epsilon_tott(i);
%                 epsilon_pa(i) = epsilon_tota(i) - epsilon_emax(i);
%                 epsilon_pt(i) = epsilon_tott(i) - epsilon_emax(i);
%                 epsilon_emaxeff(i) = (2/sqrt(3)) * sqrt(((epsilon_emax(i)^2) + epsilon_emax(i)*epsilon_emax(i) + (epsilon_emax(i))^2));
%                 if epsilon_pa(i) < 0
%                     epsilon_pa(i) = 0;
%                 end
%                 if epsilon_pt(i) < 0
%                     epsilon_pt(i) = 0;
%                 end
%                 if epsilon_pa(i) == 0 && epsilon_pt(i) == 0 
%                     epsilon_eeff(i) = epsilon_toteff(i);
%                 elseif epsilon_pa(i) > 0 && epsilon_pt(i) ==  0
%                     epsilon_eeff(i) = (2/sqrt(3)) * sqrt(((epsilon_tott(i)^2) + epsilon_tott(i)*epsilon_emax(i) + (epsilon_emax(i))^2));
%                 elseif epsilon_pa(i) == 0 && epsilon_pt(i) > 0 
%                     epsilon_eeff(i) = (2/sqrt(3)) * sqrt(((epsilon_emax(i)^2) + epsilon_emax(i)*epsilon_tota(i) + (epsilon_tota(i))^2));
%                 else 
%                     epsilon_eeff(i) = epsilon_emaxeff(i);
%                 end
% 
%                 epsilon_peff(i) = (2/sqrt(3)) * sqrt(((epsilon_pt(i)^2) + epsilon_pt(i)*epsilon_pa(i) + (epsilon_pa(i))^2));
%                 epsilon_peff2(i) = epsilon_toteff(i) - epsilon_emaxeff(i);
%                 epsilon_cs(i) = ((elong(i)/2)*((N)^(-1/2)));
%                 %epsilon_cs_tot(i) = epsilon_cs(i) + epsilon_emax(i); % Total strain after N # of hotfires
% 
%                 epsilon_cs_tot(i) = epsilon_cs(i) + epsilon_eeff(i)/2; % Total strain after N # of hotfires
%                 MS_lowcycle(i) = epsilon_cs(i) / (epsilon_peff(i));
%                 MS(i) = epsilon_cs_tot(i) / (epsilon_toteff(i));
%                 num_fires(i) = 1/4 * ((2/elong(i))* (epsilon_toteff(i)- (epsilon_eeff(i)/2)))^(-2);
%                 num_fires_lowcycle(i) = 1/4* ((2*epsilon_peff(i))/elong(i))^(-2);
%                 epsilon_cs_spacex(i) = (3.5 * epsilon_eeff(i)*((N)^(-.12))) + (elong(i)/N)^(.6);
%                 MS_spacex(i) = epsilon_cs_spacex(i) / epsilon_toteff(i);
% 
%                 sigma_eff(i) = E_current(i) * epsilon_toteff(i);
%                 sigma_a(i) = E_current(i) * epsilon_tota(i);
%                 sigma_t2(i) = E_current(i) * epsilon_tott(i);
% 
%                 i = i+1;
% 
%             end

%% Plots

figure;
plot(x_bounds, r_bounds.*u.M2IN, 'bo-', 'LineWidth', 1);
hold on;
plot(x_midpoints, r_midpoints.*u.M2IN, 'bx')
plot(x_bounds, drdx_bounds, 'ro-', 'LineWidth', 1)
plot(x_midpoints, drdx_midpoints, 'rx')
legend('Geometry', 'Derivative')
xlabel('Inches X');
ylabel('Inches Y');
%axis equal;
grid on;
title('Engine Contour and derivative');

% Temperature Plot
figure('Name', 'Temperature Plot');
hold on;
set(gca, 'FontName', 'Times New Roman')
yyaxis left
plot(x_midpoints, Twg_c, 'red', 'LineStyle', '-', 'Marker', 'none');
plot(x_midpoints, Tm_c, 'magenta', 'LineStyle', '-', 'Marker', 'none');
plot(x_midpoints, Twl_c, 'blue', 'LineStyle', '-', 'Marker','none');
plot(x_midpoints, Twg_l, 'red', 'LineStyle', '--', 'Marker','none');
plot(x_midpoints, Tm_l, 'magenta', 'LineStyle', '--', 'Marker', 'none');
plot(x_midpoints, Tfinbase, 'blue', 'LineStyle', '--', 'Marker', 'none');
%plot(x_midpoints, T_l, 'green', 'LineStyle', '-', 'Marker', 'none')
ylabel('Temperature [K]')
set(gca, 'Ycolor', 'k')
grid on
yyaxis right
plot(x_midpoints, r_midpoints, 'black', 'LineStyle', '-');
ylabel('Radius [in]')
set(gca, 'Ycolor', 'k')
axis auto;
legend('Twg (Channel)', 'Tm (Channel)', 'Twl (channel)', 'Twg (land)','Tm (land)','T finbase', 'Chamber Contour', 'Location', 'southoutside', 'Orientation', 'horizontal','Location','best')
title('Temperature Distribution')
xlabel('Location [in]')

% Coolant Flow Plot
figure('Name','Coolant Flow')   
subplot(2,2,[1,2])
hold on 
set(gca, 'FontName', 'Times New Roman')
plot(x_bounds, P_l .* u.PA2PSI)
plot(x_midpoints, P_g .* u.PA2PSI, 'y')
title("Coolant Pressure")
xlabel("Location [in]")
ylabel("Pressure [PSI]")
yyaxis right
plot(x_bounds, r_bounds, 'black', 'LineStyle', '-');
ylabel('Radius [in]')
set(gca, 'Ycolor', 'k')
axis equal;
legend("Pressure Curve"," Chamber Contour",'Location','best')
grid on
subplot(2,2,3)
hold on 
set(gca, 'FontName', 'Times New Roman')
plot(x_midpoints, vel_coolant);
title("Coolant Velocity [m/s]")
xlabel("Location [in]")
grid on
subplot(2,2,4)
hold on 
set(gca, 'FontName', 'Times New Roman')
if coolant == "isopropyl alcohol"
    plot(x_bounds, hmass_l)
else
    plot(liquid_phasex, hmass_liq, 'g', 'DisplayName','liquid' )
    plot(supercrit_phasex, hmass_crit, 'p','DisplayName','supercritical') 
    plot(supercritgas_phasex, hmass_gascrit,'b','DisplayName', 'supercritical gas') 
    plot(supercritliq_phasex, hmass_liqcrit, 'o', 'DisplayName', 'supercritical liquid')
    plot(critpoint_phasex, hmass_critpoint,'br', 'DisplayName', 'critical point')
    plot(gas_phasex, hmass_gas, 'y', 'DisplayName', 'gas')
    plot(two_phasex, hmass_twophase, 'r','DisplayName','twophase');
end
legend


%legend('liquid','supercritical','supercritical gas','supercritical liquid','critical point', 'gas', 'two phase')
title("Coolant Enthalpy [J/kg]");
grid on

% Channel Geometry plot
figure('Name','Channel Geometry');
subplot(2,2,[1,2]);
hold on
set(gca, 'FontName', 'Times New Roman')
plot(x_bounds, w_c_f .*1000);
plot(x_bounds, h_c_f .*1000);
plot(x_bounds, t_w_f .* 1000);
plot(x_bounds, w_land_f .*1000);
%plot(x_contour.*1000, rib_thickness2 .*1000);
title("Channel Dimensions");
xlabel("Location [in]");
ylabel("Channel Dimensions [in]")
yyaxis right
plot(x_bounds, r_bounds, 'black', 'LineStyle', '-');
ylabel('Chamber Contour [in]')
set(gca, 'Ycolor', 'k')
axis equal;
legend('Channel Width', 'Channel Height', 'Wall Thickness', 'Fin Thickness', 'Chamber Contour','Location','northwest')
grid on
subplot(2,2,3);
hold on 
set(gca, 'FontName', 'Times New Roman')
plot(x_midpoints, AR_fin);
title("Fin Aspect Ratio");
xlabel("Location [in]");
grid on
subplot(2,2,4);
hold on 
set(gca, 'FontName', 'Times New Roman')
plot(x_midpoints, AR_channel);
title("Channel Aspect Ratio");
xlabel("Location [in]");
grid on

% Coolant Pressure Plot
figure('Name','Coolant Pressure')   
subplot(2,2,1)
hold on 
set(gca, 'FontName', 'Times New Roman')
plot(x_bounds, P_l_tot .* u.PA2PSI)
title("Coolant Total Pressure")
xlabel("Location [in]")
ylabel("Pressure [PSI]")
grid on
subplot(2,2,2)
hold on 
set(gca, 'FontName', 'Times New Roman')
plot(x_bounds, deltaP_major_cum .* u.PA2PSI);
title("Major Pressure Loss (cumulative)")
xlabel("Location [in]")
grid on
subplot(2,2,3)
hold on 
set(gca, 'FontName', 'Times New Roman')
plot(x_midpoints, deltaP_accel .* u.PA2PSI)
title("Acceleration Pressure Change")
xlabel("Location [in]")
grid on
subplot(2,2,4);
hold on 
set(gca, 'FontName', 'Times New Roman')
plot(x_midpoints, deltaP_minor .* u.PA2PSI)
title("Minor Pressure Loss")
xlabel("Location [in]")
grid on

% Heat Load Plot
figure('Name', 'Heat Transfer Plots');
hold on;
set(gca, 'FontName', 'Times New Roman')
yyaxis left
plot(x_midpoints, q1 ./ 1000, 'red', 'LineStyle', '-', 'Marker', 'none');
plot(x_midpoints, q2 ./ 1000, 'red', 'LineStyle', '--', 'Marker', 'none');
plot(x_midpoints, q3 ./ 1000, 'magenta', 'LineStyle', '-', 'Marker', 'none')
plot(x_midpoints, q4 ./ 1000, 'magenta', 'LineStyle', '--', 'Marker', 'none');
plot(x_midpoints, q5 ./ 1000, 'Color', '#ffa500','LineStyle', '-', 'Marker', 'none');
plot(x_midpoints, q6 ./ 1000, 'green', 'LineStyle', '-', 'Marker', 'none');
plot(x_midpoints, q7 ./ 1000, 'green', 'LineStyle', '--', 'Marker', 'none');
plot(x_midpoints, q8 ./ 1000, 'blue', 'LineStyle', '-', 'Marker', 'none');
plot(x_midpoints, q9 ./ 1000, 'blue', 'LineStyle', '--', 'Marker', 'none');
plot(x_midpoints, q10 ./ 1000, 'Color','#F72F05', 'LineStyle', '-', 'Marker', 'none');
plot(x_midpoints, q11 ./ 1000, 'Color', '#F72F05', 'LineStyle', '-', 'Marker', 'none');
plot(x_midpoints, q12 ./ 1000, 'Color', '#F72F05','LineStyle', '--', 'Marker', 'none');
plot(x_midpoints, q13 ./ 1000, 'Color', '#F72F05','LineStyle', '--', 'Marker', 'none');
ylabel('Heat Flow per hotwall arc length [kW/m]')
set(gca, 'Ycolor', 'k')
grid on
yyaxis right
plot(x_midpoints, r_midpoints, 'black', 'LineStyle', '-');
ylabel('Radius [in]')
set(gca, 'Ycolor', 'k')
legend('q1','q2','q3','q4','q5','q6','q7','q8','q9','q10','q11','q12','q13', 'Chamber Contour','Location','best')
title('Heat Load Distribution')
xlabel('Location [in]')

% Heat Flux Plot
figure('Name', 'Heat Transfer Plots');
hold on;
set(gca, 'FontName', 'Times New Roman')
yyaxis left
plot(x_midpoints, q1_flux ./ 1000, 'red', 'LineStyle', '-', 'Marker', 'x');
plot(x_midpoints, q2_flux ./ 1000, 'red', 'LineStyle', '--', 'Marker', 'x');
plot(x_midpoints, q3_flux ./ 1000, 'magenta', 'LineStyle', '-', 'Marker', 'x')
plot(x_midpoints, q4_flux ./ 1000, 'magenta', 'LineStyle', '--', 'Marker', 'x');
plot(x_midpoints, q5_flux ./ 1000, 'Color', '#ffa500','LineStyle', '-', 'Marker', 'x');
plot(x_midpoints, q6_flux ./ 1000, 'green', 'LineStyle', '-', 'Marker', 'x');
plot(x_midpoints, q7_flux ./ 1000, 'green', 'LineStyle', '--', 'Marker', 'x');
plot(x_midpoints, q8_flux ./ 1000, 'blue', 'LineStyle', '-', 'Marker', 'x');
plot(x_midpoints, q9_flux ./ 1000, 'blue', 'LineStyle', '--', 'Marker', 'x');
ylabel('Heat Flux [kW/m^2]')
set(gca, 'Ycolor', 'k')
grid on
yyaxis right
plot(x_midpoints, r_midpoints, 'black', 'LineStyle', '-');
ylabel('Radius [in]')
set(gca, 'Ycolor', 'k')
legend('q1','q2','q3','q4','q5','q6','q7','q8','q9', 'Chamber Contour','Location','best')
title('Heat Flux Distribution')
xlabel('Location [in]')

%function(x,zeta,eta,sector,r,drdx,kappa,tw,dtwdx,hc,wc)


