%% sizeRegen2 throat calc only
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
addpath("cea\");
addpath("IPA Data\");
addpath("Material Data\");
Data = LoadData();
u = convertUnits;
%% Engine Inputs 

% Propellants
P_c = 250; 
OF = 1.2;
mdot_coolant = 1.5 * u.LB2KG; %lbm/s
%mdot_coolant = 0.997305579515822;
fuel = 0; % Choose fuel propellant
oxidizer = 0; % Choose Ox propellant 
fuel_weight = 100;%[75,25]; % weight of each fuel used (if only 1 fuel, set to 100)  
fuel_temp = 293.15;%[293.15,293.15]; % [K]
oxidizer_temp = 90.17; % [K]
oxidizer_weight = 100; %

% Coolant Parameters
%coolant = 'methanol'; % Choose coolant fluid
coolant = 'water';
%coolant = "isopropyl alcohol"
coolant_direction = 0;
inlet_temp = 273.15; % Initial Coolant Temperature [K]
inlet_pressure = 500; % Initial Coolant Pressure [psi]

% Chamber/Nozzle Geometry
R_t = 0.73217;
exp_ratio = 2.9525;
con_ratio = 6.558;
conv_angle = 37.5;
conical_half_angle = 15; 
L_star = 45; 
R_fillets = [1.5,1.5,.5];

% Channel Geometry
num_channels = 62;
t_w = [.001 .00075 .00075];
h_c = [.002 .0017 .0023]; 
w_c = [.002 .0005 .0012]; %arc length at the channel's mid radius 



% Material & Manufacturing Method
material = 1; %0: 6061-RAM2 1: 718 inco

%% Material Properties
if material == 0
    properties = readmatrix(pwd + "/Material Data/al6061_RAM2.xlsx");
    k_w = properties(16:24,1:2); % thermal conductivity [W/m-K]
    E = [properties(16:21, 8) properties(16:21,9)];
    CTE = [properties(16:20,4) properties(16:20,5)]; % [m/m*K]
    yield_strength = properties(1:11,1:2);
    elongation_break = [properties(1:11,1) properties(1:11,4)];
elseif material == 1 
    properties = readmatrix(pwd + "/Material Data/Inconel718.xlsx");
    k_w = properties(13:end,1:2); % thermal conductivity [W/m-K]
    yield_strength = properties(1:7,1:2);
    elongation_break = [properties(1:7,1) properties(1:7,5)];
    E = [properties(1:6, 9) properties(1:6,10)]; %youngs modulus
    CTE = [properties(1:7,1) properties(1:7,3)]; % [ppm]
elseif material == 2
    properties = readmatrix(pwd + "/Material Data/GrCop42.xlsx");
    k_w = properties(13:end,1:2); % thermal conductivity [W/m-K] 
    E = [properties(1:13, 9) properties(1:13,10)];
    CTE = [properties(1:5,1) properties(1:5,3)]; % [ppm] 
    yield_strength = properties(1:8,1:2);
    elongation_break = [properties(1:8,1) properties(1:8,5)];
elseif material == 3
    properties = readmatrix(pwd + "/Material Data/Inconel625.xlsx");
    k_w = properties(8:19,1:2); % thermal conductivity [W/m-K] 
    E = [properties(1:5, 1) properties(1:5,5)];
    CTE = [properties(22:end,1:2)]; % [ppm] 
    yield_strength = properties(1:5,1:2);
    elongation_break = [properties(1:5,1) properties(1:5,3)];
elseif material == 4
    properties = readmatrix(pwd + "/Material Data/C101.xlsx");
    k_w = properties(13:21,1:2); % thermal conductivity [W/m-K] 
    E = [properties(1:5, 9:10)];
    CTE = [properties(1:7,1) properties(1:7,3)]; % [ppm] 
    yield_strength = properties(1:7, 1:2);
    elongation_break = [properties(1:7,1) properties(1:7,5)];
elseif material == 5
    properties = readmatrix(pwd + "/Material Data/AlSi10Mg.xlsx");
    k_w = properties(13:15,1:2); % thermal conductivity [W/m-K] 
    E = [properties(1:6, 9:10)];
    CTE = [properties(1:4,1) properties(1:4,3)]; % [ppm] 
    yield_strength = properties(1:4,1:2);
    elongation_break = [properties(1:4,1) properties(1:4,5)];
end

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


D_t = 2* R_t * u.IN2M; % Throat diameter [m]
R_of_curve_up = R_fillets(2) * D_t/2;
R_of_curve_down = R_fillets(3) * D_t/2;
R_of_curve =  (R_of_curve_up + R_of_curve_down) / 2; % Avg Radius of  curvature of the throat [m] **This is what is used in the Huzel & Huang example problems**

%% Generate Chamber / Nozzle Contour

% Define the bounds of each cell
geom = engineContour2(R_t, exp_ratio, con_ratio, conv_angle, conical_half_angle, L_star, R_fillets);
x_chamber_bounds = linspace(geom.bounds.chamber(1),geom.bounds.chamber(2),25);
x_converging_bounds = linspace(geom.bounds.converging(1),geom.bounds.converging(2),100);
x_diverging_bounds = linspace(geom.bounds.diverging(1),geom.bounds.diverging(2),100);
x_bounds = [x_chamber_bounds, x_converging_bounds(2:end), x_diverging_bounds(2:end)];
r_bounds = geom.r(x_bounds);
drdx_bounds = geom.drdx(x_bounds);

% Define the midpoints of each cell 
i = 1; 
x_midpoints = ones(1, length(x_bounds)-1);
r_midpoints = ones(1, length(x_bounds)-1);
drdx_midpoints = ones(1,length(x_bounds)-1);
while i < length(x_bounds)
    x1 = x_bounds(i); % Front Boundary of Cell
    x2 = x_bounds(i+1); % Back Boundary of Cell
    ds = integral(@(x) sqrt(1 + geom.drdx(x).^2), x1, x2); % Arc Length of Cell
    midpointResidual = @(xm) integral(@(x) sqrt(1 + geom.drdx(x).^2),x1,xm) - ds/2; % Equation whose root is the arc-length midpoint
    x_mid = fzero(midpointResidual, [x1, x2]);
    r_mid = geom.r(x_mid);
    drdx_mid = geom.drdx(x_mid);

    x_midpoints(i) = x_mid; 
    r_midpoints(i) = r_mid;
    drdx_mid(i) = drdx_mid;

    i = i+1;
end

%% Get Combustion Properties
% Use CEA for this 
% Currently uses Kamons RunCEA wrapper. Might want to change to yash's
% callCEA wrapper
% Get throat Properties
CEA_input_name = 'regenCEAThroat';
%[~, ~, ~, M_t, gamma_t, P_g_t, T_g_t, rho_g_t, mu_g_t, Pr_g_t, ~, ~, son_g_t, cp_g_t] = RunCEA(P_c, P_e, fuel, fuel_weight, fuel_temp, oxidizer, oxidizer_temp, OF_tot, 0, 1, 2, 0, 0, CEA_input_name);
[CEAdata] = callCEA('fr', P_c, 'psi', 'o/f', OF, 'sup', 1, fuel_def, 'K', fuel_temp, fuel_weight, oxidizer_def, 'K', oxidizer_temp, oxidizer_weight);
T_g_out = squeeze(CEAdata('t')); 
T_g_t = T_g_out(2);
T_g_tot = T_g_out(1);
cp_g_out = squeeze(CEAdata('cp')); 
cp_g_t = cp_g_out(2);
cp_g_tot = cp_g_out(1);
k_g_out = squeeze(CEAdata('k')); 
k_g_t = k_g_out(2); 
M_out = squeeze(CEAdata('mach'));
M_t = M_out(2);
cf_g_out = squeeze(CEAdata('cf')); 
cf_g_t = cf_g_out(2);
Pr_g_out = squeeze(CEAdata('prandtl')); 
Pr_g_t = Pr_g_out(2); 
Pr_g_tot = Pr_g_out(1);
gamma_out = squeeze(CEAdata('gammas'));
gamma_t = gamma_out(2);
mu_g_out = squeeze(CEAdata('visc'));
mu_g_t = mu_g_out(2); 
mu_g_tot = mu_g_out(1);
cstar_out = squeeze(CEAdata('cstar'));
c_star = cstar_out(1);

%% Throat Calc w/ Picard solution

P_c = P_c * u.PSI2PA;

w_c_t = w_c(2) ; 
h_c_t = h_c(2) ; % THis needs to be calculated in the geo section based on the arc length
t_w_t = t_w(2); 

r_wg_t = D_t/2;
r_wl = r_wg_t + t_w_t;
r_outer = r_wl + h_c_t;
r_midchan = 0.5* (r_wl + r_outer);
theta_channel = w_c_t /r_midchan; % radians
A_channel_t = theta_channel/2  * (r_outer^2 - r_wl^2);
theta_p = 2*pi()/ num_channels; 
theta_land = theta_p - theta_channel; 
r_m_t = (r_wg_t + r_wl)/2;
mdot_channel = mdot_coolant / num_channels;
P_wet_t = (r_wl + r_outer)*theta_channel + 2*(r_outer - r_wl);
D_hydrualic_t = 4 * A_channel_t / P_wet_t;



T_l_t = 290; % Coolant temp [K]
P_l_t = 400 * u.PSI2PA; % Coolant Pressure [psi to pa]
heatflux_factor = 1;
fin_efficiency_factor = 1;
coolant_heatflux_factor = 1;
nolatconductance = 0;
Twg_c = 700;
Twl_c = 300;
Tm_c = (Twg_c + Twl_c) ./ 2;  
Twg_l = 700;
Tfinbase = 300;
Tm_l = (Twg_l + Tfinbase) ./ 2; 
T = [700; 500; 300; 700; 500; 300]; %Initial Temperature Guesses

r = Pr_g_t ^ (1 / 3);
T_r = T_g_t * (1 + (gamma_t - 1) / 2 * r * M_t ^ 2); % recovery temperature [K] - corrects for compressible boundry layers (Heister EQ 6.15). 

% Try out the numerical solution at the throat 
converged = 0;
tolerance = .1;
while ~converged
    % Old guesses
    T_old = T;
    Twg_c = T_old(1);
    Tm_c = T_old(2);
    Twl_c = T_old(3);  
    Twg_l = T_old(4);
    Tm_l = T_old(5);
    Tfinbase = T_old(6); 

    % Resistance 1, Hot gas convection to channel hot wall
    sigma_c = (.5 * Twg_c / T_g_tot * (1 + (gamma_t - 1) / 2 * M_t ^ 2) + .5) ^ -.68 * (1 + (gamma_t - 1) / 2 * M_t ^ 2) ^ -.12; % film coefficient correction factor [N/A] (Huzel & Huang 86).
    h_g_c = heatflux_factor * (0.026 / D_t ^ 0.2) * (mu_g_tot ^ 0.2 * cp_g_tot / Pr_g_tot ^ 0.6) * (P_c / c_star) ^ 0.8 * (D_t / R_of_curve) ^ 0.1 * (1 / 1) ^ .9 * sigma_c; % gas film coefficient [W/m^2-K] - bartz equation (Huzel & Huang 86).
    %R1 = 1 / ( h_g_c * w_c_t);
    R1 = 1 / (h_g_c * theta_channel * r_wg_t);


    % Resistance 2: Hot gas convection to land hot wall
    sigma_l = (.5 * Twg_l / T_g_tot * (1 + (gamma_t - 1) / 2 * M_t ^ 2) + .5) ^ -.68 * (1 + (gamma_t - 1) / 2 * M_t ^ 2) ^ -.12; % film coefficient correction factor [N/A] (Huzel & Huang 86).
    h_g_l = heatflux_factor * (0.026 / D_t ^ 0.2) * (mu_g_tot ^ 0.2 * cp_g_tot / Pr_g_tot ^ 0.6) * (P_c / c_star) ^ 0.8 * (D_t / R_of_curve) ^ 0.1 * (1 / 1) ^ .9 * sigma_l; % gas film coefficient [W/m^2-K] - bartz equation (Huzel & Huang 86).
    %R2 = 1 / ( h_g_l * w_l_t); 
    R2 = 1 / (h_g_l * theta_land * r_wg_t);

    % Resistance 3: 
    k3 = interp1(k_w(:,1), k_w(:,2), (Twg_c+Tm_c)/2, 'linear', 'extrap');
    %R3 = (t_w_t/2) / (k3 * w_c_t);
    R3 = log(r_m_t / r_wg_t) / (k3 * theta_channel);

    % Resisitance 4: 
    k4 = interp1(k_w(:,1), k_w(:,2), (Twg_l+Tm_l)/2, 'linear', 'extrap');
    %R4 = (t_w_t/2) / (k4 * w_l_t);
    R4 = log(r_m_t / r_wg_t) / (k4 * theta_land);

    % Resistance 5:
    k5 = interp1(k_w(:,1), k_w(:,2), (Tm_c+Tm_l)/2, 'linear', 'extrap');
    %R4 = ((w_c_t + w_l_t)/2) / (k5 * t_w_t);
    R5 = (theta_p/2) / (2*k5 * log(r_wl/r_wg_t));
    
    % Resistance 6:
    k6 = interp1(k_w(:,1), k_w(:,2), (Tm_c+Twl_c)/2, 'linear', 'extrap');
    %R6 = (t_w_t/2) / (k6 * w_c_t);
    R6 = log(r_wl / r_m_t) / (k6 * theta_channel);

    % Resistance 7:
    k7 = interp1(k_w(:,1), k_w(:,2), (Tm_l+Tfinbase)/2, 'linear', 'extrap');
    %R7 = (t_w_t/2) / (k7 * w_l_t);
    R7 = log(r_wl / r_m_t) / (k7 * theta_land);

    % Resistance 8: 
    % Cold wall convection 
    if coolant == "isopropyl alcohol"
        P_kPa_t = P_l_t / 1000;
        T_K_t = T_l_t;
        
        rho_l_t = Data.F_rho(T_K_t, P_kPa_t);
        cp_l_t  = Data.F_cp(T_K_t, P_kPa_t) * 1000; 
        mu_data = Data.F_mu_data(T_K_t);
        k_data  = Data.F_k_data(T_K_t);
        
        if (T_K_t - 273.15) < 36.5
            mu_fit = 4.5054 * exp(-0.031 * (T_K_t - 273.15)) / 1000;
        else
            mu_fit = 3.3724 * exp(-0.023 * (T_K_t - 273.15)) / 1000;
        end
            
        mu_avg = (mu_fit + mu_data) / 2;
        mu_l_t = mu_avg * exp(5.4574e-9 * (P_l_t - 101325));
        k_l_t  = k_data * (1 + 2.9531e-9 * (P_l_t - 101325));
    else
        mu_l_t = double(py.CoolProp.CoolProp.PropsSI('V','T', T_l_t, 'P', P_l_t, coolant)); % viscosity of bulk coolant [Pa-s]
        cp_l_t = double(py.CoolProp.CoolProp.PropsSI('C' , 'T', T_l_t, 'P', P_l_t, coolant)); % specific heat of coolant [J/kg-k] 
        k_l_t = double(py.CoolProp.CoolProp.PropsSI('L', 'T', T_l_t, 'P', P_l_t, coolant)); % thermal conductivity of coolant [W/m-K]
        rho_l_t = double(py.CoolProp.CoolProp.PropsSI('D','T', T_l_t,'P', P_l_t, coolant)); % density of the coolant [kg/m^3]
    end
    
    vel_coolant = mdot_channel / rho_l_t / A_channel_t;
    Re_l_t = (rho_l_t * vel_coolant * D_hydrualic_t) / mu_l_t; % reynolds number for channel flow [N/A] (Huzel and Huang , pg 90)
    Pr_l_t = (cp_l_t * mu_l_t) / k_l_t; % prantl number [N/A] (Huzel and Huang, pg 90) 

    % Use moody diagram to find coefficient of friction
    ed = e(2)/(D_hydrualic_t*1000); % relative roughness
    f = moody(ed, Re_l_t); % friction factor

    %Nu_l(i) = 0.023 * (Re_l(i) ^ .8) * (Pr_l(i) ^ .4) * (T_wl(i) / T_l(i)) ^ -.3; % nusselt number [N/A] - applicable for Re > 10,000, .7 < Pr < 160 (Heister EQ 6.19). ****
    Nu_l_t = (f / 8) * (Re_l_t - 1000) * Pr_l_t / (1 + 12.7 * (f / 8) ^ 0.5 * (Pr_l_t ^ (2/3) - 1)); % Gnielinksy correlation nusselt number [N/A] - 0.5 < Pr < 2000, 3000 < Re < 5e6
    h_l_t = coolant_heatflux_factor *  (Nu_l_t * k_l_t) / D_hydrualic_t; % liquid film coefficient [W/m^2-K] (Heister EQ 6.19)

    R8 = 1 /(h_l_t * theta_channel * r_wl);

    % Resistance 9 
    k9 = interp1(k_w(:,1), k_w(:,2), Tfinbase, 'linear', 'extrap');
    r_midfin = 0.5*(r_wl + r_outer);
    w_land_midfin = r_midfin * theta_land;
    m_fin = sqrt(2*h_l_t / (k9*w_land_midfin)); % Constant fin width approximation use midwidth of fin
    eta_fin = fin_efficiency_factor * tanh(m_fin * (r_outer - r_wl)) / (m_fin*(r_outer - r_wl));  % Fin efficiency
    R9 = 1/ ( eta_fin * h_l_t * 2 * (r_outer-r_wl)); 

    % Convert resistances to Conductances
    G1 = 1/R1; 
    G2 = 1/R2;
    G3 = 1/R3;
    G4 = 1/R4; 
    G5 = 1/R5;
    G6 = 1/R6; 
    G7 = 1/R7; 
    G8 = 1/R8; 
    G9 = 1/R9; 
    if nolatconductance 
        G5 = 0;
    end

    % Construct Conductance Matrix
    A = [-(G1+G3),          G3,        0,        0,           0,        0;
               G3, -(G3+G5+G6),       G6,        0,          G5,        0; 
                0,          G6, -(G6+G8),        0,           0,        0;
                0,           0,        0, -(G2+G4),          G4,        0;
                0,          G5,        0,       G4, -(G4+G5+G7),       G7;
                0,           0,        0,        0,          G7, -(G7+G9) ];
    b = [-G1*T_r; 0; -G8*T_l_t; -G2*T_r; 0; -G9*T_l_t];

    % Solve for T 
    T = A \ b; 
    

    % Check convergence 
    error_max = max(abs(T-T_old))
    if error_max < tolerance
        converged = true;
    end

end

h = T
q1 = (T_r - T(1)) / R1
q2 = (T_r - T(4)) / R2
q3 = (T(1) - T(2)) / R3
q4 = (T(4) - T(5)) / R4
q5 = (T(2) - T(5)) / R5
q6 = (T(2) - T(3)) / R6
q7 = (T(5) - T(6)) / R7
q8 = (T(3) - T_l_t) / R8
q9 = (T(6) - T_l_t) / R9 
res_Twg_c    = q1 - q3
res_Tm_c     = q3 - q6 - q5
res_Twl_c    = q6 - q8

res_Twg_l    = q2 - q4
res_Tm_l     = q4 + q5 - q7
res_Tfinbase = q7 - q9

res_global = q1 + q2 - q8 - q9
%% Plots
figure;
plot(x_bounds, r_bounds, 'bo-', 'LineWidth', 1);
hold on;
plot(x_midpoints, r_midpoints, 'bx')
plot(x_bounds, drdx_bounds, 'ro-', 'LineWidth', 1)
plot(x_midpoints, drdx_midpoints, 'rx')
legend('Geometry', 'Derivative')
xlabel('Inches X');
ylabel('Inches Y');
axis equal;
grid on;
title('Engine Contour and derivative');
