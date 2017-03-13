% PET [unit?] and gs (canopy conductance) [unit?] calculation
% I used the supplement info in Novick et al 2016
% I used eq. 11 in Knauer et al 2015 JGR Biogeosciences to infer canopy
% conductance

% I use Penman Monteith to infer PET and gs. 

% PET = (s.*Rn+cp.*pa.*ga.*vpd)./(psv.*(s+ps.*(1+(ga./gs)))) % Unit: [mm s-1]
% parameters are described below

% Input
% P, air pressure, [kPa]
% ea, vapor pressure of the air, [kPa]
% Rn, net radiation, [W/m2]
% t_air, air temperature, [°C]
% ws, wind speed, [m/s]
% vpd, vapor pressure deficit, [kPa]

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Load input data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear; % initialise workspace
% Load processed data (netCDF file)
% Location of processed input file 
source_processed = 'Input\CumberlandPlain_2014_L6_EP_moderate.nc';
% Information on file, including variable name
finfo = ncinfo(source_processed);
% Read variable from .cd file
year_t = ncread(source_processed,'Year'); year_t = reshape(year_t,[],1);
month_t = ncread(source_processed,'Month'); month_t = reshape(month_t,[],1);
day_t = ncread(source_processed,'Day'); day_t = reshape(day_t,[],1);
hour_t = ncread(source_processed,'Hour'); hour_t = reshape(hour_t,[],1);
minute_t = ncread(source_processed,'Minute'); minute_t = reshape(minute_t,[],1);
second_t = ncread(source_processed,'Second'); second_t = reshape(second_t,[],1);
DateTime_CUP = datetime(year_t,month_t,day_t,hour_t,minute_t,second_t);
% SWC_s = ncread(source_processed,'Sws'); SWC_s = reshape(SWC_s,[],1); SWC_s = SWC_s*100;
vpd = ncread(source_processed,'VPD'); vpd = reshape(vpd,[],1); 
% tsoil = ncread(source_processed,'Ts'); tsoil = reshape(tsoil,[],1); 
t_air = ncread(source_processed,'Ta'); t_air = reshape(t_air,[],1); 
ws = ncread(source_processed,'Ws'); ws = reshape(ws,[],1);
Rn = ncread(source_processed,'Fn'); Rn = reshape(Rn,[],1);
ea = ncread(source_processed,'e'); ea = reshape(ea,[],1);
P = ncread(source_processed,'ps'); P = reshape(P,[],1);
ET = ncread(source_processed,'ET'); ET = reshape(ET,[],1);
Fe = ncread(source_processed,'Fe'); Fe = reshape(Fe,[],1);
Fh  = ncread(source_processed,'Fh'); Fh = reshape(Fh,[],1);
% clear unused variables
clearvars source_processed finfo;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Data processing and analysis
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% s is the temperature-dependent slope of the saturation-vapor pressure
% curve, [kPa/K)
% saturation vapor pressure, [kPa]
es = 0.61078.*exp((17.269.*t_air)./(237.3+t_air));
% s, [kPa/K]
s = es.*(17.269./(237.3+t_air)).*(1-(t_air./(237.3+t_air)));

% Rn is net radiation, loaded from measured met data, [W/m2]

% pa is the density of dry air, [kg/m3]
pa = 1.2041; % at 20 °C and 101 kPa -do I need a t_air dependent pa?

% P is air pressure, loaded from met data, [kPa]

% psv is the temperature dependent latent heat of vaporization, [J/kg]
psv = 2.5005.*(10.^6) - 2.359.*(10.^3).*(t_air+237.15);

% cp is the specific heat capacity for dry air, [J kg K-1]
cp = 1004.7; 

% ps is the temperature-dependent psychrometric constant, [kPa K-1]
ps = (1.61.*cp.*P)./psv;

% ga is the aerodynamic conductance, [m/s]
k = 0.40; % Von Karman Constant [unitless]
h = 20; % canopy height [m]
z = 30; % measurement height [m]
ga = (ws.*(k.^2))./((log((z-(0.67).*h)./((0.1).*h))).^2);

% Canopy conductance [m s-1] from Penman-Monteith [Monteith, 1965]
% Eq. found in Knauer 2015 -epsilom is s
gc = (ps.*Fe.*ga)./(s.*Rn + pa.*cp.*vpd.*ga - Fe.*(s+psv));


% Novick 2016 supplement: "zm is the measurement height,
% zd is the zero plane displacement, and zo is the momentum roughness length. The zd and zo where taken as
% 0.67h and 0.1h, respectively, where h is canopy height, as is common practice in the absence of other
% information about these parameters
% gs is the highest value of the reference,
% well-watered surface conductance rate observed
% across the study domain", 
gs = 0.0148; % [m/s] for now... I need to know how to calculate it
% PET calculation, penman monteith
PET = (s.*Rn+cp.*pa.*ga.*vpd)./(psv.*(s+ps.*(1+(ga./gs)))); % [unit?]
PET = PET*3600; % from [mm s-1] to [mm]
% PET_pt, Priestley-Taylor
PET_pt = 1.26.*(((s.*Rn)./(s+ps)).*(1./psv));
PET_pt = PET_pt*3600; % from [mm s-1] to [mm]

% if PET is [mm s-1], to go to [mm] we simply do *3600 second

PET_m_ET = PET-ET;
DI = sum(PET)/sum(ET);
% figure; scatter(vpd,PET_m_ET);
