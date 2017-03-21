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
Year_t = ncread(source_processed,'Year'); Year_t = reshape(Year_t,[],1);
Month_t = ncread(source_processed,'Month'); Month_t = reshape(Month_t,[],1);
Day_t = ncread(source_processed,'Day'); Day_t = reshape(Day_t,[],1);
Hour_t = ncread(source_processed,'Hour'); Hour_t = reshape(Hour_t,[],1);
Minute_t = ncread(source_processed,'Minute'); Minute_t = reshape(Minute_t,[],1);
Second_t = ncread(source_processed,'Second'); Second_t = reshape(Second_t,[],1);
DateTime_CUP = datetime(Year_t,Month_t,Day_t,Hour_t,Minute_t,Second_t);
vpd = ncread(source_processed,'VPD'); vpd = reshape(vpd,[],1);  
t_air = ncread(source_processed,'Ta'); t_air = reshape(t_air,[],1); 
ws = ncread(source_processed,'Ws'); ws = reshape(ws,[],1);
Rn = ncread(source_processed,'Fn'); Rn = reshape(Rn,[],1);
ea = ncread(source_processed,'e'); ea = reshape(ea,[],1);
P = ncread(source_processed,'ps'); P = reshape(P,[],1);
ET = ncread(source_processed,'ET'); ET = reshape(ET,[],1);
Fe = ncread(source_processed,'Fe'); Fe = reshape(Fe,[],1);
Fh  = ncread(source_processed,'Fh'); Fh = reshape(Fh,[],1);
SWC  = ncread(source_processed,'Sws'); SWC = reshape(SWC,[],1);
Precip  = ncread(source_processed,'Precip'); Precip = reshape(Precip,[],1);
% clear unused variables
clearvars source_processed finfo;

% Load raw data (.csv file)
% Location of raw input file
source_raw = 'Input\CUP_EddyPro_QC_170207.csv';
% Create datastore to access collection of data
ds_raw = datastore(source_raw);
% Select variable of interest
ds_raw.SelectedVariableNames = {'DateTime','NEE_c','qc_Sc','AGC_c','daytime','qc_co2_flux','qc_h2o_flux','u_','h2o_flux'}; % for example
% Read selected variables, save it in the workspace as a table
Raw_Data = readall(ds_raw);
% Get data from the table, change format if necessary
NEE_c = Raw_Data.NEE_c;
qc = Raw_Data.qc_co2_flux;
qc_h2o_flux = Raw_Data.qc_h2o_flux;
h2o_flux = Raw_Data.h2o_flux;
qc_Sc = Raw_Data.qc_Sc;
AGC_c = Raw_Data.AGC_c;
daytime = Raw_Data.daytime;
u = Raw_Data.u_;
DateTime_CUP_cell = Raw_Data.DateTime;
t_s = datetime(DateTime_CUP_cell,'InputFormat','dd/MM/yyyy HH:mm');
% clear unused variables
clearvars DateTime_CUP_cell Raw_Data ds_raw source_raw; 

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
pa = (P.*1000)./(287.058.*(t_air+237.3)); 

% P is air pressure, loaded from met data, [kPa]

% psv is the temperature dependent latent heat of vaporization, [J/kg]
psv = (2.501-0.00237.*t_air).*1000000;

% cp is the specific heat capacity for dry air, [J kg-1 K-1]
cp = 1004.834; 

% ps is the temperature-dependent psychrometric constant, [kPa K-1]
ps = (cp.*P)./(0.622.*psv);

% ga is the aerodynamic conductance, [m/s]
k = 0.40; % Von Karman Constant [unitless]
h = 20; % canopy height [m]
z = 30; % measurement height [m]
ga = (ws.*(k.^2))./((log((z-(0.67).*h)./((0.1).*h))).^2);

% Canopy conductance [m s-1] from Penman-Monteith [Monteith, 1965]
% Eq. found in Knauer 2015 -epsilom is s
gc = (ps.*Fe.*ga)./(s.*Rn + pa.*cp.*vpd.*ga - Fe.*(s+ps));


% Novick 2016 supplement: "zm is the measurement height,
% zd is the zero plane displacement, and zo is the momentum roughness length. The zd and zo where taken as
% 0.67h and 0.1h, respectively, where h is canopy height, as is common practice in the absence of other
% information about these parameters
% gs is the highest value of the reference,
% well-watered surface conductance rate observed
% across the study domain", 

% filter after rain event
% Rain in the past 2 days, 1 day and 12 hour
n = length(Precip);
Precip_2daysago = NaN(n,1); % 96 half-hour (2 days)
Precip_1daysago = NaN(n,1); % 48 half-hour (1 day)
Precip_halfdayago = NaN(n,1); % 24 half-hour (12 hour)
for i = 97:n
   Precip_2daysago(i) = sum(Precip(i-96:i)); 
end
for i = 49:n
   Precip_1daysago(i) = sum(Precip(i-48:i)); 
end
for i = 25:n
   Precip_halfdayago(i) = sum(Precip(i-24:i)); 
end

% rain in the past 2 days lower than 1mm
% rain in the past 24 hour lower than 0.5 mm
% rain in the past 12 hour lower than 0.2 mm
use = find(Precip_2daysago < 1 & Precip_1daysago < 0.5 & Precip_halfdayago < 0.2); 
% 2014-2016 at CUP, this is 60% of the data

% 6 quantiles of SWC for rain filtered data
binedges_f = quantile(SWC(use),0:1/6:1);
% 6 quantiles of SWC for all data
binedges_nof = quantile(SWC,0:1/6:1);
% ref,ww condition for rain filtered data
use_f = find(Precip_2daysago < 1 & Precip_1daysago < 0.5 & Precip_halfdayago < 0.2 & vpd > 0.9 & vpd < 1.3 & qc_h2o_flux == 0 & AGC_c == 0 & u > 0.2 & daytime == 1 & SWC > binedges_f(6) & gc > 0);
% ref,ww condition for all data
use_nof = find(vpd > 0.9 & vpd < 1.3 & qc_h2o_flux == 0 & AGC_c == 0 & u > 0.2 & daytime == 1 & SWC > binedges_nof(6) & gc > 0);
% gs,ref,ww for rain filtered data
gs_f = mean(gc(use_f));
% gs,ref,ww for all data
gs_nof = mean(gc(use_nof));


% gs = gc when AGC is good, qc is good, VPD > 0.9 & VPD < 1.1, wet soil
% (quantile > 0.8 maybe), 
% gs_f: days with rainfall and the two subsequent
% dayswere excluded if precipitation exceeded 0.2mm(day with rainfall), 0.5mm(day before), or 1mm
% (2 days before).

% PET calculation, penman monteith
PET = (s.*Rn+cp.*pa.*ga.*vpd)./(psv.*(s+ps.*(1+(ga./gs_nof)))); % [unit?]
PET = PET*3600; % from [mm s-1] to [mm]
% PET_pt, Priestley-Taylor
PET_pt = 1.26.*(((s.*Rn)./(s+ps)).*(1./psv));
PET_pt = PET_pt*3600; % from [mm s-1] to [mm]

% if PET is [mm s-1], to go to [mm] we simply do *3600 second

PET_m_ET = PET-ET;
DI = sum(PET)/sum(Precip);
% figure; scatter(vpd,PET_m_ET);

figure; hold on;
for s = 1:2
    if s == 1
        this_season = find(Month_t >= 5 & Month_t <= 8 & Precip_2daysago < 1 & Precip_1daysago < 0.5 & Precip_halfdayago < 0.2 & vpd > 0.9 & qc_h2o_flux == 0 & AGC_c == 0 & u > 0.2 & daytime == 1 & gc > 0); % winter daytime
    elseif s == 2
        this_season = find((Month_t >= 10 | Month_t <= 2) & Precip_2daysago < 1 & Precip_1daysago < 0.5 & Precip_halfdayago < 0.2 & vpd > 0.9 & qc_h2o_flux == 0 & AGC_c == 0 & u > 0.2 & daytime == 1 & gc > 0); % summer daytime
    end
    SWC_ss = SWC(this_season);
    VPD_s = vpd(this_season);
    gc_relative_s = gc(this_season)/gs_f;
    SWC_s_binedges = quantile(SWC_ss,0:1/3:1); % 3 quantiles of SWC
    c = 0.8; % dark blue/red is wet soil
    for j = 1:3 % SWC
        if s == 1
            colors = [c c 1];
        elseif s == 2
            colors = [1 c c];
        end
        c = c - 0.4;
        this_bin = find(SWC_ss >= SWC_s_binedges(j) & SWC_ss <= SWC_s_binedges(j+1)); %* iteration change SWC
        binplot(VPD_s(this_bin),gc_relative_s(this_bin),4,colors);
    end
end
ax = gca; ax.FontSize = 16; xlim([0.5 4.5]); ylim([0.1 1.1]); ax.YTick = 0.1:0.2:1.1; ax.XTick = 0.5:1:4.5;
xlabel('VPD (kPa)');
ylabel('Relative surface conductance');


clearvars n i 
