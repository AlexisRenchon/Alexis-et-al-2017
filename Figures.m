% All figures

% Load input

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
tsoil = ncread(source_processed,'Ts'); tsoil = reshape(tsoil,[],1); 
ws = ncread(source_processed,'Ws'); ws = reshape(ws,[],1);
Rn = ncread(source_processed,'Fn'); Rn = reshape(Rn,[],1);
ea = ncread(source_processed,'e'); ea = reshape(ea,[],1);
P = ncread(source_processed,'ps'); P = reshape(P,[],1);
GEP = ncread(source_processed,'GPP_SOLO_EP'); GEP = reshape(GEP,[],1); 
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

% PAR data
% Load raw data (.csv file)
% Location of raw input file
source_PAR = 'Input\FACELawn_FACE_diffPAR_20142017_clean.csv';
% Create datastore to access collection of data
ds_PAR = datastore(source_PAR);
% Select variable of interest
ds_PAR.SelectedVariableNames = {'DateTime','PAR_Avg_mean'}; 
% Read selected variables, save it in the workspace as a table
PAR_Data = readall(ds_PAR);
% Get data from the table, change format if necessary
PAR = PAR_Data.PAR_Avg_mean;
DateTime_PAR = PAR_Data.DateTime;
DateTime_PAR = datetime(DateTime_PAR,'InputFormat','dd/MM/yyyy HH:mm');
% clear unused variables
clearvars PAR_Data ds_PAR source_PAR;

% Daily data from Summary .xls
source_summary = 'Input\CumberlandPlain_2014_L6_EP_moderate_Summary.xlsx';
Data_sum = readtable(source_summary,'Sheet','Daily (all)','Range','A3:AL1097','ReadVariableNames',false);
Date_daily = Data_sum.Var1; Date_daily = datetime(Date_daily,'InputFormat','dd/MM/yyyy hh:mm:ss a');
ET_daily = Data_sum.Var8;
ER_d = Data_sum.Var7;
GPP_d = Data_sum.Var20;
NEE_d = Data_sum.Var23;
Month_d = Data_sum.Var36;
clearvars Data_sum source_summary;


% LAI data
% Load raw data (.csv file)
% Location of raw input file
source_LAI = 'Input\CUP_LAI_20170127.csv';
% Create datastore to access collection of data
ds_LAI = datastore(source_LAI);
% Select variable of interest
ds_LAI.SelectedVariableNames = {'Date','LAI'}; 
% Read selected variables, save it in the workspace as a table
LAI_Data = readall(ds_LAI);
% Get data from the table, change format if necessary
LAI = LAI_Data.LAI;
DateTime_LAI = LAI_Data.Date;
DateTime_LAI = datetime(DateTime_LAI,'InputFormat','yyyy-MM-dd');
% clear unused variables
clearvars LAI_Data ds_LAI source_LAI; 

% Script

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

% PAR saturated (> 800) ET, GEP and WUE vs VPD, for SWC quantiles (n quantiles)
figure; hold on;
for s = 1:2 % 2 seasons, winter or summer
    if s == 1
        %this_season = find(Precip_2daysago < 1 & Precip_1daysago < 0.5 & Precip_halfdayago < 0.2 & Month_t >= 5 & Month_t <= 8 & daytime == 1 & qc == 0 & qc_h20_flux == 0 & u > 0.2 & AGC_c == 0 & qc_Sc == 0 & NEE_c > -20 & NEE_c < 15 & PAR > 1000); % winter daytime
        this_season = find(Precip_2daysago < 1 & Precip_1daysago < 0.5 & Precip_halfdayago < 0.2 & Month_t >= 5 & Month_t <= 8 & daytime == 1 & qc == 0 & qc_h2o_flux == 0 & u > 0.2 & AGC_c == 0 & qc_Sc == 0 & NEE_c > -20 & NEE_c < 15 & PAR > 800); % winter daytime
    elseif s == 2
        %this_season = find(Precip_2daysago < 1 & Precip_1daysago < 0.5 & Precip_halfdayago < 0.2 & (Month_t >= 10 | Month_t <= 2) & daytime == 1 & qc == 0  & qc_h20_flux == 0 & u > 0.2 & AGC_c == 0 & qc_Sc == 0 & NEE_c > -20 & NEE_c < 15 & PAR > 1000); % summer daytime
        this_season = find(Precip_2daysago < 1 & Precip_1daysago < 0.5 & Precip_halfdayago < 0.2 & (Month_t >= 10 | Month_t <= 2) & daytime == 1 & qc == 0  & qc_h2o_flux == 0 & u > 0.2 & AGC_c == 0 & qc_Sc == 0 & NEE_c > -20 & NEE_c < 15 & PAR > 800); % summer daytime
    end
    GEP_s = GEP(this_season); 
    ET_s = ET(this_season);
    WUE_s = GEP_s./ET_s;
    PAR_s = PAR(this_season);
    SWCs = SWC(this_season);
    VPD_s = vpd(this_season);
    gc_s = gc(this_season);
    %VPD_binedges = quantile(VPD_s,0:1/n:1); %* n bins of driver class, by quantile (n bins with same amount of data) 
    SWC_binedges = quantile(SWCs,0:1/3:1); % 3 quantiles of SWC
    c = 0.8; % dark blue/red is wet soil
    for j = 1:3 % SWC
        if s == 1
            colors = [c c 1];
        elseif s == 2
            colors = [1 c c];
        end
        c = c - 0.4;
        this_bin = find(SWCs >= SWC_binedges(j) & SWCs <= SWC_binedges(j+1)); %* iteration change SWC
        subplot(2,2,1); 
        binplot_v2(VPD_s(this_bin),GEP_s(this_bin),5,colors,6);
        subplot(2,2,2); 
        binplot_v2(VPD_s(this_bin),ET_s(this_bin),5,colors,6);
        subplot(2,2,3);
        binplot_v2(VPD_s(this_bin),WUE_s(this_bin),5,colors,6);
        subplot(2,2,4);
        binplot_v2(VPD_s(this_bin),gc_s(this_bin),5,colors,6);
    end
end

subplot(2,2,1); hold on; ylabel('GEP (\mumol m^-^2 s^-^1)'); box off;
ax = gca; ax.XTick = 0:1:5; ax.YTick = 4:4:16; ax.FontSize = 12;
set(gca,'xticklabel',[]); ylim([4 16]); xlim([0 5]);
subplot(2,2,2); hold on; ylabel('ET (mm)'); box off;
ax = gca; ax.XTick = 0:1:5; ax.YTick = 0.05:0.05:0.2; ax.FontSize = 12;
set(gca,'xticklabel',[]); ylim([0.05 0.2]); xlim([0 5]);
subplot(2,2,3); hold on; ylabel('WUE (\mumolm^-^2 mm^-^1)'); xlabel('VPD (kPa)'); box off;
ax = gca; ax.XTick = 0:1:5; ax.YTick = 50:25:150; ax.FontSize = 12;
ylim([50 150]); xlim([0 5]);
subplot(2,2,4); hold on; ylabel('gs (m s^-^1)'); xlabel('VPD (kPa)'); box off;
ax = gca; ax.XTick = 0:1:5; %ax.YTick = 50:50:150; %ylim([50 150]); 
ax.FontSize = 12; xlim([0 5]);
labelp = {'(a)','(b)','(c)','(d)'};
for i = 1:4
    subplot(2,2,i); hold on;
    text(0.90,0.98,labelp(i),'Units', 'Normalized', 'VerticalAlignment', 'Top');
    %sub_pos = get(gca,'position'); % get subplot axis position
    %set(gca,'position',sub_pos.*[1.2 1.2 1 1]) % stretch its width and height
end

subplot(2,2,4); hold on; 
h = legend('0-33%', '33-67%', '67-100%');
set(h,'FontSize',6);legend('boxoff');


% next figure


n = 4; % n drivers quantile/bin
LR_param = cell(2,3,n); LR_confint = cell(2,3,n);
driver_bin_median = nan(2,3,n);
Alpha = nan(2,3,n); Beta = nan(2,3,n); Rd = nan(2,3,n);
Alpha_sup = nan(2,3,n); Beta_sup = nan(2,3,n); Rd_sup = nan(2,3,n);
Alpha_inf = nan(2,3,n); Beta_inf = nan(2,3,n); Rd_inf = nan(2,3,n);
for s = 1:2 % 2 seasons, winter or summer
    if s == 1
        this_season = find(Month_t >= 5 & Month_t <= 8 & daytime == 1 & qc == 0 & u > 0.2 & AGC_c == 0 & qc_Sc == 0 & NEE_c > -20 & NEE_c < 15); % winter daytime
    elseif s == 2
        this_season = find((Month_t >= 10 | Month_t <= 2) & daytime == 1 & qc == 0 & u > 0.2 & AGC_c == 0 & qc_Sc == 0 & NEE_c > -20 & NEE_c < 15); % summer daytime
    end
    NEE_s = NEE_c(this_season); 
    PAR_s = PAR(this_season);
    for d = 1:3 % 3 drivers
        if d == 1 
            driver = tsoil(this_season);
        elseif d == 2
            driver = SWC(this_season);
        elseif d == 3
            driver = vpd(this_season);
        end
        driver_binedges = quantile(driver,0:1/n:1); %* n bins of driver class, by quantile (n bins with same amount of data) 
        for i=1:n %* n bins/quantile
            % driver_bin_middle(s,d,i) = (driver_binedges(i) + driver_binedges(i+1))/2; 
            this_driverbin=find(driver>=driver_binedges(i) & driver<=driver_binedges(i+1)); %* iteration change driver class (n bins/quantile)
            driver_bin_median(s,d,i) = median(driver(this_driverbin)); % x position for figure, makes more sense to be median of driver bin
            [LR_param{s,d,i},LR_confint{s,d,i}] = LRFIT_param_confint(PAR_s,NEE_s,this_driverbin);
            if isnan(LR_param{s,d,i}(1)) == 0 
                Alpha(s,d,i) = LR_param{s,d,i}(1); Beta(s,d,i) = LR_param{s,d,i}(2); Rd(s,d,i) = LR_param{s,d,i}(3);
                Alpha_sup(s,d,i) = LR_confint{s,d,i}(2,1); Beta_sup(s,d,i) = LR_confint{s,d,i}(2,2); Rd_sup(s,d,i) = LR_confint{s,d,i}(2,3);
                Alpha_inf(s,d,i) = LR_confint{s,d,i}(1,1); Beta_inf(s,d,i) = LR_confint{s,d,i}(1,2); Rd_inf(s,d,i) = LR_confint{s,d,i}(1,3);
            end
        end
    end
end

param_i = nan(1,n);
conf_sup_i = nan(1,n);
conf_inf_i = nan(1,n);
driver_bin_middle_i = nan(1,n);
figure;
for s = 1:2 % winter and summer
    if s == 1
        colors = 'b';
    elseif s == 2
        colors = 'r';
    end
    for p = 1:3 % alpha, beta, rd
        if p == 1 
            param = Alpha;
            conf_sup = Alpha_sup;
            conf_inf = Alpha_inf;
        elseif p == 2
            param = Beta;
            conf_sup = Beta_sup;
            conf_inf = Beta_inf;
        elseif p == 3
            param = Rd;
            conf_sup = Rd_sup;
            conf_inf = Rd_inf;
        end
        for d = 1:3 % tsoil, swc, vpd
            for i = 1:n
                param_i(1,i) = param(s,d,i);
                conf_sup_i(1,i) = conf_sup(s,d,i);
                conf_inf_i(1,i) = conf_inf(s,d,i);
                driver_bin_middle_i(1,i) = driver_bin_median(s,d,i);
            end
            subplot(3,3,(3*p-3)+d); hold on; % 1st subplot is Alpha, tsoil, 2nd is Alpha, SWC, ...
            plot(driver_bin_middle_i,param_i,'LineWidth',3,'Marker','o','MarkerEdgeColor', ...
    'none','MarkerFaceColor',colors,'Color',colors,'MarkerSize',10);
            x = driver_bin_middle_i;
            X = [x fliplr(x)];
            y1 = conf_inf_i; y2 = conf_sup_i;
            Y = [y1 fliplr(y2)];
            fill(X,Y,colors,'FaceAlpha',0.15,'EdgeAlpha',0.0); 
        end
    end
end

subplot(3,3,1); hold on; ylabel('\alpha (10^-^2 \mumol \mumol^-^1)');
ax = gca;  ax.YTick = 0:0.01:0.04; set(gca,'yticklabel',[0 1 2 3 4]);
subplot(3,3,2); hold on; set(gca,'yticklabel',[]);
subplot(3,3,3); hold on; set(gca,'yticklabel',[]);
subplot(3,3,4); hold on; ylabel('N_e_s (\mumol m^-^2 s^-^1)');
ax = gca;  ax.YTick = 0:5:20;
subplot(3,3,5); hold on; set(gca,'yticklabel',[]);
subplot(3,3,6); hold on; set(gca,'yticklabel',[]);
subplot(3,3,7); hold on; ylabel('R_d (\mumol m^-^2 s^-^1)'); xlabel('Tsoil_5_c_m (°C)');
subplot(3,3,8); hold on; xlabel('SWC_0_-_8_c_m (%)'); set(gca,'yticklabel',[]);
subplot(3,3,9); hold on; xlabel('VPD_3_0_m (kPa)');set(gca,'yticklabel',[]);
for i = 1:3
    subplot(3,3,i); hold on;
    ylim([0 0.04]);
    set(gca,'xticklabel',[]);
    set(gca,'FontSize',12);
end
for i = 4:6
    subplot(3,3,i); hold on;
    ylim([0 20]);
    set(gca,'xticklabel',[]);
    set(gca,'FontSize',12);
end
for i = 7:9
    subplot(3,3,i); hold on;
    ylim([0 8]);
    set(gca,'FontSize',12);
end
labelp = {'(a)','(b)','(c)','(d)','(e)','(f)','(g)','(h)','(i)'};
for i = 1:9
    subplot(3,3,i); hold on;
    text(0.90,0.98,labelp(i),'Units', 'Normalized', 'VerticalAlignment', 'Top');
    sub_pos = get(gca,'position'); % get subplot axis position
    set(gca,'position',sub_pos.*[1 1 1.2 1.2]) % stretch its width and height
end



% next figure



% NEE, ET and VPD diurnal trend, + PAR and tsoil, winter and summer
% Alexis, 02/12/2016

% Create half-hour vector (0-47), starting at half-hour = 1
i = 1;
j = 1;
n = length(NEE_c);
hHour_t = NaN(n,1);
while i < n + 1
    if j < 47
        hHour_t(i) = j;
        j = j + 1;
    else
        hHour_t(i) = 47;
        j = 0;
    end
    i = i + 1;
end

% Define color
light_blue = [0.8 0.8 1]; light_red = [1 0.8 0.8]; light_green = [0.8 1 0.8];

Diurn_quantiles_PAR_L = nan(47,3);
Diurn_quantiles_PAR_NL = nan(47,3);
Diurn_quantiles_Ta_L = nan(47,3);
Diurn_quantiles_Ta_NL = nan(47,3);
Diurn_quantiles_VPD_L = nan(47,3);
Diurn_quantiles_VPD_NL = nan(47,3);
Diurn_quantiles_NEE_NL = nan(47,3);
Diurn_quantiles_NEE_L = nan(47,3);
Diurn_quantiles_h2o_L = nan(47,3);
Diurn_quantiles_h2o_NL = nan(47,3);
for i = 2:48 % start 1st half-hour of the day, hHour_t = 1 corresponding to flux 00:00:00 to 00:30:00, end at hHour_t = 47, corresponding to 23:00:00 to 23:30:00, i.e. we don't look 23:30:30 to 00:00:00 
    these_L = find(hHour_t == i-1 & (month(t_s) >= 10 | month(t_s) <= 2)); % 1 Nov through 30 Apr
    these_L_qc = find(hHour_t==i-1 & (month(t_s) >= 10 | month(t_s) <= 2) & qc ~= 2 & qc_Sc ~= 2 & AGC_c == 0 & u > 0.2);
    these_NL = find(hHour_t==i-1 & month(t_s) >= 5 & month(t_s) <= 8); % 1 May through 31 Oct
    these_NL_qc = find(hHour_t==i-1 & month(t_s) >= 5 & month(t_s) <= 8 & qc ~= 2 & qc_Sc ~= 2 & AGC_c == 0 & u > 0.2);
    these_h2o_L_qc = find(hHour_t==i-1 & (month(t_s) >= 10 | month(t_s) <= 2) & qc_h2o_flux == 0 & AGC_c == 0);
    these_h2o_NL_qc = find(hHour_t==i-1 & month(t_s) >= 5 & month(t_s) <= 8 & qc_h2o_flux == 0 & AGC_c == 0);
    Diurn_quantiles_VPD_L(i-1,[1 2 3]) = quantile(vpd(these_L),[0.25, 0.5, 0.75]);
    Diurn_quantiles_VPD_NL(i-1,[1 2 3]) = quantile(vpd(these_NL),[0.25, 0.5, 0.75]);
    Diurn_quantiles_PAR_L(i-1,[1 2 3]) = quantile(PAR(these_L),[0.25, 0.5, 0.75]);
    Diurn_quantiles_PAR_NL(i-1,[1 2 3]) = quantile(PAR(these_NL),[0.25, 0.5, 0.75]);
    Diurn_quantiles_Ta_L(i-1,[1 2 3]) = quantile(t_air(these_L),[0.25, 0.5, 0.75]);
    Diurn_quantiles_Ta_NL(i-1,[1 2 3]) = quantile(t_air(these_NL),[0.25, 0.5, 0.75]);
    Diurn_quantiles_NEE_NL(i-1,[1 2 3]) = quantile(NEE_c(these_NL_qc),[0.25, 0.5, 0.75]);
    Diurn_quantiles_NEE_L(i-1,[1 2 3]) = quantile(NEE_c(these_L_qc),[0.25, 0.5, 0.75]);
    Diurn_quantiles_h2o_L(i-1,[1 2 3]) = quantile(h2o_flux(these_h2o_L_qc),[0.25, 0.5, 0.75]);
    Diurn_quantiles_h2o_NL(i-1,[1 2 3]) = quantile(h2o_flux(these_h2o_NL_qc),[0.25, 0.5, 0.75]);
end
NEE_ff_day = nan(length(NEE_c),1); % ff = full filter, daytime, 
for i = 1:length(NEE_c)
    if qc(i) == 0 && qc_Sc(i) == 0 && AGC_c(i) == 0 && u(i) > 0.15 && NEE_c(i) > -25 && NEE_c(i) < 10 && daytime(i) == 1
        NEE_ff_day(i) = NEE_c(i);
    end
end

x = 0.25:0.5:23.25;
figure; 
subplot(1,2,1); hold on;

x1=x;
x2=x;
y1 = zeros(1,47);
y2 = Diurn_quantiles_NEE_L(:,2)';
X = [x1 fliplr(x2)]; % !! x1,x2,y1 and y2 must be dim 1-n and 1-m (not n-1 and m-1)
Y = [y1 fliplr(y2)];
fill(X,Y,Y,'FaceAlpha',1,'EdgeAlpha',0.0);  caxis([-12,8.3267]); %colormap(map_NEE);


X = [x fliplr(x)];
y1 = Diurn_quantiles_NEE_L(:,1)'; y2 = Diurn_quantiles_NEE_L(:,3)';
y3 = Diurn_quantiles_h2o_L(:,1)'; y4 = Diurn_quantiles_h2o_L(:,3)';
Y = [y1 fliplr(y2)];
Y2 = [y3 fliplr(y4)];
% fill(X,Y,'r','FaceAlpha',0.4,'EdgeAlpha',0.0);
fill(X,Y,'k','FaceAlpha',0.2,'EdgeAlpha',0.0); % colormap(map_NEE); caxis([-12,8.3267]);
fill(X,Y2,'c','FaceAlpha',0.4,'EdgeAlpha',0.0);

x_p = [0 6 6 0];
y_p = [-15 -15 15 15];
patch(x_p,y_p,'b','FaceAlpha',.3,'EdgeColor','none');
x_p = [18 24 24 18];
y_p = [-15 -15 15 15];
patch(x_p,y_p,'b','FaceAlpha',.3,'EdgeColor','none');
plot(x,Diurn_quantiles_NEE_L(:,2),'k','LineWidth',2);
plot(x,Diurn_quantiles_h2o_L(:,2),'b','LineWidth',2);
plot([0 24],[0 0],'k'); axis([0 24 -15 15]); ax = gca; ax.YTick = -15:5:15;
set(ax,'FontSize',16); box off; ax.XTick = 0:2:24; %set(gca,'xticklabel',[]);
plot([12 12],[-15 15],'k-.');
xlabel ('Hour of the day'); ylabel('NEE or ET (\mu or mmolm^-^2s^-^1)');
text(0.90,0.98,'(a)','Units', 'Normalized', 'VerticalAlignment', 'Top');
text(0.07,0.57,'+ 1.0','Units', 'Normalized', 'VerticalAlignment', 'Top','Color','k');
text(0.30,0.48,'- 1.7','Units', 'Normalized', 'VerticalAlignment', 'Top','Color','k');
% text(0.57,0.27,'- 3','Units', 'Normalized', 'VerticalAlignment', 'Top','Color','b');
text(0.8,0.57,'+ 1.2','Units', 'Normalized', 'VerticalAlignment', 'Top','Color','k');
text(0.05,0.95,'+ 0.5 gCm^-^2d^-^1','Units', 'Normalized', 'VerticalAlignment', 'Top','Color','r','FontSize',12);
title('Summer (1^s^t Oct - 28^t^h Feb)');
ax1_pos = ax.Position; % position of first axes
ax2 = axes('Position',ax1_pos,...
'XAxisLocation','top',...
'YAxisLocation','right',...
'Color','none');
ax2.YColor = 'r';
line(x,Diurn_quantiles_VPD_L(:,2),'Parent',ax2,'Color','r','LineStyle','-','LineWidth',2);
line(x,Diurn_quantiles_VPD_L(:,1),'Parent',ax2,'Color',[1 0.8 0.8],'LineStyle','--');
line(x,Diurn_quantiles_VPD_L(:,3),'Parent',ax2,'Color',[1 0.8 0.8],'LineStyle','--');
axis([0 24 0 3]); set(ax2,'FontSize',16); ax2.YTick = 0:0.5:3;
set(gca,'xticklabel',[]); set(gca,'yticklabel',[]); % ylabel('VPD (kPa)');

subplot(1,2,2); hold on;
x1=x;
x2=x;
y1 = zeros(1,47);
y2 = Diurn_quantiles_NEE_NL(:,2)';
X = [x1 fliplr(x2)]; % !! x1,x2,y1 and y2 must be dim 1-n and 1-m (not n-1 and m-1)
Y = [y1 fliplr(y2)];
fill(X,Y,Y,'FaceAlpha',1,'EdgeAlpha',0.0); caxis([-12,8.3267]);

X = [x fliplr(x)];
y1 = Diurn_quantiles_NEE_NL(:,1)'; y2 = Diurn_quantiles_NEE_NL(:,3)';
y3 = Diurn_quantiles_h2o_NL(:,1)'; y4 = Diurn_quantiles_h2o_NL(:,3)';
Y = [y1 fliplr(y2)];
Y2 = [y3 fliplr(y4)];
% fill(X,Y,'r','FaceAlpha',0.4,'EdgeAlpha',0.0);
fill(X,Y,'k','FaceAlpha',0.2,'EdgeAlpha',0.0); % colormap(map_NEE); caxis([-12,8.3267]);
fill(X,Y2,'c','FaceAlpha',0.4,'EdgeAlpha',0.0);

x_p = [0 8 8 0];
y_p = [-15 -15 15 15];
patch(x_p,y_p,'b','FaceAlpha',.3,'EdgeColor','none');
x_p = [16 24 24 16];
y_p = [-15 -15 15 15];
patch(x_p,y_p,'b','FaceAlpha',.3,'EdgeColor','none');
plot(x,Diurn_quantiles_NEE_NL(:,2),'k','LineWidth',2);
plot(x,Diurn_quantiles_h2o_NL(:,2),'b','LineWidth',2);
plot([0 24],[0 0],'k'); axis([0 24 -15 15]); ax = gca; ax.YTick = -15:5:15;
set(ax,'FontSize',16); box off; ax.XTick = 0:2:24; %set(gca,'xticklabel',[]);
plot([12 12],[-15 15],'k-.');
set(gca,'yticklabel',[]);
xlabel ('Hour of the day'); % ylabel('NEE (\mumolm^-^2s^-^1)');
text(0.90,0.98,'(b)','Units', 'Normalized', 'VerticalAlignment', 'Top');
text(0.07,0.57,'+ 0.5','Units', 'Normalized', 'VerticalAlignment', 'Top');
text(0.37,0.48,'- 2.3','Units', 'Normalized', 'VerticalAlignment', 'Top');
text(0.8,0.57,'+ 0.6','Units', 'Normalized', 'VerticalAlignment', 'Top');
text(0.07,0.95,'- 1.2 gCm^-^2d^-^1','Units', 'Normalized', 'VerticalAlignment', 'Top','Color','b','FontSize',12);
title('Winter (1^s^t May - 31^t^h Aug)');
ax1_pos = ax.Position; % position of first axes
ax2 = axes('Position',ax1_pos,...
'XAxisLocation','top',...
'YAxisLocation','right',...
'Color','none');
ax2.YColor = 'r';
line(x,Diurn_quantiles_VPD_NL(:,2),'Parent',ax2,'Color','r','LineStyle','-','LineWidth',2);
line(x,Diurn_quantiles_VPD_NL(:,1),'Parent',ax2,'Color',[1 0.8 0.8],'LineStyle','--');
line(x,Diurn_quantiles_VPD_NL(:,3),'Parent',ax2,'Color',[1 0.8 0.8],'LineStyle','--');
axis([0 24 0 3]); set(ax2,'FontSize',16); ax2.YTick = 0:0.5:3;
set(gca,'xticklabel',[]); ylabel('VPD (kPa)');

% subplot(3,2,5); hold on;
% i = find((Month_t >= 10 | Month_t <= 2) & hour(t_s) >= 12 & daytime == 1);
% scatter(PAR(i),NEE_ff_day(i),15,'MarkerFaceColor',light_red,'MarkerEdgeColor','none'); 
% i = find((Month_t >= 10 | Month_t <= 2) & hour(t_s) < 12 & daytime == 1); 
% scatter(PAR(i),NEE_ff_day(i),15,'MarkerFaceColor',light_blue,'MarkerEdgeColor','none'); 

% axis([0 2500 -12 6]);
% ylabel('NEE (\mumol m^-^2 s^-^1)');
% xlabel('PAR (\mumol m^-^2 s^-^1)');
% set(gca,'FontSize',16);
% % PAR diurnal quantiles (left y-axis) and Ta quantiles (right y-axis),
% % summer
% X = [x fliplr(x)];
% y1 = Diurn_quantiles_PAR_L(:,1)'; y2 = Diurn_quantiles_PAR_L(:,3)';
% Y = [y1 fliplr(y2)];
% % fill(X,Y,'r','FaceAlpha',0.4,'EdgeAlpha',0.0);
% fill(X,Y,'k','FaceAlpha',0.2,'EdgeAlpha',0.0); % colormap(map_NEE); caxis([-12,8.3267]);
% 
% x_p = [0 6 6 0];
% y_p = [0 0 2500 2500];
% patch(x_p,y_p,'b','FaceAlpha',.3,'EdgeColor','none');
% x_p = [18 24 24 18];
% y_p = [0 0 2500 2500];
% patch(x_p,y_p,'b','FaceAlpha',.3,'EdgeColor','none');
% % X = [x fliplr(x)];
% % y1 = Diurn_quantiles_PAR_L(:,1)'; y2 = Diurn_quantiles_PAR_L(:,3)';
% % Y = [y1 fliplr(y2)];
% % fill(X,Y,'y','FaceAlpha',1,'EdgeAlpha',0.0);
% plot(x,Diurn_quantiles_PAR_L(:,2),'k','LineWidth',2);
% %plot(x,Diurn_quantiles_PAR_L(:,1),'k:','LineWidth',1);
% %plot(x,Diurn_quantiles_PAR_L(:,3),'k:','LineWidth',1);
% axis([0 24 0 2500]); ax = gca; ax.XTick = 0:2:24; ylabel('PAR (\mumolm^-^2s^-^1)');
% set(ax,'FontSize',12);
% plot([12 12],[0 2500],'k-.');
% text(0.90,0.98,'(c)','Units', 'Normalized', 'VerticalAlignment', 'Top');
% ax1_pos = ax.Position; % position of first axes
% ax2 = axes('Position',ax1_pos,...
% 'XAxisLocation','top',...
% 'YAxisLocation','right',...
% 'Color','none');
% ax2.YColor = 'r';
% line(x,Diurn_quantiles_Ta_L(:,2),'Parent',ax2,'Color','r','LineStyle','-');
% line(x,Diurn_quantiles_Ta_L(:,1),'Parent',ax2,'Color','r','LineStyle',':');
% line(x,Diurn_quantiles_Ta_L(:,3),'Parent',ax2,'Color','r','LineStyle',':');
% axis([0 24 5 35]); set(ax2,'FontSize',12); %ax2.YTick = 0:0.5:3;
% set(gca,'xticklabel',[]); set(gca,'yticklabel',[]); %ylabel('Ta (C)');
% 
% subplot(3,2,6); hold on;
% i = find(Month_t >= 5 & Month_t <= 8 & hour(t_s) >= 12 & daytime == 1);
% scatter(PAR(i),NEE_ff_day(i),15,'MarkerFaceColor',light_red,'MarkerEdgeColor','none'); 
% i = find(Month_t >= 5 & Month_t <= 8 & hour(t_s) < 12 & daytime == 1);
% scatter(PAR(i),NEE_ff_day(i),15,'MarkerFaceColor',light_blue,'MarkerEdgeColor','none');

% axis([0 2500 -12 6]);
% set(gca,'yticklabel',[]);
% xlabel('PAR (\mumol m^-^2 s^-^1)');
% legend('Afternoon','Morning');
% set(gca,'FontSize',16);
% % PAR diurnal quantiles (left y-axis) and Ta quantiles (right y-axis), winter
% X = [x fliplr(x)];
% y1 = Diurn_quantiles_PAR_NL(:,1)'; y2 = Diurn_quantiles_PAR_NL(:,3)';
% Y = [y1 fliplr(y2)];
% % fill(X,Y,'r','FaceAlpha',0.4,'EdgeAlpha',0.0);
% fill(X,Y,'k','FaceAlpha',0.2,'EdgeAlpha',0.0); % colormap(map_NEE); caxis([-12,8.3267]);
% 
% x_p = [0 8 8 0];
% y_p = [0 0 2500 2500];
% patch(x_p,y_p,'b','FaceAlpha',.3,'EdgeColor','none');
% x_p = [16 24 24 16];
% y_p = [0 0 2500 2500];
% patch(x_p,y_p,'b','FaceAlpha',.3,'EdgeColor','none');
% % X = [x fliplr(x)];
% % y1 = Diurn_quantiles_PAR_NL(:,1)'; y2 = Diurn_quantiles_PAR_NL(:,3)';
% % Y = [y1 fliplr(y2)];
% % fill(X,Y,'y','FaceAlpha',1,'EdgeAlpha',0.0);
% plot(x,Diurn_quantiles_PAR_NL(:,2),'k','LineWidth',2);
% % plot(x,Diurn_quantiles_PAR_NL(:,1),'k:','LineWidth',1);
% % plot(x,Diurn_quantiles_PAR_NL(:,3),'k:','LineWidth',1);
% axis([0 24 0 2500]); ax = gca; ax.XTick = 0:2:24; set(gca,'yticklabel',[]);
% plot([12 12],[0 2500],'k-.');
% set(ax,'FontSize',12);
% text(0.90,0.98,'(d)','Units', 'Normalized', 'VerticalAlignment', 'Top');
% ax1_pos = ax.Position; % position of first axes
% ax2 = axes('Position',ax1_pos,...
% 'XAxisLocation','top',...
% 'YAxisLocation','right',...
% 'Color','none');
% ax2.YColor = 'r';
% line(x,Diurn_quantiles_Ta_NL(:,2),'Parent',ax2,'Color','r','LineStyle','-');
% line(x,Diurn_quantiles_Ta_NL(:,1),'Parent',ax2,'Color','r','LineStyle',':');
% line(x,Diurn_quantiles_Ta_NL(:,3),'Parent',ax2,'Color','r','LineStyle',':');
% axis([0 24 5 35]); set(ax2,'FontSize',12); %ax2.YTick = 0:0.5:3;
% set(gca,'xticklabel',[]); ylabel('Ta (°C)');

% for i = 1:4
%     subplot(3,3,i); hold on;
%     sub_pos = get(gca,'position'); % get subplot axis position
%     set(gca,'position',sub_pos.*[1 1 1.2 1.2]) % stretch its width and height
% end


% subplot(3,2,5); hold on;
% these1 = find((month(t_s) >= 10 | month(t_s) <= 2) & qc == 0 & qc_Sc == 0 & AGC_c == 0 & u > 0.2 & daytime == 1 & vpd < 1.5);
% these2 = find((month(t_s) >= 10 | month(t_s) <= 2) & qc == 0 & qc_Sc == 0 & AGC_c == 0 & u > 0.2 & daytime == 1 & vpd >= 1.5 & vpd < 3);
% these3 = find((month(t_s) >= 10 | month(t_s) <= 2) & qc == 0 & qc_Sc == 0 & AGC_c == 0 & u > 0.2 & daytime == 1 & vpd >= 3);
% scatter(PAR(these1),NEE_c(these1),10,'b','filled'); %caxis([nanmin(vpd) nanmax(vpd)]); colormap(customap);
% scatter(PAR(these2),NEE_c(these2),10,'g','filled');
% scatter(PAR(these3),NEE_c(these3),10,'r','filled');
% axis([0 2500 -20 10]); ax = gca; set(ax,'FontSize',12); box off;
% % set(gca,'Position',[0 0 1 1]);
% ylabel('NEE (\mumolm^-^2s^-^1)');
% xlabel('PAR (\mumolm^-^2s^-^1)');
% hold on; plot([0 2500],[0 0],'k');
% text(0.90,0.98,'(e)','Units', 'Normalized', 'VerticalAlignment', 'Top');
% % legend('VPD < 1.5','VPD 1.5-3','VPD > 3');
% 
% subplot(3,2,6); hold on;
% these1 = find(month(t_s) >= 5 & month(t_s) <= 8 & qc == 0 & qc_Sc == 0 & AGC_c == 0 & u > 0.2 & daytime == 1 & vpd < 1.5);
% these2 = find(month(t_s) >= 5 & month(t_s) <= 8 & qc == 0 & qc_Sc == 0 & AGC_c == 0 & u > 0.2 & daytime == 1 & vpd >= 1.5);
% %these3 = find(month(t_s) >= 5 & month(t_s) <= 8 & qc == 0 & qc_Sc == 0 & AGC_c == 0 & u > 0.2 & daytime == 1);
% scatter(PAR(these1),NEE_c(these1),10,'b','filled'); %colormap(customap); caxis([nanmin(vpd) nanmax(vpd)]);
% scatter(PAR(these2),NEE_c(these2),10,'g','filled');
% axis([0 2500 -20 10]); ax = gca; set(ax,'FontSize',12); box off; set(gca,'yticklabel',[]);
% % set(gca,'Position',[0 0 1 1]);
% xlabel('PAR (\mumolm^-^2s^-^1)');
% hold on; plot([0 2500],[0 0],'k');
% %caxis([nanmin(vpd) nanmax(vpd)]); %colorbar; h = colorbar; ylabel(h, 'VPD (kPa)');
% text(0.90,0.98,'(f)','Units', 'Normalized', 'VerticalAlignment', 'Top');

% h = colorbar;
% ylabel(h, 'VPD (kPa)');


% NEE vs. PAR, VPD color, summer and winter panels

lowVPD = [0.7569    0.8667    0.7765];
mediumVPD = [0.2314    0.4431    0.3373];
highVPD = [0.0706    0.2118    0.1412];
morningc = [0.6 0.6 1];
afternoonc = [1 0.8 0.8];


figure
subplot(1,2,1); hold on;
these1 = find((month(DateTime_CUP) >= 10 | month(DateTime_CUP) <= 2) & qc == 0 & qc_Sc == 0 & AGC_c == 0 & u > 0.2 & daytime == 1 & vpd < 1.5);
these2 = find((month(DateTime_CUP) >= 10 | month(DateTime_CUP) <= 2) & qc == 0 & qc_Sc == 0 & AGC_c == 0 & u > 0.2 & daytime == 1 & vpd >= 1.5 & vpd < 3);
these3 = find((month(DateTime_CUP) >= 10 | month(DateTime_CUP) <= 2) & qc == 0 & qc_Sc == 0 & AGC_c == 0 & u > 0.2 & daytime == 1 & vpd >= 3);
scatter(PAR(these1),NEE_c(these1),20,lowVPD,'filled'); %caxis([nanmin(vpd) nanmax(vpd)]); colormap(customap);
scatter(PAR(these2),NEE_c(these2),20,mediumVPD,'filled');
scatter(PAR(these3),NEE_c(these3),20,highVPD,'filled');
i = find((Month_t >= 10 | Month_t <= 2) & hour(t_s) < 12 & daytime == 1);
binplot(PAR(i),NEE_ff_day(i),15,morningc,10);
i = find((Month_t >= 10 | Month_t <= 2) & hour(t_s) >= 12 & daytime == 1);
binplot(PAR(i),NEE_ff_day(i),15,afternoonc,10);
plot([0 2500],[0 0],'k');
% binplot(PAR(these1),NEE_c(these1),4,dark_lowVPD,10);
% binplot(PAR(these2),NEE_c(these2),4,dark_mediumVPD,10);
% binplot(PAR(these3),NEE_c(these3),4,dark_highVPD,10);
axis([0 2500 -20 10]); ax = gca; set(ax,'FontSize',16); box off;
ax.XTick = 0:500:2500;
% set(gca,'Position',[0 0 1 1]);
ylabel('NEE (\mumol m^-^2 s^-^1)');
xlabel('PAR (\mumol m^-^2 s^-^1)');
%title('Summer (1^s^t Oct - 28^t^h Feb)');
text(0.90,0.98,'(a)','Units', 'Normalized', 'VerticalAlignment', 'Top');


subplot(1,2,2); hold on;
these1 = find(month(DateTime_CUP) >= 5 & month(DateTime_CUP) <= 8 & qc == 0 & qc_Sc == 0 & AGC_c == 0 & u > 0.2 & daytime == 1 & vpd < 1.5);
these2 = find(month(DateTime_CUP) >= 5 & month(DateTime_CUP) <= 8 & qc == 0 & qc_Sc == 0 & AGC_c == 0 & u > 0.2 & daytime == 1 & vpd >= 1.5);
%these3 = find(month(t_s) >= 5 & month(t_s) <= 8 & qc == 0 & qc_Sc == 0 & AGC_c == 0 & u > 0.2 & daytime == 1);
dotsa = scatter(PAR(these1),NEE_c(these1),20,lowVPD,'filled'); %colormap(customap); caxis([nanmin(vpd) nanmax(vpd)]);
dotsb = scatter(PAR(these2),NEE_c(these2),20,mediumVPD,'filled');
dotsc = scatter(3000,2,20,highVPD,'filled'); % for the legend
dotsd = plot(3000,2,'LineWidth',3,'Marker','o','MarkerEdgeColor', ...
    'none','MarkerFaceColor',morningc,'Color',morningc,'MarkerSize',10);
dotse = plot(3000,2,'LineWidth',3,'Marker','o','MarkerEdgeColor', ...
    'none','MarkerFaceColor',afternoonc,'Color',afternoonc,'MarkerSize',10);
i = find(Month_t >= 5 & Month_t <= 8 & hour(t_s) < 12 & daytime == 1);
binplot(PAR(i),NEE_ff_day(i),15,morningc,10);
i = find(Month_t >= 5 & Month_t <= 8 & hour(t_s) >= 12 & daytime == 1);
binplot(PAR(i),NEE_ff_day(i),15,afternoonc,10);
plot([0 2500],[0 0],'k');
% binplot(PAR(these1),NEE_c(these1),4,dark_lowVPD,10);
% binplot(PAR(these2),NEE_c(these2),4,dark_mediumVPD,10);
axis([0 2500 -20 10]); ax = gca; set(ax,'FontSize',16); box off;
set(gca,'yticklabel',[]); ax.XTick = 0:500:2500;
% set(gca,'Position',[0 0 1 1]);
xlabel('PAR (\mumolm^-^2s^-^1)');
%caxis([nanmin(vpd) nanmax(vpd)]); %colorbar; h = colorbar; ylabel(h, 'VPD (kPa)');
%title('Winter (1^s^t May - 31^t^h Aug)');
text(0.90,0.98,'(b)','Units', 'Normalized', 'VerticalAlignment', 'Top');
L = legend([dotsa dotsb dotsc dotsd dotse],'VPD < 1.5','VPD [1.5 3]','VPD > 3','Morning','Afternoon');
L.Location = 'southeast';
L.FontSize = 8;


% next figure


% Monthly "PC" as median GEP when PAR ~ 1000 and VPD ~ 1
% !! contains all gap-filled GEP. 
% Monthly surface conductance, "Gs_m", as median Gs when PAR ~ 1000 and VPD ~ 1
% !! contains all gap-filled gs, no rain filter. 
PC = nan(12*3,3);
gs_m = nan(12*3,3);
date_m = datetime;
LAI_m = nan(36,1);
start_d = datenum(DateTime_LAI(1));
end_d = datenum(DateTime_LAI(297));
st_LAI = find(datenum(DateTime_CUP) == start_d);
en_LAI = find(datenum(DateTime_CUP) == end_d);
hh_num = datenum(DateTime_CUP(2)) - datenum(DateTime_CUP(1));
xq = start_d:hh_num:end_d;
% smooth LAI
LAIdatenum = datenum(DateTime_LAI);
[xData, yData] = prepareCurveData( LAIdatenum, LAI );
% Set up fittype and options.
ft = fittype( 'smoothingspline' );
excludedPoints = excludedata( xData, yData, 'Indices', [8 129 131 132] );
opts = fitoptions( 'Method', 'SmoothingSpline' );
opts.SmoothingParam = 4.46646226829198e-05;
opts.Exclude = excludedPoints;
% Fit model to data.
[fitresult, gof] = fit( xData, yData, ft, opts );
s_LAI = feval(fitresult,xq);
i = 1;
for y = 1:3
    for m = 1:12
        LAI_m(i) = mean(s_LAI(Year_t(st_LAI:en_LAI) == 2013+y & Month_t(st_LAI:en_LAI) == m & Day_t(st_LAI:en_LAI) == 15));
        use = find(Year_t == 2013+y & Month_t == m & PAR > 800 & PAR < 1200 & vpd > 1 & vpd < 1.5);
        PC(i,1) = quantile(GEP(use),0.25);
        PC(i,2) = quantile(GEP(use),0.5);
        PC(i,3) = quantile(GEP(use),0.75);
        gs_m(i,1) = quantile(gc(use),0.25);
        gs_m(i,2) = quantile(gc(use),0.5);
        gs_m(i,3) = quantile(gc(use),0.75);
        date_m(i) = datetime(2013+y,m,15);
        i = i+1;
    end
end

figure; subplot(4,1,1); plot(DateTime_LAI,LAI,'LineStyle','none','Marker','o', ...
    'MarkerFaceColor',[0.8 0.8 0.8],'MarkerEdgeColor','none','MarkerSize',2); 
hold on; plot(xq,s_LAI,'LineWidth',2,'Color','k');
subplot(4,1,2); plot(date_m,PC(:,2),'LineWidth',2,'Color','k');
hold on; plot(date_m,PC(:,1),'LineWidth',1,'Color',[0.8 0.8 0.8]);
plot(date_m,PC(:,3),'LineWidth',1,'Color',[0.8 0.8 0.8]);
subplot(4,1,3); plot(date_m,gs_m(:,2),'LineWidth',2,'Color','k');
hold on; plot(date_m,gs_m(:,1),'LineWidth',1,'Color',[0.8 0.8 0.8]);
plot(date_m,gs_m(:,3),'LineWidth',1,'Color',[0.8 0.8 0.8]);
subplot(4,1,1); ylim([0.6 1.1]);
xlim([datenum(datetime(2014,01,01)) datenum(datetime(2017,01,01))]);
set(gca,'FontSize',12); ylabel('LAI (m^2 m^-^2)');
subplot(4,1,2); ylim([5 20]);
xlim([datenum(datetime(2014,01,01)) datenum(datetime(2017,01,01))]);
set(gca,'FontSize',12); ylabel('PC (\mumol m^-^2 s^-^1)');
subplot(4,1,3); ylim([0.002 0.008]);
xlim([datenum(datetime(2014,01,01)) datenum(datetime(2017,01,01))]);
set(gca,'FontSize',12); ylabel('Gs_o_p_t (m s^-^1)');
subplot(4,1,4); plot(DateTime_CUP,SWC,'LineWidth',2,'Color','k');
xlim([datenum(datetime(2014,01,01)) datenum(datetime(2017,01,01))]);
set(gca,'FontSize',12); ylabel('SWC _0_-_8_c_m');

mdl = fitlm(LAI_m,PC(:,2));
figure; plot(mdl,'Color','k');
xlabel('Leaf area index'); ylabel('Photosynthetic capacity');
box off; title(''); h = legend; set(h,'FontSize',6);
set(h,'Location','northwest');