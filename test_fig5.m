% Clear the workspace
clear;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Load input data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Inputs:
% daytime (half-hourly 0 for night, 1 for day)
% qc (half-hourly qc Mauder and Foken 2004)
% qc_h2o_flux (half-hourly qc for ET)
% h2o_flux (half-hourly ET)
% AGC_c (half-hourly, 0 = maximum signal strength quality value, else NaN)
% u (half-hourly u* at 30m)
% qc_Sc (half-hourly qc for storage flux, 0 is best quality)
% NEE_c (CO2 turbulent flux (Fc) + CO2 storage flux (Sc))
% PAR (half-hourly PAR at 30m)
% month_t (half-hourly month vector)
% SWC_s (half-hourly soil water content 0-8cm)
% VPD_30 (half-hourly VPD at 30m)
% tsoil
% tair

% LRFIT_param_confint() function

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
SWC_s = ncread(source_processed,'Sws'); SWC_s = reshape(SWC_s,[],1); SWC_s = SWC_s*100;
VPD_30 = ncread(source_processed,'VPD'); VPD_30 = reshape(VPD_30,[],1); 
tsoil = ncread(source_processed,'Ts'); tsoil = reshape(tsoil,[],1); 
GEP = ncread(source_processed,'GPP_SOLO_EP'); GEP = reshape(GEP,[],1); 
ET = ncread(source_processed,'ET'); ET = reshape(ET,[],1);
Precip  = ncread(source_processed,'Precip'); Precip = reshape(Precip,[],1);
% clear unused variables
clearvars source_processed finfo;

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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Figure script
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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

% PAR saturated (> 800) ET, GEP and WUE vs VPD, for SWC quantiles (n quantiles)
figure; hold on;
for s = 1:2 % 2 seasons, winter or summer
    if s == 1
        %this_season = find(Precip_2daysago < 1 & Precip_1daysago < 0.5 & Precip_halfdayago < 0.2 & Month_t >= 5 & Month_t <= 8 & daytime == 1 & qc == 0 & qc_h20_flux == 0 & u > 0.2 & AGC_c == 0 & qc_Sc == 0 & NEE_c > -20 & NEE_c < 15 & PAR > 1000); % winter daytime
        this_season = find(Month_t >= 5 & Month_t <= 8 & daytime == 1 & qc == 0 & qc_h2o_flux == 0 & u > 0.2 & AGC_c == 0 & qc_Sc == 0 & NEE_c > -20 & NEE_c < 15 & PAR > 800); % winter daytime
    elseif s == 2
        %this_season = find(Precip_2daysago < 1 & Precip_1daysago < 0.5 & Precip_halfdayago < 0.2 & (Month_t >= 10 | Month_t <= 2) & daytime == 1 & qc == 0  & qc_h20_flux == 0 & u > 0.2 & AGC_c == 0 & qc_Sc == 0 & NEE_c > -20 & NEE_c < 15 & PAR > 1000); % summer daytime
        this_season = find((Month_t >= 10 | Month_t <= 2) & daytime == 1 & qc == 0  & qc_h2o_flux == 0 & u > 0.2 & AGC_c == 0 & qc_Sc == 0 & NEE_c > -20 & NEE_c < 15 & PAR > 800); % summer daytime
    end
    GEP_s = GEP(this_season); 
    ET_s = ET(this_season);
    WUE_s = GEP_s./ET_s;
    PAR_s = PAR(this_season);
    SWC_ss = SWC_s(this_season);
    VPD_s = VPD_30(this_season);
    %VPD_binedges = quantile(VPD_s,0:1/n:1); %* n bins of driver class, by quantile (n bins with same amount of data) 
    SWC_binedges = quantile(SWC_ss,0:1/3:1); % 3 quantiles of SWC
    c = 0.8; % dark blue/red is wet soil
    for j = 1:3 % SWC
        if s == 1
            colors = [c c 1];
        elseif s == 2
            colors = [1 c c];
        end
        c = c - 0.4;
        this_bin = find(SWC_ss >= SWC_binedges(j) & SWC_ss <= SWC_binedges(j+1)); %* iteration change SWC
        subplot(3,1,1); 
        binplot_v2(VPD_s(this_bin),GEP_s(this_bin),6,colors,6);
        subplot(3,1,2); 
        binplot_v2(VPD_s(this_bin),ET_s(this_bin),6,colors,6);
        subplot(3,1,3);
        binplot_v2(VPD_s(this_bin),WUE_s(this_bin),6,colors,6);
    end
end

subplot(3,1,1); hold on; ylabel('GEP (\mumol m^-^2 s^-^1)'); box off;
ax = gca; ax.XTick = 0:1:4; ax.YTick = 4:2:16; ax.FontSize = 12; set(gca,'xticklabel',[]); ylim([4 16]);
subplot(3,1,2); hold on; ylabel('ET (mm)'); box off;
ax = gca; ax.XTick = 0:1:4; ax.YTick = 0.05:0.05:0.2; ax.FontSize = 12; set(gca,'xticklabel',[]); ylim([0.05 0.2]);
subplot(3,1,3); hold on; ylabel('WUE (\mumolm^-^2 mm^-^1)'); xlabel('VPD (kPa)'); box off;
ax = gca; ax.XTick = 0:1:4; ax.YTick = 50:50:150; ax.FontSize = 12; ylim([50 150]);
labelp = {'(a)','(b)','(c)'};
for i = 1:3
    subplot(3,1,i); hold on;
    text(0.90,0.98,labelp(i),'Units', 'Normalized', 'VerticalAlignment', 'Top');
    %sub_pos = get(gca,'position'); % get subplot axis position
    %set(gca,'position',sub_pos.*[1.2 1.2 1 1]) % stretch its width and height
end
