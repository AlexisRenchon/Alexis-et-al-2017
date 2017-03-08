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
year_t = ncread(source_processed,'Year'); year_t = reshape(year_t,[],1);
month_t = ncread(source_processed,'Month'); month_t = reshape(month_t,[],1);
day_t = ncread(source_processed,'Day'); day_t = reshape(day_t,[],1);
hour_t = ncread(source_processed,'Hour'); hour_t = reshape(hour_t,[],1);
minute_t = ncread(source_processed,'Minute'); minute_t = reshape(minute_t,[],1);
second_t = ncread(source_processed,'Second'); second_t = reshape(second_t,[],1);
DateTime_CUP = datetime(year_t,month_t,day_t,hour_t,minute_t,second_t);
SWC_s = ncread(source_processed,'Sws'); SWC_s = reshape(SWC_s,[],1); SWC_s = SWC_s*100;
VPD_30 = ncread(source_processed,'VPD'); VPD_30 = reshape(VPD_30,[],1); 
tsoil = ncread(source_processed,'Ts'); tsoil = reshape(tsoil,[],1); 
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

% PAR saturated ET and NEE vs VPD, for dry and wet soil (50% quantile), in summer
these_wet_NEE_s = find(qc_h2o_flux == 0 & qc == 0 & u > 0.15 & qc_Sc == 0 & AGC_c == 0 & (month(t_s) >= 10 | month(t_s) <= 2) ...
    & daytime == 1 & PAR > 1000 & SWC_s > nanmedian(SWC_s));
these_wet_NEE_w = find(qc_h2o_flux == 0 & qc == 0 & u > 0.15 & qc_Sc == 0 & AGC_c == 0 & month(t_s) >= 5 & month(t_s) <= 8 ...
    & daytime == 1 & PAR > 1000 & SWC_s > nanmedian(SWC_s));
these_wet_h2o_s = find(qc_h2o_flux == 0 & AGC_c == 0 & (month(t_s) >= 10 | month(t_s) <= 2) ...
    & daytime == 1 & PAR > 1000 & SWC_s > nanmedian(SWC_s));
these_wet_h2o_w = find(qc_h2o_flux == 0 & AGC_c == 0 & month(t_s) >= 5 & month(t_s) <= 8 ...
    & daytime == 1 & PAR > 1000 & SWC_s > nanmedian(SWC_s));
these_dry_NEE_s = find(qc_h2o_flux == 0 & qc == 0 & u > 0.15 & qc_Sc == 0 & AGC_c == 0 & (month(t_s) >= 10 | month(t_s) <= 2) ...
    & daytime == 1 & PAR > 1000 & SWC_s < nanmedian(SWC_s));
these_dry_NEE_w = find(qc_h2o_flux == 0 & qc == 0 & u > 0.15 & qc_Sc == 0 & AGC_c == 0 & month(t_s) >= 5 & month(t_s) <= 8 ...
    & daytime == 1 & PAR > 1000 & SWC_s < nanmedian(SWC_s));
these_dry_h2o_s = find(qc_h2o_flux == 0 & AGC_c == 0 & (month(t_s) >= 10 | month(t_s) <= 2) ...
    & daytime == 1 & PAR > 1000 & SWC_s < nanmedian(SWC_s));
these_dry_h2o_w = find(qc_h2o_flux == 0 & AGC_c == 0 & month(t_s) >= 5 & month(t_s) <= 8 ...
    & daytime == 1 & PAR > 1000 & SWC_s < nanmedian(SWC_s));
WUE_NEE_ET = nan(length(t_s),1);
for i = 1:length(t_s)
if qc(i) == 0 && qc_Sc(i) == 0 && daytime(i) == 1 && AGC_c(i) == 0 && u(i) > 0.15 && qc_h2o_flux(i) == 0
WUE_NEE_ET(i) = NEE_c(i)/h2o_flux(i);
end
end
WUE_NEE_ET = -WUE_NEE_ET;
these_dry_s = find(PAR > 1000 & SWC_s < nanmedian(SWC_s) & (month(t_s) >= 10 | month(t_s) <= 2));
these_dry_w = find(PAR > 1000 & SWC_s < nanmedian(SWC_s) & month(t_s) >= 5 & month(t_s) <= 8);
these_wet_s = find(PAR > 1000 & SWC_s > nanmedian(SWC_s) & (month(t_s) >= 10 | month(t_s) <= 2));
these_wet_w = find(PAR > 1000 & SWC_s > nanmedian(SWC_s) & month(t_s) >= 5 & month(t_s) <= 8);

figure; subplot(3,2,1); hold on;
binplot(VPD_30(these_wet_NEE_s),-NEE_c(these_wet_NEE_s),5,'b');
binplot(VPD_30(these_dry_NEE_s),-NEE_c(these_dry_NEE_s),5,'r');
ylabel('NEP (\mumol m^-^2 s^-^1)'); title('Summer (1^s^t Oct - 28^t^h Feb)');
axis([0 5 -2 12]); ax = gca; set(ax,'FontSize',16); set(gca,'xticklabel',[]);
subplot(3,2,2); hold on;
binplot(VPD_30(these_wet_NEE_w),-NEE_c(these_wet_NEE_w),5,'b');
binplot(VPD_30(these_dry_NEE_w),-NEE_c(these_dry_NEE_w),5,'r');
title('Winter (1^s^t May - 31^t^h Aug)'); ax = gca; set(ax,'FontSize',16);
axis([0 5 -2 12]); set(gca,'yticklabel',[]); set(gca,'xticklabel',[]);
subplot(3,2,3);
binplot(VPD_30(these_dry_h2o_s),h2o_flux(these_dry_h2o_s),5,'r');
binplot(VPD_30(these_wet_h2o_s),h2o_flux(these_wet_h2o_s),5,'b');
box off;
ylabel('ET (mmol m^-^2 s^-^1)');
axis([0 5 1 8]); ax = gca; set(ax,'FontSize',16); set(gca,'xticklabel',[]);
subplot(3,2,4);
binplot(VPD_30(these_dry_h2o_w),h2o_flux(these_dry_h2o_w),5,'r');
binplot(VPD_30(these_wet_h2o_w),h2o_flux(these_wet_h2o_w),5,'b');
box off; axis([0 5 1 8]); ax = gca; set(ax,'FontSize',16); set(gca,'yticklabel',[]);
set(gca,'xticklabel',[]);
subplot(3,2,5); hold on; 
binplot(VPD_30(these_wet_s),WUE_NEE_ET(these_wet_s),5,'b');
binplot(VPD_30(these_dry_s),WUE_NEE_ET(these_dry_s),5,'r');
ylabel('WUE (\mumol mmol^-^1)'); axis([0 5 -1 8]); xlabel('VPD (kPa)');
ax = gca; set(ax,'FontSize',16);
subplot(3,2,6); hold on; 
binplot(VPD_30(these_wet_w),WUE_NEE_ET(these_wet_w),5,'b');
binplot(VPD_30(these_dry_w),WUE_NEE_ET(these_dry_w),5,'r');
axis([0 5 -1 8]); xlabel('VPD (kPa)'); ax = gca; set(ax,'FontSize',16);
set(gca,'yticklabel',[]);

labelp = {'(a)','(b)','(c)','(d)','(e)','(f)','(g)','(h)','(i)'};
for i = 1:6
    subplot(3,2,i); hold on;
    text(0.90,0.98,labelp(i),'Units', 'Normalized', 'VerticalAlignment', 'Top');
    sub_pos = get(gca,'position'); % get subplot axis position
    set(gca,'position',sub_pos.*[1 1 1.17 1.17]) % stretch its width and height
end


