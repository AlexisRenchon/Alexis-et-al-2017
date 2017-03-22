% Clear the workspace
clear;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Load input data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Inputs to load
% Daily_ER daily_NEE daily_GPP daily_date SWC_s VPD Ta Ts DateTime_CUP
% Precip_daily DateTime_LAI DateTime_PAR LAI PAR

% Load processed data (netCDF file)
% Location of processed input file 
source_processed = 'Input\CumberlandPlain_2014_L6_EP_moderate.nc';
% Information on file, including variable name
finfo = ncinfo(source_processed);
% Read variable from .cd file
Yeari = ncread(source_processed,'Year'); Yeari = reshape(Yeari,[],1);
Monthi = ncread(source_processed,'Month'); Monthi = reshape(Monthi,[],1);
Dayi = ncread(source_processed,'Day'); Dayi = reshape(Dayi,[],1);
Houri = ncread(source_processed,'Hour'); Houri = reshape(Houri,[],1);
Minutei = ncread(source_processed,'Minute'); Minutei = reshape(Minutei,[],1);
Secondi = ncread(source_processed,'Second'); Secondi = reshape(Secondi,[],1);
DateTime_CUP = datetime(Yeari,Monthi,Dayi,Houri,Minutei,Secondi);
Ts = ncread(source_processed,'Ts'); Ts = reshape(Ts,[],1); 
Ta = ncread(source_processed,'Ta'); Ta = reshape(Ta,[],1); 
SWC_s = ncread(source_processed,'Sws'); SWC_s = reshape(SWC_s,[],1); SWC_s = SWC_s*100;
VPD = ncread(source_processed,'VPD'); VPD = reshape(VPD,[],1); 
% clear unused variables
clearvars source_processed finfo;


% Load raw data (.csv file)
% Location of raw input file
source_raw = 'Input\CUP_EddyPro_QC_170207.csv';
% Create datastore to access collection of data
ds_raw = datastore(source_raw);
% Select variable of interest
ds_raw.SelectedVariableNames = {'DateTime','NEE_c','qc_Sc','AGC_c','daytime','qc_co2_flux','u_'}; % for example
% Read selected variables, save it in the workspace as a table
Raw_Data = readall(ds_raw);
% Get data from the table, change format if necessary
NEE_c = Raw_Data.NEE_c;
qc = Raw_Data.qc_co2_flux;
qc_Sc = Raw_Data.qc_Sc;
AGC_c = Raw_Data.AGC_c;
daytime = Raw_Data.daytime;
u = Raw_Data.u_;
DateTime_CUP_cell = Raw_Data.DateTime;
DateTime_CUP_raw = datetime(DateTime_CUP_cell,'InputFormat','dd/MM/yyyy HH:mm');
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Figure script
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% NEE vs. PAR, VPD color, summer and winter panels
dark_lowVPD = [0.2314    0.4431    0.3373];
light_lowVPD = [0.8941    0.9412    0.9020];
dark_mediumVPD = [0.1647    0.3843    0.2745];
light_mediumVPD =[0.8392    0.9098    0.8510];
dark_highVPD = [0.0706    0.2118    0.1412];
light_highVPD = [0.7569    0.8667    0.7765];


figure;
subplot(1,2,1); hold on;
these1 = find((month(DateTime_CUP) >= 10 | month(DateTime_CUP) <= 2) & qc == 0 & qc_Sc == 0 & AGC_c == 0 & u > 0.2 & daytime == 1 & VPD < 1.5);
these2 = find((month(DateTime_CUP) >= 10 | month(DateTime_CUP) <= 2) & qc == 0 & qc_Sc == 0 & AGC_c == 0 & u > 0.2 & daytime == 1 & VPD >= 1.5 & VPD < 3);
these3 = find((month(DateTime_CUP) >= 10 | month(DateTime_CUP) <= 2) & qc == 0 & qc_Sc == 0 & AGC_c == 0 & u > 0.2 & daytime == 1 & VPD >= 3);
scatter(PAR(these1),NEE_c(these1),20,light_lowVPD,'filled'); %caxis([nanmin(VPD_30) nanmax(VPD_30)]); colormap(customap);
scatter(PAR(these2),NEE_c(these2),20,light_mediumVPD,'filled');
scatter(PAR(these3),NEE_c(these3),20,light_highVPD,'filled');
plot([0 2500],[0 0],'k');
binplot(PAR(these1),NEE_c(these1),4,dark_lowVPD,10);
binplot(PAR(these2),NEE_c(these2),4,dark_mediumVPD,10);
binplot(PAR(these3),NEE_c(these3),4,dark_highVPD,10);
axis([0 2500 -20 10]); ax = gca; set(ax,'FontSize',12); box off;
ax.XTick = 0:500:2500;
% set(gca,'Position',[0 0 1 1]);
ylabel('NEE (\mumol m^-^2 s^-^1)');
xlabel('PAR (\mumol m^-^2 s^-^1)');
title('Summer (1^s^t Oct - 28^t^h Feb)');
text(0.90,0.98,'(a)','Units', 'Normalized', 'VerticalAlignment', 'Top');
L = legend('VPD < 1.5','VPD [1.5 3]','VPD > 3');
L.Location = 'southwest';

subplot(1,2,2); hold on;
these1 = find(month(DateTime_CUP) >= 5 & month(DateTime_CUP) <= 8 & qc == 0 & qc_Sc == 0 & AGC_c == 0 & u > 0.2 & daytime == 1 & VPD < 1.5);
these2 = find(month(DateTime_CUP) >= 5 & month(DateTime_CUP) <= 8 & qc == 0 & qc_Sc == 0 & AGC_c == 0 & u > 0.2 & daytime == 1 & VPD >= 1.5);
%these3 = find(month(t_s) >= 5 & month(t_s) <= 8 & qc == 0 & qc_Sc == 0 & AGC_c == 0 & u > 0.2 & daytime == 1);
scatter(PAR(these1),NEE_c(these1),20,light_lowVPD,'filled'); %colormap(customap); caxis([nanmin(VPD_30) nanmax(VPD_30)]);
scatter(PAR(these2),NEE_c(these2),20,light_mediumVPD,'filled');
plot([0 2500],[0 0],'k');
binplot(PAR(these1),NEE_c(these1),4,dark_lowVPD,10);
binplot(PAR(these2),NEE_c(these2),4,dark_mediumVPD,10);
axis([0 2500 -20 10]); ax = gca; set(ax,'FontSize',12); box off;
set(gca,'yticklabel',[]); ax.XTick = 0:500:2500;
% set(gca,'Position',[0 0 1 1]);
xlabel('PAR (\mumolm^-^2s^-^1)');
%caxis([nanmin(VPD_30) nanmax(VPD_30)]); %colorbar; h = colorbar; ylabel(h, 'VPD (kPa)');
title('Winter (1^s^t May - 31^t^h Aug)');
text(0.90,0.98,'(b)','Units', 'Normalized', 'VerticalAlignment', 'Top');
