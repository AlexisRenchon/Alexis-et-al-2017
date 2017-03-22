% Clear the workspace
clear;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Load input data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Inputs:
% Month_d 
% ta_30 (half-hourly air temperature at 30m)
% VPD_30 (half-hourly VPD at 30m)
% PAR (half-hourly PAR at 30m)
% SWC_s (half-hourly soil water content 0-8cm)
% Date_daily 
% day_t (half-hourly day vector)
% month_t (half-hourly month vector)
% year_t (half-hourly year vector)
% NEE_d (daily NEE)

% Daily data from Summary .xls
source_summary = 'Input\CumberlandPlain_2014_L6_EP_moderate_Summary.xlsx';
Data_sum = readtable(source_summary,'Sheet','Daily (all)','Range','A3:AL1097','ReadVariableNames',false);
Date_daily = Data_sum.Var1; Date_daily = datetime(Date_daily,'InputFormat','dd/MM/yyyy hh:mm:ss a');
ET_daily = Data_sum.Var8;
NEE_d = Data_sum.Var23;
Month_d = Data_sum.Var36;
clearvars Data_sum source_summary;

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
ta_30 = ncread(source_processed,'Ta'); ta_30 = reshape(ta_30,[],1); 
SWC_s = ncread(source_processed,'Sws'); SWC_s = reshape(SWC_s,[],1); SWC_s = SWC_s*100;
VPD_30 = ncread(source_processed,'VPD'); VPD_30 = reshape(VPD_30,[],1); 
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

% Daily min max and mean of variables
n = length(Date_daily);
Max_VPD = nan(n,1);
Min_Ta = nan(n,1);
Mean_PAR = nan(n,1);
Mean_SWC = nan(n,1);
d=1;
for yeari = 2013:2016
    for monthi = 1:12
        for dayi = 1:31
            if isempty(find(Year_t == yeari & Month_t == monthi & Day_t == dayi,1)) == 0 && d <= n
                Max_VPD(d,1) = max(VPD_30(Year_t == yeari & Month_t == monthi & Day_t == dayi));
                Min_Ta(d,1) = min(ta_30(Year_t == yeari & Month_t == monthi & Day_t == dayi));
                Mean_PAR(d,1) = mean(PAR(Year_t == yeari & Month_t == monthi & Day_t == dayi));
                Mean_SWC(d,1) = mean(SWC_s(Year_t == yeari & Month_t == monthi & Day_t == dayi));
                d = d+1;
            end
        end
    end
end

wintercolor = [0.7 0.7 1];
summercolor = [1 0.7 0.7];

%plots
these_winter = find(Month_d >= 5 & Month_d <= 8);
these_summer = find(Month_d >= 10 | Month_d <= 2);

figure; subplot(1,4,1); hold on; 
plot([0 25],[0 0],'k','LineWidth',1);
scatter(Min_Ta(these_winter),NEE_d(these_winter),'Marker','o','MarkerEdgeColor',wintercolor,'LineWidth',0.5,'SizeData',10);
scatter(Min_Ta(these_summer),NEE_d(these_summer),'Marker','o','MarkerEdgeColor',summercolor,'LineWidth',0.5,'SizeData',10);
xlabel('Min Ta (°C)'); ylabel ('NEE (gC m^-^2 d^-^1)'); ax = gca; set(ax,'FontSize',12);
binplot(Min_Ta(these_winter),NEE_d(these_winter),12,'b');
binplot(Min_Ta(these_summer),NEE_d(these_summer),12,'r');
axis([0 25 -6 6]); ax.XTick = 0:5:25; ax.YTick = -6:2:6;

subplot(1,4,2); hold on; 
plot([0 10],[0 0],'k','LineWidth',1); 
scatter(Max_VPD(these_winter),NEE_d(these_winter),'Marker','o','MarkerEdgeColor',wintercolor,'LineWidth',0.5,'SizeData',10);
scatter(Max_VPD(these_summer),NEE_d(these_summer),'Marker','o','MarkerEdgeColor',summercolor,'LineWidth',0.5,'SizeData',10);
xlabel('Max VPD (kPa)'); ylabel ('NEE (gC m^-^2 d^-^1)'); ax = gca; set(ax,'FontSize',12);
binplot(Max_VPD(these_winter),NEE_d(these_winter),12,'b');
binplot(Max_VPD(these_summer),NEE_d(these_summer),12,'r');
axis([0 10 -6 6]); ax.XTick = 0:2:10; ax.YTick = -6:2:6; ax.YTickLabel = []; ylabel([]);

subplot(1,4,3); hold on;
plot([0 40000],[0 0],'k','LineWidth',1); 
scatter(Mean_PAR(these_winter),NEE_d(these_winter),'Marker','o','MarkerEdgeColor',wintercolor,'LineWidth',0.5,'SizeData',10);
scatter(Mean_PAR(these_summer),NEE_d(these_summer),'Marker','o','MarkerEdgeColor',summercolor,'LineWidth',0.5,'SizeData',10);
binplot(Mean_PAR(these_winter),NEE_d(these_winter),12,'b');
binplot(Mean_PAR(these_summer),NEE_d(these_summer),12,'r');
xlabel('Mean PAR (\mumol m^-^2 s^-^1)'); ylabel ('NEE (gC m^-^2 d^-^1)'); ax = gca; set(ax,'FontSize',12);
axis([0 1000 -6 6]); ax.XTick = 0:200:1000; ax.YTick = -6:2:6; ax.YTickLabel = []; ylabel([]);

subplot(1,4,4); hold on; 
plot([0 40],[0 0],'k','LineWidth',1);
scatter(Mean_SWC(these_winter),NEE_d(these_winter),'Marker','o','MarkerEdgeColor',wintercolor,'LineWidth',0.5,'SizeData',10);
scatter(Mean_SWC(these_summer),NEE_d(these_summer),'Marker','o','MarkerEdgeColor',summercolor,'LineWidth',0.5,'SizeData',10);
binplot(Mean_SWC(these_winter),NEE_d(these_winter),12,'b');
binplot(Mean_SWC(these_summer),NEE_d(these_summer),12,'r');
xlabel('Mean SWC (%)'); ylabel ('NEE (gC m^-^2 d^-^1)'); ax = gca; set(ax,'FontSize',12);
axis([0 40 -6 6]); ax.XTick = 0:10:40; ax.YTick = -6:2:6; ax.YTickLabel = []; ylabel([]);

labelp = {'(a)','(b)','(c)','(d)'};
for i = 1:4
    subplot(1,4,i); hold on;
    text(0.90,0.98,labelp(i),'Units', 'Normalized', 'VerticalAlignment', 'Top');
    sub_pos = get(gca,'position'); % get subplot axis position
    set(gca,'position',sub_pos.*[1 1 1.2 1]) % stretch its width and height
end
