% Clear the workspace
clear;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Load input data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Inputs:
% Month_d 
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
Mean_PAR = nan(n,1);
Mean_SWC = nan(n,1);
d=1;
for yeari = 2013:2016
    for monthi = 1:12
        for dayi = 1:31
            if isempty(find(Year_t == yeari & Month_t == monthi & Day_t == dayi,1)) == 0 && d <= n
                Max_VPD(d,1) = max(VPD_30(Year_t == yeari & Month_t == monthi & Day_t == dayi));
                Mean_PAR(d,1) = mean(PAR(Year_t == yeari & Month_t == monthi & Day_t == dayi));
                Mean_SWC(d,1) = mean(SWC_s(Year_t == yeari & Month_t == monthi & Day_t == dayi));
                d = d+1;
            end
        end
    end
end

% Define condition 
these_summer_dry_L = find((Month_d >= 10 | Month_d <= 2) & Mean_SWC < nanmedian(Mean_SWC) & Max_VPD < 1.5);
these_summer_dry_M = find((Month_d >= 10 | Month_d <= 2) & Mean_SWC < nanmedian(Mean_SWC) & Max_VPD >= 1.5 & Max_VPD < 3);
these_summer_dry_H = find((Month_d >= 10 | Month_d <= 2) & Mean_SWC < nanmedian(Mean_SWC) & Max_VPD >= 3);
these_summer_wet_L = find((Month_d >= 10 | Month_d <= 2) & Mean_SWC > nanmedian(Mean_SWC) & Max_VPD < 1.5);
these_summer_wet_M = find((Month_d >= 10 | Month_d <= 2) & Mean_SWC > nanmedian(Mean_SWC) & Max_VPD >= 1.5 & Max_VPD < 3);
these_summer_wet_H = find((Month_d >= 10 | Month_d <= 2) & Mean_SWC > nanmedian(Mean_SWC) & Max_VPD >= 3);
these_winter_dry_L = find(Month_d >= 5 & Month_d <= 8 & Mean_SWC < nanmedian(Mean_SWC) & Max_VPD < 1.5);
these_winter_dry_M = find(Month_d >= 5 & Month_d <= 8 & Mean_SWC < nanmedian(Mean_SWC) & Max_VPD >= 1.5 & Max_VPD < 3);
these_winter_dry_H = find(Month_d >= 5 & Month_d <= 8 & Mean_SWC < nanmedian(Mean_SWC) & Max_VPD >= 3);
these_winter_wet_L = find(Month_d >= 5 & Month_d <= 8 & Mean_SWC > nanmedian(Mean_SWC) & Max_VPD < 1.5);
these_winter_wet_M = find(Month_d >= 5 & Month_d <= 8 & Mean_SWC > nanmedian(Mean_SWC) & Max_VPD >= 1.5 & Max_VPD < 3);
these_winter_wet_H = find(Month_d >= 5 & Month_d <= 8 & Mean_SWC > nanmedian(Mean_SWC) & Max_VPD >= 3);

% Figure
light_blue = [0.8 0.8 1]; light_red = [1 0.8 0.8]; light_green = [0.8 1 0.8];
figure; subplot(2,2,1); hold on; % Summer, NEE vs. PAR, daily, wet
scatter(Mean_PAR(these_summer_wet_L),NEE_d(these_summer_wet_L),20,'MarkerFaceColor',light_blue,'MarkerEdgeColor','none'); % colormap(maptest); axis([0 9 -2 4]);
scatter(Mean_PAR(these_summer_wet_M),NEE_d(these_summer_wet_M),20,'MarkerFaceColor',light_green,'MarkerEdgeColor','none');
scatter(Mean_PAR(these_summer_wet_H),NEE_d(these_summer_wet_H),20,'MarkerFaceColor',light_red,'MarkerEdgeColor','none');
plot([0 1000],[0 0],'k','LineWidth',1);
binplot(Mean_PAR(these_summer_wet_L),NEE_d(these_summer_wet_L),4,'b');
binplot(Mean_PAR(these_summer_wet_M),NEE_d(these_summer_wet_M),4,'g');
binplot(Mean_PAR(these_summer_wet_H),NEE_d(these_summer_wet_H),4,'r');
ylabel ('Daily NEE (gC m^-^2 d^-^1)'); set(gca,'FontSize',16);
set(gca,'xticklabel',[]); title('Summer (1^s^t Oct - 28^t^h Feb)');
axis([0 1000 -5 5]);
subplot(2,2,3); hold on; % Summer, NEE vs. PAR, daily, dry
scatter(Mean_PAR(these_summer_dry_L),NEE_d(these_summer_dry_L),20,'MarkerFaceColor',light_blue,'MarkerEdgeColor','none'); % colormap(maptest); axis([0 9 -2 4]);
scatter(Mean_PAR(these_summer_dry_M),NEE_d(these_summer_dry_M),20,'MarkerFaceColor',light_green,'MarkerEdgeColor','none');
scatter(Mean_PAR(these_summer_dry_H),NEE_d(these_summer_dry_H),20,'MarkerFaceColor',light_red,'MarkerEdgeColor','none');
plot([0 1000],[0 0],'k','LineWidth',1);
binplot(Mean_PAR(these_summer_dry_L),NEE_d(these_summer_dry_L),4,'b');
binplot(Mean_PAR(these_summer_dry_M),NEE_d(these_summer_dry_M),4,'g');
binplot(Mean_PAR(these_summer_dry_H),NEE_d(these_summer_dry_H),4,'r');
ylabel ('Daily NEE (gC m^-^2 d^-^1)');
xlabel ('Mean PAR (\mumol m^-^2 s^-^1)');
set(gca,'FontSize',16);
axis([0 1000 -5 5]);
subplot(2,2,2); hold on; % winter, NEE vs. PAR, daily, wet
scatter(Mean_PAR(these_winter_wet_L),NEE_d(these_winter_wet_L),20,'MarkerFaceColor',light_blue,'MarkerEdgeColor','none'); % colormap(maptest); axis([0 9 -2 4]);
scatter(Mean_PAR(these_winter_wet_M),NEE_d(these_winter_wet_M),20,'MarkerFaceColor',light_green,'MarkerEdgeColor','none');
scatter(Mean_PAR(these_winter_wet_H),NEE_d(these_winter_wet_H),20,'MarkerFaceColor',light_red,'MarkerEdgeColor','none');
plot([0 1000],[0 0],'k','LineWidth',1);
binplot(Mean_PAR(these_winter_wet_L),NEE_d(these_winter_wet_L),4,'b');
binplot(Mean_PAR(these_winter_wet_M),NEE_d(these_winter_wet_M),4,'g');
binplot(Mean_PAR(these_winter_wet_H),NEE_d(these_winter_wet_H),4,'r');
set(gca,'yticklabel',[]); set(gca,'FontSize',16);
set(gca,'xticklabel',[]); title('Winter (1^s^t May - 31^t^h Aug)');
axis([0 1000 -5 5]);
subplot(2,2,4); hold on; % winter, NEE vs. PAR, daily, dry
scatter(Mean_PAR(these_winter_dry_L),NEE_d(these_winter_dry_L),20,'MarkerFaceColor',light_blue,'MarkerEdgeColor','none'); % colormap(maptest); axis([0 9 -2 4]);
scatter(Mean_PAR(these_winter_dry_M),NEE_d(these_winter_dry_M),20,'MarkerFaceColor',light_green,'MarkerEdgeColor','none');
scatter(Mean_PAR(these_winter_dry_H),NEE_d(these_winter_dry_H),20,'MarkerFaceColor',light_red,'MarkerEdgeColor','none');
plot([0 1000],[0 0],'k','LineWidth',1);
binplot(Mean_PAR(these_winter_dry_L),NEE_d(these_winter_dry_L),4,'b');
binplot(Mean_PAR(these_winter_dry_M),NEE_d(these_winter_dry_M),4,'g');
binplot(Mean_PAR(these_winter_dry_H),NEE_d(these_winter_dry_H),4,'r');
xlabel ('Mean PAR (\mumol m^-^2 s^-^1)');
set(gca,'FontSize',16);
set(gca,'yticklabel',[]); 
axis([0 1000 -5 5]);

for i = 1:4
    if i == 1 || i == 2 
        subplot(2,2,i); hold on;
        text(0.70,0.98,'Wet soil','Units', 'Normalized', 'VerticalAlignment', 'Top');
    end
    if i == 3 || i == 4 
        subplot(2,2,i); hold on;
        text(0.70,0.98,'Dry soil','Units', 'Normalized', 'VerticalAlignment', 'Top');
    end
end

labelp = {'(a)','(b)','(c)','(d)'};
for i = 1:4
    subplot(2,2,i); hold on;
    text(0.90,0.98,labelp(i),'Units', 'Normalized', 'VerticalAlignment', 'Top');
    sub_pos = get(gca,'position'); % get subplot axis position
    set(gca,'position',sub_pos.*[1 1 1.1 1.1]) % stretch its width and height
end

