% Clear the workspace
clear;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Load input data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% InpuDateTime_CUP to load
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

% Daily data from Summary .xls
source_summary = 'Input\CumberlandPlain_2014_L6_EP_moderate_Summary.xlsx';
Data_sum = readtable(source_summary,'Sheet','Daily (all)','Range','A3:AI1097','ReadVariableNames',false);
daily_date_cell = Data_sum.Var1;
daily_date = datetime(daily_date_cell,'InputFormat','dd/MM/yyyy hh:mm:ss a');
ER_solo_daily = Data_sum.Var7;
ET_daily = Data_sum.Var8;
GPP_solo_daily = -Data_sum.Var20;
NEE_solo_daily = Data_sum.Var23;
Precip_daily = Data_sum.Var27;
clearvars daily_date_cell Data_sum source_summary;

% LAI data
% Load raw data (.csv file)
% Location of raw input file
source_LAI = 'Input\CUP_LAI_20161021.csv';
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

% Time series plot, C flux, Drivers and canopy dynamic

% Alexis, 25/11/2016

% !! For some weird reason, need to run the script in 2 part

% 1. fluxes bigger subplot, try different symbol instead of color (triangle,
% circle, square for example), smaller marker size
% 2. PAR in a smaller subplot?
LF_avg = [0.22 0.26 0.57 0.26 0.38 0.38 0.55 0.03  0.03 0.15 0.07 0.07];
LF_std = [0.04 0.06 0.11 0.04 0.08 0.08 0.13 0.01 0.01 0.05 0.02 0.02];
for i = 1:5
    LF_date(i) = datetime(2015,i+7,15);
end
for i = 1:7
    LF_date(i+5) = datetime(2016,i,15);
end
Panel_num = {'a','b','c','d','e'};
yaxismin = [-8 0 0 0 0.6];
yaxismax = [8 120 40 3000 1.1];

figure;
for i = 1:5
    subplot(5,1,i); hold on;
    for yeari = 2013:2017
        x1 = datenum(datetime(yeari,05,01));
        x2 = datenum(datetime(yeari,08,31));
        x = [x1 x2 x2 x1];
        y = [yaxismin(i) yaxismin(i) yaxismax(i) yaxismax(i)];
        patch(x,y,'k','FaceAlpha',.2,'EdgeColor','none');
        x1 = datenum(datetime(yeari,10,01));
        x2 = datenum(datetime(yeari+1,02,30));
        x = [x1 x2 x2 x1];
        y = [yaxismin(i) yaxismin(i) yaxismax(i) yaxismax(i)];
        patch(x,y,'k','FaceAlpha',.5,'EdgeColor','none');
    end
    text(0.90,0.98,Panel_num(i),'Units', 'Normalized', 'VerticalAlignment', 'Top');
end

subplot(5,1,1); hold on; 
ERplot = plot(daily_date,ER_solo_daily,'LineStyle','none','Marker','o', ...
    'MarkerFaceColor','r','MarkerEdgeColor','none','MarkerSize',4);
GPPplot = plot(daily_date,GPP_solo_daily,'LineStyle','none','Marker','o', ...
    'MarkerFaceColor','g','MarkerEdgeColor','none','MarkerSize',4);
NEEplot = plot(daily_date,NEE_solo_daily,'LineStyle','none','Marker','o', ...
    'MarkerFaceColor','k','MarkerEdgeColor','none','MarkerSize',4);
ylabel('Flux (gCm^-^2d^-^1)');
%legend([NEEplot ERplot GPPplot],'NEE','ER','GPP','Location','northouDateTime_CUPide','Orientation','horizontal');
title('Daily C flux');

subplot(5,1,2); 
hold on;
plot(daily_date,Precip_daily,'k'); ylabel('Rain (mm)');
title('Drivers'); 

subplot(5,1,3);
hold on;
plot(DateTime_CUP,Ta,'k'); ylabel('T_a_i_r (°C)');

subplot(5,1,4);
hold on;
plot(DateTime_CUP,PAR,'k'); ylabel('PAR (\mumolm^-^2s^-^1)');

subplot(5,1,5);
hold on;
binplot(datenum(DateTime_LAI),LAI,6,'k');
title('Canopy dynamic'); ylabel('LAI (m^-^2m^-^2)');

i = 1;
for y = 1:3
    for m = 1:12
        month_dt(i) = datetime(2013+y,0+m,1);
        i = i+1;
    end
end
month_dt(37) = datetime(2017,1,1);

for i = 1:5
    subplot(5,1,i); hold on;
    plot([min(datenum(DateTime_CUP)) max(datenum(DateTime_CUP))], [0 0],'k');
    plot([min(datenum(DateTime_CUP)) max(datenum(DateTime_CUP))], [0 0],'k');
    axis([datenum(datetime(2014,1,1)) datenum(datetime(2017,1,1)) yaxismin(i) yaxismax(i)]);
    ax = gca; set(ax,'FontSize',14); box on;
    ax.XTick = datenum(month_dt);
    datetick('x','myy','keeplimits','keepticks');
    if i ~= 5
    	set(gca,'xticklabel',[]);
    end
end

% For some weird reason, I have to run the following part AFTER the first
% part

subplot(5,1,2); hold on;
ax = gca;
ax1_pos = ax.Position; % position of first axes
ax2 = axes('Position',ax1_pos,...
'XAxisLocation','top',...
'YAxisLocation','right',...
'Color','none');
ax2.YColor = 'r';
line(datenum(DateTime_CUP),SWC_s,'Parent',ax2,'Color','r','LineStyle','-'); 
axis([min(datenum(DateTime_CUP)) max(datenum(DateTime_CUP)) 0 40]); set(ax2,'FontSize',16); % ax2.YTick = 0:0.5:3;
set(gca,'xticklabel',[]); %axes(ax); % set(gca,'yticklabel',[]); % ylabel('VPD (kPa)');
ylabel('SWC _0_-_8_c_m (%)');

subplot(5,1,3); hold on;
ax = gca;
ax1_pos = ax.Position; % position of first axes
ax2 = axes('Position',ax1_pos,...
'XAxisLocation','top',...
'YAxisLocation','right',...
'Color','none');
ax2.YColor = 'r';
line(datenum(DateTime_CUP),VPD,'Parent',ax2,'Color','r','LineStyle','-');
axis([min(datenum(DateTime_CUP)) max(datenum(DateTime_CUP)) 0 10]); set(ax2,'FontSize',16); % ax2.YTick = 0:0.5:3;
set(gca,'xticklabel',[]); %axes(ax); % set(gca,'yticklabel',[]); % ylabel('VPD (kPa)');
ylabel('VPD (kPa)');

subplot(5,1,5); hold on;
ax = gca;
ax1_pos = ax.Position; % position of first axes
ax2 = axes('Position',ax1_pos,...
'XAxisLocation','top',...
'YAxisLocation','right',...
'Color','none');
ax2.YColor = 'r';
line(datenum(LF_date),LF_avg,'Parent',ax2,'Color','r','LineStyle','-');
line(datenum(LF_date),LF_avg + LF_std,'Parent',ax2,'Color','r','LineStyle','--');
line(datenum(LF_date),LF_avg - LF_std,'Parent',ax2,'Color','r','LineStyle','--');
axis([min(datenum(DateTime_CUP)) max(datenum(DateTime_CUP)) 0 0.7]); set(ax2,'FontSize',16); % ax2.YTick = 0:0.5:3;
set(gca,'xticklabel',[]); %axes(ax); % set(gca,'yticklabel',[]); % ylabel('VPD (kPa)');
ylabel('Litter fall (gb^-^1d^-^1)');



