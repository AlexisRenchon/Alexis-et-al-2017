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
% ta_30 (half-hourly air temperature)

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
ta_30 = ncread(source_processed,'Ta'); ta_30 = reshape(ta_30,[],1); 
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

% NEE, ET and VPD diurnal trend, + PAR and tsoil, winter and summer
% Alexis, 02/12/2016

% Create half-hour vector (0-47), starting at half-hour = 1
i = 1;
j = 1;
n = length(NEE_c);
hhour_t = NaN(n,1);
while i < n + 1
    if j < 47
        hhour_t(i) = j;
        j = j + 1;
    else
        hhour_t(i) = 47;
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
for i = 2:48 % start 1st half-hour of the day, hhour_t = 1 corresponding to flux 00:00:00 to 00:30:00, end at hhour_t = 47, corresponding to 23:00:00 to 23:30:00, i.e. we don't look 23:30:30 to 00:00:00 
    these_L = find(hhour_t == i-1 & (month(t_s) >= 10 | month(t_s) <= 2)); % 1 Nov through 30 Apr
    these_L_qc = find(hhour_t==i-1 & (month(t_s) >= 10 | month(t_s) <= 2) & qc ~= 2 & qc_Sc ~= 2 & AGC_c == 0 & u > 0.2);
    these_NL = find(hhour_t==i-1 & month(t_s) >= 5 & month(t_s) <= 8); % 1 May through 31 Oct
    these_NL_qc = find(hhour_t==i-1 & month(t_s) >= 5 & month(t_s) <= 8 & qc ~= 2 & qc_Sc ~= 2 & AGC_c == 0 & u > 0.2);
    these_h2o_L_qc = find(hhour_t==i-1 & (month(t_s) >= 10 | month(t_s) <= 2) & qc_h2o_flux == 0 & AGC_c == 0);
    these_h2o_NL_qc = find(hhour_t==i-1 & month(t_s) >= 5 & month(t_s) <= 8 & qc_h2o_flux == 0 & AGC_c == 0);
    Diurn_quantiles_VPD_L(i-1,[1 2 3]) = quantile(VPD_30(these_L),[0.25, 0.5, 0.75]);
    Diurn_quantiles_VPD_NL(i-1,[1 2 3]) = quantile(VPD_30(these_NL),[0.25, 0.5, 0.75]);
    Diurn_quantiles_PAR_L(i-1,[1 2 3]) = quantile(PAR(these_L),[0.25, 0.5, 0.75]);
    Diurn_quantiles_PAR_NL(i-1,[1 2 3]) = quantile(PAR(these_NL),[0.25, 0.5, 0.75]);
    Diurn_quantiles_Ta_L(i-1,[1 2 3]) = quantile(ta_30(these_L),[0.25, 0.5, 0.75]);
    Diurn_quantiles_Ta_NL(i-1,[1 2 3]) = quantile(ta_30(these_NL),[0.25, 0.5, 0.75]);
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
subplot(3,2,[1 3]); hold on;

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
line(x,Diurn_quantiles_VPD_L(:,2),'Parent',ax2,'Color','r','LineStyle','-');
line(x,Diurn_quantiles_VPD_L(:,1),'Parent',ax2,'Color','r','LineStyle',':');
line(x,Diurn_quantiles_VPD_L(:,3),'Parent',ax2,'Color','r','LineStyle',':');
axis([0 24 0 3]); set(ax2,'FontSize',16); ax2.YTick = 0:0.5:3;
set(gca,'xticklabel',[]); set(gca,'yticklabel',[]); % ylabel('VPD (kPa)');

subplot(3,2,[2 4]); hold on;
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
line(x,Diurn_quantiles_VPD_NL(:,2),'Parent',ax2,'Color','r','LineStyle','-');
line(x,Diurn_quantiles_VPD_NL(:,1),'Parent',ax2,'Color','r','LineStyle',':');
line(x,Diurn_quantiles_VPD_NL(:,3),'Parent',ax2,'Color','r','LineStyle',':');
axis([0 24 0 3]); set(ax2,'FontSize',16); ax2.YTick = 0:0.5:3;
set(gca,'xticklabel',[]); ylabel('VPD (kPa)');

% subplot(3,2,5); hold on;
% i = find((month_t >= 10 | month_t <= 2) & hour(t_s) >= 12 & daytime == 1);
% scatter(PAR(i),NEE_ff_day(i),15,'MarkerFaceColor',light_red,'MarkerEdgeColor','none'); 
% i = find((month_t >= 10 | month_t <= 2) & hour(t_s) < 12 & daytime == 1); 
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
% i = find(month_t >= 5 & month_t <= 8 & hour(t_s) >= 12 & daytime == 1);
% scatter(PAR(i),NEE_ff_day(i),15,'MarkerFaceColor',light_red,'MarkerEdgeColor','none'); 
% i = find(month_t >= 5 & month_t <= 8 & hour(t_s) < 12 & daytime == 1);
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
% these1 = find((month(t_s) >= 10 | month(t_s) <= 2) & qc == 0 & qc_Sc == 0 & AGC_c == 0 & u > 0.2 & daytime == 1 & VPD_30 < 1.5);
% these2 = find((month(t_s) >= 10 | month(t_s) <= 2) & qc == 0 & qc_Sc == 0 & AGC_c == 0 & u > 0.2 & daytime == 1 & VPD_30 >= 1.5 & VPD_30 < 3);
% these3 = find((month(t_s) >= 10 | month(t_s) <= 2) & qc == 0 & qc_Sc == 0 & AGC_c == 0 & u > 0.2 & daytime == 1 & VPD_30 >= 3);
% scatter(PAR(these1),NEE_c(these1),10,'b','filled'); %caxis([nanmin(VPD_30) nanmax(VPD_30)]); colormap(customap);
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
% these1 = find(month(t_s) >= 5 & month(t_s) <= 8 & qc == 0 & qc_Sc == 0 & AGC_c == 0 & u > 0.2 & daytime == 1 & VPD_30 < 1.5);
% these2 = find(month(t_s) >= 5 & month(t_s) <= 8 & qc == 0 & qc_Sc == 0 & AGC_c == 0 & u > 0.2 & daytime == 1 & VPD_30 >= 1.5);
% %these3 = find(month(t_s) >= 5 & month(t_s) <= 8 & qc == 0 & qc_Sc == 0 & AGC_c == 0 & u > 0.2 & daytime == 1);
% scatter(PAR(these1),NEE_c(these1),10,'b','filled'); %colormap(customap); caxis([nanmin(VPD_30) nanmax(VPD_30)]);
% scatter(PAR(these2),NEE_c(these2),10,'g','filled');
% axis([0 2500 -20 10]); ax = gca; set(ax,'FontSize',12); box off; set(gca,'yticklabel',[]);
% % set(gca,'Position',[0 0 1 1]);
% xlabel('PAR (\mumolm^-^2s^-^1)');
% hold on; plot([0 2500],[0 0],'k');
% %caxis([nanmin(VPD_30) nanmax(VPD_30)]); %colorbar; h = colorbar; ylabel(h, 'VPD (kPa)');
% text(0.90,0.98,'(f)','Units', 'Normalized', 'VerticalAlignment', 'Top');

% h = colorbar;
% ylabel(h, 'VPD (kPa)');





% 



% NEE vs. PAR, VPD color, summer and winter panels

lowVPD = [0.7569    0.8667    0.7765];
mediumVPD = [0.2314    0.4431    0.3373];
highVPD = [0.0706    0.2118    0.1412];
morningc = [0.6 0.6 1];
afternoonc = [1 0.8 0.8];



subplot(3,2,5); hold on;
these1 = find((month(DateTime_CUP) >= 10 | month(DateTime_CUP) <= 2) & qc == 0 & qc_Sc == 0 & AGC_c == 0 & u > 0.2 & daytime == 1 & VPD_30 < 1.5);
these2 = find((month(DateTime_CUP) >= 10 | month(DateTime_CUP) <= 2) & qc == 0 & qc_Sc == 0 & AGC_c == 0 & u > 0.2 & daytime == 1 & VPD_30 >= 1.5 & VPD_30 < 3);
these3 = find((month(DateTime_CUP) >= 10 | month(DateTime_CUP) <= 2) & qc == 0 & qc_Sc == 0 & AGC_c == 0 & u > 0.2 & daytime == 1 & VPD_30 >= 3);
scatter(PAR(these1),NEE_c(these1),20,lowVPD,'filled'); %caxis([nanmin(VPD_30) nanmax(VPD_30)]); colormap(customap);
scatter(PAR(these2),NEE_c(these2),20,mediumVPD,'filled');
scatter(PAR(these3),NEE_c(these3),20,highVPD,'filled');
i = find((month_t >= 10 | month_t <= 2) & hour(t_s) < 12 & daytime == 1);
binplot(PAR(i),NEE_ff_day(i),15,morningc,10);
i = find((month_t >= 10 | month_t <= 2) & hour(t_s) >= 12 & daytime == 1);
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


subplot(3,2,6); hold on;
these1 = find(month(DateTime_CUP) >= 5 & month(DateTime_CUP) <= 8 & qc == 0 & qc_Sc == 0 & AGC_c == 0 & u > 0.2 & daytime == 1 & VPD_30 < 1.5);
these2 = find(month(DateTime_CUP) >= 5 & month(DateTime_CUP) <= 8 & qc == 0 & qc_Sc == 0 & AGC_c == 0 & u > 0.2 & daytime == 1 & VPD_30 >= 1.5);
%these3 = find(month(t_s) >= 5 & month(t_s) <= 8 & qc == 0 & qc_Sc == 0 & AGC_c == 0 & u > 0.2 & daytime == 1);
dotsa = scatter(PAR(these1),NEE_c(these1),20,lowVPD,'filled'); %colormap(customap); caxis([nanmin(VPD_30) nanmax(VPD_30)]);
dotsb = scatter(PAR(these2),NEE_c(these2),20,mediumVPD,'filled');
dotsc = scatter(3000,2,20,highVPD,'filled'); % for the legend
dotsd = plot(3000,2,'LineWidth',3,'Marker','o','MarkerEdgeColor', ...
    'none','MarkerFaceColor',morningc,'Color',morningc,'MarkerSize',10);
dotse = plot(3000,2,'LineWidth',3,'Marker','o','MarkerEdgeColor', ...
    'none','MarkerFaceColor',afternoonc,'Color',afternoonc,'MarkerSize',10);
i = find(month_t >= 5 & month_t <= 8 & hour(t_s) < 12 & daytime == 1);
binplot(PAR(i),NEE_ff_day(i),15,morningc,10);
i = find(month_t >= 5 & month_t <= 8 & hour(t_s) >= 12 & daytime == 1);
binplot(PAR(i),NEE_ff_day(i),15,afternoonc,10);
plot([0 2500],[0 0],'k');
% binplot(PAR(these1),NEE_c(these1),4,dark_lowVPD,10);
% binplot(PAR(these2),NEE_c(these2),4,dark_mediumVPD,10);
axis([0 2500 -20 10]); ax = gca; set(ax,'FontSize',16); box off;
set(gca,'yticklabel',[]); ax.XTick = 0:500:2500;
% set(gca,'Position',[0 0 1 1]);
xlabel('PAR (\mumolm^-^2s^-^1)');
%caxis([nanmin(VPD_30) nanmax(VPD_30)]); %colorbar; h = colorbar; ylabel(h, 'VPD (kPa)');
%title('Winter (1^s^t May - 31^t^h Aug)');
text(0.90,0.98,'(b)','Units', 'Normalized', 'VerticalAlignment', 'Top');
L = legend([dotsa dotsb dotsc dotsd dotse],'VPD < 1.5','VPD [1.5 3]','VPD > 3','Morning','Afternoon');
L.Location = 'southeast';
L.FontSize = 8;

