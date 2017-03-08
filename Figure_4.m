% Clear the workspace
clear;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Load input data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Inputs:
% daytime (half-hourly 0 for night, 1 for day)
% qc (half-hourly qc Mauder and Foken 2004)
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


% Figure 4, v1
% Alexis, 29/11/2016

% 2 seasons, winter and summer
% 10 quantiles/bins of drivers; VPD, SWC and tsoil 
% parameters of Mitscherlich for each quantile/bin, along with CI
% plot

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
            driver = SWC_s(this_season);
        elseif d == 3
            driver = VPD_30(this_season);
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
subplot(3,3,4); hold on; ylabel('N_e_s (\mumolm^-^2s^-^1)');
ax = gca;  ax.YTick = 0:5:20;
subplot(3,3,5); hold on; set(gca,'yticklabel',[]);
subplot(3,3,6); hold on; set(gca,'yticklabel',[]);
subplot(3,3,7); hold on; ylabel('R_d (\mumolm^-^2s^-^1)'); xlabel('Tsoil_5_c_m (°C)');
subplot(3,3,8); hold on; xlabel('SWC_0_-_8_c_m (%)'); set(gca,'yticklabel',[]);
subplot(3,3,9); hold on; xlabel('VPD_3_0_m (kPa)');set(gca,'yticklabel',[]);
for i = 1:3
    subplot(3,3,i); hold on;
    ylim([0 0.04]);
    set(gca,'xticklabel',[]);
    set(gca,'FontSize',16);
end
for i = 4:6
    subplot(3,3,i); hold on;
    ylim([0 20]);
    set(gca,'xticklabel',[]);
    set(gca,'FontSize',16);
end
for i = 7:9
    subplot(3,3,i); hold on;
    ylim([0 8]);
    set(gca,'FontSize',16);
end
labelp = {'(a)','(b)','(c)','(d)','(e)','(f)','(g)','(h)','(i)'};
for i = 1:9
    subplot(3,3,i); hold on;
    text(0.90,0.98,labelp(i),'Units', 'Normalized', 'VerticalAlignment', 'Top');
    sub_pos = get(gca,'position'); % get subplot axis position
    set(gca,'position',sub_pos.*[1 1 1.2 1.2]) % stretch its width and height
end

