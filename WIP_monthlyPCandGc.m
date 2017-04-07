% Monthly "PC" as median GEP when PAR ~ 1000 and VPD ~ 1
% !! contains all gap-filled GEP. 
% Monthly surface conductance, "Gs_m", as median Gs when PAR ~ 1000 and VPD ~ 1
% !! contains all gap-filled gs, no rain filter. 
PC = nan(12*3,3);
date_m = datetime;
gs_m = nan(12*3,3);
LAI_m = nan(36,1);
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
start_d = datenum(DateTime_LAI(1));
end_d = datenum(DateTime_LAI(297));
hh_num = datenum(DateTime_CUP(2)) - datenum(DateTime_CUP(1));
xq = start_d:hh_num:end_d;
s_LAI = feval(fitresult,xq);

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





use_gc = find(Precip_2daysago < 1 & Precip_1daysago < 0.5 & Precip_halfdayago < 0.2 & qc_h2o_flux == 0 & AGC_c == 0 ...
     & vpd > 1 & vpd < 1.5 & SWC > 0.11 & PAR > 800);
use_pc = find(PAR > 800 & PAR < 1200 & vpd > 1 & vpd < 1.5 & SWC > 0.11);

figure; plot(DateTime_LAI,LAI,'LineStyle','none','Marker','o','MarkerEdgeColor','none','MarkerFaceColor','k','MarkerSize',4); hold on;
plot(DateTime_CUP(use_pc),(GEP(use_pc)/1000)*100,'LineStyle','none','Marker','o','MarkerEdgeColor','none','MarkerFaceColor','r','MarkerSize',4);
plot(DateTime_CUP(use_gc),gc(use_gc)*250,'LineStyle','none','Marker','o','MarkerEdgeColor','none','MarkerFaceColor','b','MarkerSize',4);

use = find(Precip_2daysago < 1 & Precip_1daysago < 0.5 & Precip_halfdayago < 0.2 & qc_h2o_flux == 0 & AGC_c == 0 ...
     & vpd > 1 & vpd < 1.5 & SWC > 0.11 & PAR > 800 & PAR < 1200 & qc == 0 & qc_Sc == 0);

figure; scatter(gc(use),GEP(use));
% TO DO 
% GET DAILY (or weekly?) PC and DAILY GC instead of half hourly, as median over a day
% (48 half-hour), get the value if at least 5 or 6 (?) data points




st_LAI = find(datenum(DateTime_CUP) == start_d);
en_LAI = find(datenum(DateTime_CUP) == end_d);

GEP_sLAI = GEP(st_LAI:en_LAI);
gc_sLAI = gc(st_LAI:en_LAI);
vpd_sLAI = vpd(st_LAI:en_LAI);
SWC_sLAI = SWC(st_LAI:en_LAI);
PAR_sLAI = PAR(st_LAI:en_LAI);
qc_sLAI = qc(st_LAI:en_LAI);
qc_Sc_sLAI = qc_Sc(st_LAI:en_LAI);
AGC_c_sLAI = AGC_c(st_LAI:en_LAI);
qc_h2o_flux_sLAI = qc_h2o_flux(st_LAI:en_LAI);
Precip_halfdayago_sLAI = Precip_halfdayago(st_LAI:en_LAI);
Precip_1daysago_sLAI = Precip_1daysago(st_LAI:en_LAI);
Precip_2daysago_sLAI = Precip_2daysago(st_LAI:en_LAI);


use_sLAI = find(Precip_2daysago_sLAI < 1 & Precip_1daysago_sLAI < 0.5 & Precip_halfdayago_sLAI < 0.2 & qc_h2o_flux_sLAI == 0 & AGC_c_sLAI == 0 ...
     & vpd_sLAI > 1 & vpd_sLAI < 1.5 & SWC_sLAI > 0.11 & PAR_sLAI > 800 & PAR_sLAI < 1200 & qc_sLAI == 0 & qc_Sc_sLAI == 0);

test_s_LAI_sLAI = s_LAI(use_sLAI);
test_gc_sLAI = gc_sLAI(use_sLAI);
test_GEP_sLAI = GEP_sLAI(use_sLAI);



figure; binplot_v2(test_gc_sLAI,test_GEP_sLAI,10,'k',10);
figure; binplot_v2(test_s_LAI_sLAI,test_GEP_sLAI,10,'k',10);



% 












