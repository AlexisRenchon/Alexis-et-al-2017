PC = nan(12*3,1);
date_PC = datetime;
i = 1;
for y = 1:3
    for m = 1:12
        use = find(Year_t == 2013+y & Month_t == m & qc == 0 & qc_Sc == 0 & AGC_c == 0 ...
    & PAR > 800 & PAR < 1200 & vpd > 1 & vpd < 1.5 & SWC > 0.11);
        PC(i) = median(GEP(use))/1000;
        date_PC(i) = datetime(2013+y,m,15);
        i = i+1;
    end
    i = i+1;
end

gc_m = nan(12*3,1);
date_gcm = datetime;
i = 1;
for y = 1:3
    for m = 1:12
        use = find(Year_t == 2013+y & Month_t == m & Precip_2daysago < 1 & Precip_1daysago < 0.5 & Precip_halfdayago < 0.2 & qc_h2o_flux == 0 & AGC_c == 0 ...
     & vpd > 1 & vpd < 1.5 & SWC > 0.11);
        date_gcm(i) = datetime(2013+y,m,15);
        gc_m(i) = median(gc(use));
        i = i+1;
    end
    i = i+1;
end

figure; plot(DateTime_LAI,LAI); hold on; plot(date_PC,PC*100);
plot(date_gcm,gc_m*250);

use_gc = find(Precip_2daysago < 1 & Precip_1daysago < 0.5 & Precip_halfdayago < 0.2 & qc_h2o_flux == 0 & AGC_c == 0 ...
     & vpd > 1 & vpd < 1.5 & SWC > 0.11);
use_pc = find(PAR > 800 & PAR < 1200 & vpd > 1 & vpd < 1.5 & SWC > 0.11);

figure; plot(DateTime_LAI,LAI,'LineStyle','none','Marker','o','MarkerEdgeColor','none','MarkerFaceColor','k','MarkerSize',4); hold on;
plot(DateTime_CUP(use_pc),(GEP(use_pc)/1000)*100,'LineStyle','none','Marker','o','MarkerEdgeColor','none','MarkerFaceColor','r','MarkerSize',4);
plot(DateTime_CUP(use_gc),gc(use_gc)*250,'LineStyle','none','Marker','o','MarkerEdgeColor','none','MarkerFaceColor','b','MarkerSize',4);


% TO DO 
% GET DAILY (or weekly?) PC and DAILY GC instead of half hourly, as median over a day
% (48 half-hour), get the value if at least 5 or 6 (?) data points
































