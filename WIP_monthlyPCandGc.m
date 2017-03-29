PC = nan(12*3,1);
date_PC = datetime;
i = 1;
for y = 1:3
    for m = 1:12
        use = find(Year_t == 2013+y & Month_t == m & qc == 0 & qc_Sc == 0 & AGC_c == 0 ...
    & PAR > 800 & PAR < 1200 & vpd > 1 & vpd < 2);
        PC(i) = median(GEP(use))/1000;
        date_PC(i) = datetime(2013+y,m,15);
        gc_m(i) = median(gc(use));
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
    & PAR > 800 & PAR < 1200 & vpd > 1 & vpd < 2);
        date_gcm(i) = datetime(2013+y,m,15);
        gc_m(i) = median(gc(use));
        i = i+1;
    end
    i = i+1;
end

figure; plot(DateTime_LAI,LAI); hold on; plot(date_PC,PC*100);
plot(date_PC,gc_m*250);

