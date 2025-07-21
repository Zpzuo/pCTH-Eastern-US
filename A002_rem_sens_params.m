fgrp_cd = 500; % code of forest type group of interest
%ecor_cd = 1323; % code of a specific EPA ecoregion (Level 2)
fage_th = 90; % threshold of forest age of interest
Lai_th = 2.0; % threshold of LAI of interest (options: 0.0, 2.0, 4.0)
%%
addpath('/projectnb/amazondr/data10/cliveg/zpzuo/paper2_intermediate_data/')
ht = readgeoraster('Ht_byGediL2A_250m.tif');
Wue = readgeoraster('Wue_multiYearPeakSeasonAggr_by20yrPMLv2.tif'); % in gC/gH2O
Wue = Wue * 1E3; % in gC/kgH2O
Lai = readgeoraster('Lai_multiYearPeakSeasonAggr_bySnInd.tif');
fgroup = readgeoraster('Fgroup_byFIA_reproj5070_250m.tif');
fage = readgeoraster('Fage_bilinear_250m_byBGI.tif');
ecor = readgeoraster('ecoRegion_L2_raster_byEPA.tif');

%% Create masks
mask_fgroup = (fgroup==fgrp_cd);
mask_fage = (fage>=fage_th);
mask_ecor = (ecor==ecor_cd);
mask_all = (mask_fgroup & mask_fage);

%% Fit LAI (in m2[leaf]/m2[ground])
Lai_masked = Lai(mask_all & Lai>Lai_th);
histogram(Lai(mask_all), 200, 'Normalization','pdf')
hold on
plot(0:0.1:10, ...
     gampdf(0:0.1:10, ...
            mean(Lai_masked,'all','omitnan')^2/std(Lai_masked,0,'all','omitnan')^2, ...
            std(Lai_masked,0,'all','omitnan')^2/mean(Lai_masked,'all','omitnan')), ...
     'Color','red'); % gamma distribution fit
plot(0:0.1:10, ...
     normpdf(0:0.1:10, ...
             mean(Lai_masked,'all','omitnan'), ...
             std(Lai_masked,0,'all','omitnan')), ...
     "Color",'green'); % normal distribution fit
hold off
disp(['Forest group code: ', num2str(fgrp_cd)])
disp(['mean: ', num2str(mean(Lai_masked,'all','omitnan'))])
disp(['std: ', num2str(std(Lai_masked,0,'all','omitnan'))])
disp(['a: ', num2str(mean(Lai_masked,'all','omitnan')^2/std(Lai_masked,0,'all','omitnan')^2)])
disp(['b: ', num2str(std(Lai_masked,0,'all','omitnan')^2/mean(Lai_masked,'all','omitnan'))])
fprintf('\n')

%% Fit WUE (in gC/kgH2O)
Wue_masked = Wue(mask_all);
histogram(Wue_masked, 200, 'Normalization','pdf')
hold on
plot(0:0.1:10, ...
     gampdf(0:0.1:10, ...
            mean(Wue_masked,'all','omitnan')^2/std(Wue_masked,0,'all','omitnan')^2, ...
            std(Wue_masked,0,'all','omitnan')^2/mean(Wue_masked,'all','omitnan')), ...
     'Color','red'); % gamma distribution fit
plot(0:0.1:10, ...
     normpdf(0:0.1:10, ...
             mean(Wue_masked,'all','omitnan'), ...
             std(Wue_masked,0,'all','omitnan')), ...
    'Color','green'); % normal distribution fit
hold off
disp(['Forest group code: ', num2str(fgrp_cd)])
disp(['mean: ', num2str(mean(Wue_masked,'all','omitnan'))])
disp(['std: ', num2str(std(Wue_masked,0,'all','omitnan'))])
disp(['a: ', num2str(mean(Wue_masked,'all','omitnan')^2/std(Wue_masked,0,'all','omitnan')^2)])
disp(['b: ', num2str(std(Wue_masked,0,'all','omitnan')^2/mean(Wue_masked,'all','omitnan'))])
fprintf('\n')

%% Fit h_ref (in meters)
ht_masked = ht(mask_all);
ht_tall = ht_masked(ht_masked>15);
histogram(ht_masked, 200, 'Normalization','pdf')
hold on
plot(0:0.1:50, ...
     gampdf(0:0.1:50, ...
            mean(ht_tall,'all','omitnan')^2/std(ht_tall,0,'all','omitnan')^2, ...
            std(ht_tall,0,'all','omitnan')^2/mean(ht_tall,'all','omitnan')), ...
     'Color','red'); % gamma distribution fit
plot(0:0.1:50, ...
     normpdf(0:0.1:50, ...
             mean(ht_tall,'all','omitnan'), ...
             std(ht_tall,0,'all','omitnan')), ...
    'Color','green'); % normal distribution fit
hold off

disp(['mean: ', num2str(mean(ht_tall,'all','omitnan'))])
disp(['std: ', num2str(std(ht_tall,0,'all','omitnan'))])
disp(['a: ', num2str(mean(ht_tall,'all','omitnan')^2/std(ht_tall,0,'all','omitnan')^2)])
disp(['b: ', num2str(std(ht_tall,0,'all','omitnan')^2/mean(ht_tall,'all','omitnan'))])
fprintf('\n')