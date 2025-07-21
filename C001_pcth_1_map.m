% This code runs the canopy height prediction model with one parameter
% sample. The resulting tree height map will be saved as a geotiff image.
function A001_run_asrl_1_param_44clim (fgrp_cd, sample_cd)
%outpath = "/Users/zpzuo/OneDrive - Boston University/2nd_paper/paper2_intermediate_data/";
outpath = "/projectnb/amazondr/data10/cliveg/zpzuo/paper2_intermediate_data/ImagePerTrait_1000samples/fgrp"+fgrp_cd+"/";

% Add path of supportive packages
addpath(genpath("C:\Users\zhenp\OneDrive - Boston University\tools\matlab\")) % Windows
addpath(genpath("/Users/zpzuo/Library/CloudStorage/Onedrive-BostonUniversity/tools/matlab/")) % MacOS 1
addpath(genpath("/Users/zpzuo/OneDrive - Boston University/tools/matlab/")) % MacOS 2
addpath(genpath("/usr4/ge646/zpzuo/tools/matlab/")) % Linux - SCC

% Add path of input data
addpath(genpath("C:\Users\zhenp\OneDrive - Boston University\2nd_paper\")) % Windows
addpath(genpath("/Users/zpzuo/Library/CloudStorage/Onedrive-BostonUniversity/2nd_paper/")) % MacOS 1
addpath(genpath("/Users/zpzuo/OneDrive - Boston University/2nd_paper/")) % MacOS 2
addpath(genpath("/projectnb/amazondr/data10/cliveg/zpzuo/paper2_intermediate_data/")) % Linux - SCC

% Load parameter sample
samples = load("fgrp"+fgrp_cd+"_1000samples.mat");
paramV = samples.paramV(sample_cd,:);

% Retrieve ecoregions of occurrence
% fgrp_ecor = struct( ...
%     "fgrp100", [1412,925], ...
%     "fgrp120", [1556,982,925,964], ...
%     "fgrp140", [2020], ...
%     "fgrp160", [2158,2020,1598,1583,1604,2222,2224], ...
%     "fgrp400", [2222,2224,2158,1598,1583,], ...
%     "fgrp500", [1556,982,1552,968,969,2222,2224,1607,2158,1598,1583,1604,1585,1584,1580,1412,1469], ...
%     "fgrp600", [297,2020,1598,1583], ...
%     "fgrp800", [1556,982,1552,972,971,968,969,970,1607,1580,964,1469], ...
%     "fgrp900", [964,1469,925,1412,1580]);
fgrp_ecor = struct( ...
    "fgrp100", [1412,925, ...
                1556,982,1552,968,969,2222,2224,1607,2158,1598,1583,1604,1585,1584,1580,1469,972,971,970,1351,2197,2020,2399,964,1600,1599,1581,1606,1605,1603,1602,1601,967,2228,1582], ...
    "fgrp120", [1556,982,925,964, ...
                1552,968,969,2222,2224,1607,2158,1598,1583,1604,1585,1584,1580,1412,1469,972,971,970,1351,2197,2020,2399,1600,1599,1581,1606,1605,1603,1602,1601,967,2228,1582], ...
    "fgrp140", [2020, ...
                1556,982,925,964,1552,968,969,2222,2224,1607,2158,1598,1583,1604,1585,1584,1580,1412,1469,972,971,970,1351,2197,2399,1600,1599,1581,1606,1605,1603,1602,1601,967,2228,1582], ...
    "fgrp160", [2158,2020,1598,1583,1604,2222,2224, ...
                1556,982,925,964,1552,968,969,1607,1585,1584,1580,1412,1469,972,971,970,1351,2197,2399,1600,1599,1581,1606,1605,1603,1602,1601,967,2228,1582], ...
    "fgrp400", [2222,2224,2158,1598,1583, ...
                2020,1556,982,925,964,1552,968,969,1607,1604,1585,1584,1580,1412,1469,972,971,970,1351,2197,2399,1600,1599,1581,1606,1605,1603,1602,1601,967,2228,1582], ...
    "fgrp500", [1556,982,1552,968,969,2222,2224,1607,2158,1598,1583,1604,1585,1584,1580,1412,1469, ...
                972,971,970,1351,2197,2020,2399,964,925,1600,1599,1581,1606,1605,1603,1602,1601,967,2228,1582], ...
    "fgrp600", [2020,1598,1583, ...
                1556,982,925,964,1552,968,969,2222,2224,1607,2158,1604,1585,1584,1580,1412,1469,972,971,970,1351,2197,2399,1600,1599,1581,1606,1605,1603,1602,1601,967,2228,1582], ...
    "fgrp800", [1556,982,1552,972,971,968,969,970,1607,1580,964,1469, ...
                2020,925,2222,2224,2158,1598,1583,1604,1585,1584,1412,1351,2197,2399,1600,1599,1581,1606,1605,1603,1602,1601,967,2228,1582], ...
    "fgrp900", [964,1469,925,1412,1580, ...
                2020,1556,982,1552,968,969,2222,2224,1607,2158,1598,1583,1604,1585,1584,972,971,970,1351,2197,2399,1600,1599,1581,1606,1605,1603,1602,1601,967,2228,1582]);
ecor = readgeoraster("ecoRegion_L2_raster_byEPA.tif");
mask_ecor = ismember(ecor, getfield(fgrp_ecor, "fgrp"+num2str(fgrp_cd)));
mask = mask_ecor;

% Load data 
year_list = repmat(1980:2023,1,23); year_list = year_list(end:-1:end-999);
year = year_list(sample_cd); 
[tavg, R] = readgeoraster("tair_year"+year+"PeakSeasonAggr_byAgERA5.tif"); % daytime temperature (in C)
prec = readgeoraster("prec_year"+year+"PeakSeasonAggr_byDAYMET.tif"); % precipitation (in mm/day)
dayl = readgeoraster("dayl_year"+year+"PeakSeasonAggr_byDAYMET.tif"); % day length (in seconds)
srad = readgeoraster("srad_year"+year+"PeakSeasonAggr_byDAYMET.tif"); % shortwave radiation (in W m-2)
par = srad * 0.45; % photosynthetically active radiation, 45% of shortwave radiation
twi = readgeoraster('Twi_standard_250m_conus.tif'); % topographic wetness index
input = struct('temp', tavg(mask), ...
    'prec', prec(mask), ...
    'lday', double(dayl(mask)), ...
    'I0', double(par(mask)), ...
    'twi', double(twi(mask)) );
clear tavg prec dayl par twi srad

% Run model
hm = asrl(paramV, input);

% Write output as geotiff image
hm_map = NaN([11613, 18458]);
hm_map(mask) = hm;
clear hm
geotiffwrite(outpath + ...
        "Ht_byAsrl_"+"fgrp"+fgrp_cd+"_McMC_sample"+num2str(sample_cd)+".tif", ...
    hm_map, ...
    R, "CoordRefSysCode","EPSG:5070", ...
    "TiffTags", struct("Compression", "LZW"))

end

%%
% The model for potential tree height predictions
function ht_asrl = asrl(paramV, input)
    paramS = struct( ...
        'pow_rc_ht', 1, ...
        'T_cold', -10, ...
        'T_CO2_l', 0, ...
        'T_hot', 30, ...
        'T_CO2_h', 38, ...
        'pmax', 1E-3, ...
        'alpha', 1E-5, ...
        'k', 0.5, ...
        'm', 0.1);
    ht_candidates = [2,6,10,14,18,22,26,30,34,38,42,50,58,66,82,98,114,130,146];
    ncand = length(ht_candidates);
    for i = 1:ncand
        ht = ht_candidates(i);
        Qavail(:,:,i) = calc_rootabs(ht, paramV, paramS, input);
        Qtrans(:,:,i) = calc_transpr(ht, paramV, paramS, input);
    end
    ht_asrl = find_mat_linear_intr(ht_candidates, Qavail, Qtrans);
    ht_asrl(isnan(ht_asrl)) = 0;
end

% Root water absorption module
function Qavail = calc_rootabs(ht, paramV, paramS, input)
    h_ref = paramV(1); % reference tree height in meters
    gamma = paramV(2);
    coe_rc_ht = paramV(5); % Rr (root radius) = Rc (canopy radius) = coe_rc_ht * ht ^ pow_rc_ht
    pow_rc_ht = paramS.pow_rc_ht;
    mm_to_m = 1E-3;
    prec = input.prec * mm_to_m; % Growing season average precipitation in m/day, converted from mm/day
    Twi = input.twi; % Topographic wetness index, unitless
    Rr = 0.2 * ht.^pow_rc_ht; % Root system radius in meters
    Qavail = gamma.*pi.*Rr.^2.*prec.*Twi./16.*h_ref./ht; % Max water flow rate (growing season average) in m3(H2O)/tree/day
    %Qavail = gamma.*pi.*Rr.^2.*prec.*Sigmoid(Twi/4).*h_ref./ht; % Max water flow rate (growing season average) in m3(H2O)/tree/day
    %Qavail = gamma.*pi.*Rr.^2.*prec.*Twi./16.*h_ref./ht;
end

% Transpiration module
function Qtrans = calc_transpr(ht, paramV, paramS, input)
    Wue = paramV(3); % (current) water-use efficiency in g(C)/kg(H2O)
    Gpp = calc_gpp(ht, paramV, paramS, input); % GPP in g(C)/tree/day
    Qtrans = Gpp / Wue; % transpiration rate (growing season average) in kg(H2O)/tree/day
    Qtrans = Qtrans * 1E-3; % transpiration rate (growing season average) in m3(H2O)/tree/day
end

% Sub-function: Calculation of GPP
function Gpp = calc_gpp(ht, paramV, paramS, input)
    temp = input.temp; % air temperature in degree C
    I0 = input.I0; % PAR above canopy in J(PAR)/m2(ground)/s (averaged from sunrise to sunset during the vegetation period)
    lday = input.lday; % day length in seconds

    Lai = paramV(4); % leaf area index unitless
    coe_rc_ht = paramV(5); % rc = coe_rc_ht * ht ^ pow_rc_ht
    pow_rc_ht = paramS.pow_rc_ht; 
    pmax = paramS.pmax; % maximum leaf gross photosynthesis rate in g(CO2)/m2(leaf area)/s
    alpha = paramS.alpha; % leaf-level light use efficiency in g(CO2)/J(PAR)
    k = paramS.k; % light extinction coefficient unitless
    m = paramS.m; % leaf transmittance unitless
    T_cold = paramS.T_cold; % mean temperature of the coldest month the plant can cope with
    T_CO2_l = paramS.T_CO2_l; % the lowest temperature limit for CO2 assimilation
    T_hot = paramS.T_hot; % mean temperature of the hottest month the plant can cope with
    T_CO2_h = paramS.T_CO2_h; % the highest temperature limit for CO2 assimilation

    rc = coe_rc_ht * ht.^pow_rc_ht; % radius of crown in m
    k0 = 2*log(0.01/0.99)/(T_CO2_l - T_cold);
    k1 = 0.5*(T_CO2_l+T_cold);
    k2 = log(0.99/0.01)/(T_CO2_h-T_hot);

    phi_T_l = 1./(1 + exp(k0*k1 - temp));
    phi_T_h = 1 - 0.01*exp(k2.*(temp - T_hot));
    phi_T = phi_T_l .* phi_T_h;
    phi_T(phi_T<0) = 0; % inhibition factor for GPP due to low/high temperature

    P = pmax/k * log((alpha*k*I0 + pmax*(1-m))./(alpha*k*I0*exp(-k*Lai) + pmax*(1-m))); % instantaneous photosynthesis in g(CO2)/m2(ground)/s
    P = P * 0.4*12/44 * pi*rc.^2 .* lday; % instantaneous photosynthesis in g(C)/tree/day
    Gpp = P .* phi_T; % gross primary production in g(C)/tree/day
end

% Intersection searching module
function intr = find_mat_linear_intr (X, Y1, Y2) %
    nintervals = length(X)-1;
    intersection_layers = NaN(size(Y1));
    for i = 1:nintervals
        x_i = X(i);
        y1_i = Y1(:,:,i);
        y2_i = Y2(:,:,i);
        x_ip1 = X(i+1);
        y1_ip1 = Y1(:,:,i+1);
        y2_ip1 = Y2(:,:,i+1);

        intersection_i = ((y2_i*x_ip1 - y2_ip1*x_i) - (y1_i*x_ip1 - y1_ip1*x_i))./(y2_i - y2_ip1 + y1_ip1 - y1_i);
        intersection_i(intersection_i < x_i) = NaN;
        intersection_i(intersection_i > x_ip1) = NaN;
        intersection_layers(:,:,i) = intersection_i;
    end
    intr = min(intersection_layers,[],3);
end