fgrp_cd = 500; % code of forest type group of interest
fage_th = 90; % threshold of forest age of interest
ht_th = 15; % threshold of forest height
burnin = 500; % burn-in threshold for the MCMC chain

%rng("default")
num_samples = 5E3; % number of observation samples used in calibration

% Load packages
addpath(genpath("/usr4/ge646/zpzuo/tools/matlab/"))

%% Assemble masks
% try
%     addpath(genpath("C:\Users\zhenp\OneDrive - Boston University\2nd_paper\"))
%     addpath(genpath("/Users/zpzuo/Library/CloudStorage/Onedrive-BostonUniversity/2nd_paper/"))
%     addpath(genpath("/projectnb/amazondr/data10/cliveg/zpzuo/paper2_intermediate_data/"))
% catch
% end
inpath1 = "/projectnb/amazondr/data10/cliveg/zpzuo/paper2_intermediate_data/";
fgroup = readgeoraster(inpath1 + "Fgroup_byFIA_reproj5070_250m.tif");
fage = readgeoraster(inpath1 + "Fage_bilinear_250m_byBGI.tif");
ht_gedi = readgeoraster(inpath1 + "Ht_byGediL2A_250m.tif");
twi = readgeoraster(inpath1 + "Twi_standard_250m_conus.tif"); % topographic wetness index
ecor = readgeoraster(inpath1 + "ecoRegion_L2_raster_byEPA.tif");

fgrp_ecor = struct( ...
    "fgrp100", [1412,925], ...
    "fgrp120", [1556,982,925,964], ...
    "fgrp140", [2020], ...
    "fgrp160", [2158,2020,1598,1583,1604,2222,2224], ...
    "fgrp400", [2222,2224,2158,1598,1583,], ...
    "fgrp500", [1556,982,1552,968,969,2222,2224,1607,2158,1598,1583,1604,1585,1584,1580,1412,1469], ...
    "fgrp600", [297,2020,1598,1583], ... % 
    "fgrp800", [1556,982,1552,972,971,968,969,970,1607,1580,964,1469], ...
    "fgrp900", [964,1469,925,1412,1580]);

%% Generate a mask -- `num_samples` 1-valued cells from forest group of interest, forest age of interest, and valid GEDI observations
mask_fgroup = (fgroup==fgrp_cd);
mask_fage = (fage>=fage_th);
mask_twi = (twi>0);
mask_gedi = (ht_gedi>ht_th);
mask_ecor = ismember(ecor, getfield(fgrp_ecor, "fgrp"+num2str(fgrp_cd)));

num_available = sum(mask_fgroup&mask_fage&mask_twi&mask_gedi&mask_ecor, 'all');
mask = randsp(mask_fgroup&mask_fage&mask_twi&mask_gedi&mask_ecor, num_samples); mask = logical(mask);

%% Load the mask
% num_samples = 5E3; % number of observation samples used in calibration
% load(inpath1 + "mask5k_fgrp"+fgrp_cd+".mat")
% 
%% Load multiple year input data
tavg = NaN([num_samples,44]); 
prec = tavg; 
dayl = tavg;
par = tavg;
parfor i = 1:44
    year = i + 1979;
    % Read input data
    tavg_i = readgeoraster("tair_year"+year+"PeakSeasonAggr_byAgERA5.tif"); % daytime temperature (in C)
    prec_i = readgeoraster("prec_year"+year+"PeakSeasonAggr_byDAYMET.tif"); % precipitation (in mm/day)
    dayl_i = readgeoraster("dayl_year"+year+"PeakSeasonAggr_byDAYMET.tif"); % day length (in seconds)
    srad_i = readgeoraster("srad_year"+year+"PeakSeasonAggr_byDAYMET.tif"); % shortwave radiation (in W m-2)
    par_i = srad_i * 0.45; % photosynthetically active radiation, 45% of shortwave radiation
    tavg(:,i) = tavg_i(mask);
    prec(:,i) = prec_i(mask);
    dayl(:,i) = dayl_i(mask);
    par(:,i) = par_i(mask);

    tavg_i = []; prec_i=[]; dayl_i=[]; srad_i=[]; par_i=[];
end
input_masked = struct('temp', tavg, ...
    'prec', prec, ...
    'lday', double(dayl), ...
    'I0', double(par), ...
    'twi', double(twi(mask)), ...
    'distMat', distMat);
data.xdata = input_masked;

%% (optional) Load mean input data
% [tavg, R] = readgeoraster('tavg_multiYearPeakSeasonAggr_byRTMA.tif'); % daytime temperature (in C)
% prec = readgeoraster('prec_multiYearPeakSeasonAggr_byDAYMET.tif'); % precipitation (in mm/day)
% dayl = readgeoraster('dayl_multiYearPeakSeasonAggr_byDAYMET.tif'); % day length (in seconds)
% srad = readgeoraster('srad_multiYearPeakSeasonAggr_byDAYMET.tif'); % shortwave radiation (in W m-2)
% par = srad * 0.45; % photosynthetically active radiation, 45% of shortwave radiation
% 
% input_masked = struct('temp', tavg(mask), ...
%     'prec', prec(mask), ...
%     'lday', double(dayl(mask)), ...
%     'I0', double(par(mask)), ...
%     'twi', double(twi(mask)) );
% data.xdata = input_masked;

%% Load target output data
ho = ht_gedi(mask);
data.ydata = ho;


%% Derive prior
% Define stats for variable parameters
prior_table = readtable(inpath1 + "parameter_priors_normal.xlsx");
hrefPriorMu = prior_table{prior_table{:,1}==fgrp_cd, "h_ref_mu"};
hrefPriorSig = prior_table{prior_table{:,1}==fgrp_cd, "h_ref_sig"};
gPriorMu = prior_table{prior_table{:,1}==fgrp_cd, "gamma_mu"};
gPriorSig = prior_table{prior_table{:,1}==fgrp_cd, "gamma_sig"};
WuePriorMu = prior_table{prior_table{:,1}==fgrp_cd, "Wue_mu"};
WuePriorSig = prior_table{prior_table{:,1}==fgrp_cd, "Wue_sig"};
LaiPriorMu = prior_table{prior_table{:,1}==fgrp_cd, "Lai_mu"};
LaiPriorSig = prior_table{prior_table{:,1}==fgrp_cd, "Lai_sig"};
coeRcHtPriorMu = prior_table{prior_table{:,1}==fgrp_cd, "coe_rc_ht_mu"};
coeRcHtPriorSig = prior_table{prior_table{:,1}==fgrp_cd, "coe_rc_ht_sig"};
sigmaPriorMu = 6;
sigmaPriorSig = 2;

paramV1 = {
    {'href', 25, 0, 100, hrefPriorMu, hrefPriorSig}
    {'gamma', 0.2, 0, 1, gPriorMu, gPriorSig}
    {'WUE', 5, 0, 50, WuePriorMu, WuePriorSig}
    {'LAI', 5, 0, 10, LaiPriorMu, LaiPriorSig}
    {'coeRcHt', 0.1, 0, 1, coeRcHtPriorMu, coeRcHtPriorSig}
    {'sigma', 6, 0, 50, sigmaPriorMu, sigmaPriorSig}
  };    
paramV2 = {
    {'href',  hrefPriorMu+randn*hrefPriorSig, 0, 100, hrefPriorMu, hrefPriorSig}
    {'gamma', gPriorMu+randn*gPriorSig, 0, 1, gPriorMu, gPriorSig}
    {'WUE', WuePriorMu+randn*WuePriorSig, 0, 50, WuePriorMu, WuePriorSig}
    {'LAI', LaiPriorMu+randn*LaiPriorSig, 0, 10, LaiPriorMu, LaiPriorSig}
    {'coeRcHt', coeRcHtPriorMu+randn*coeRcHtPriorSig, 0, 1, coeRcHtPriorMu, coeRcHtPriorSig}
    {'sigma', sigmaPriorMu+randn*sigmaPriorSig, 0, 50, sigmaPriorMu, sigmaPriorSig}
  };

% Logarithmic prior & likelihood function
model.N = 1E4 + burin;
model.priorfun = @(paramV, mu, sig) -2*sum(-log(sig) - .5*log(2*pi) - .5*(((paramV - mu)./sig).^2));
%model.priorfun = @(paramV, mu, sig) 0;
model.ssfun = @(paramV, data) -2*uvnLike(paramV, data);

%% MCMC sampling
options.nsimu = 5E2;
options.updatesigma = 0;
options.waitbar = 1;
options.method = 'dram';
options.adaptint = 50;

[res,chain] = mcmcrun(model, data, paramV1, options);

%%
figure
Ym = mean(asrl(mean(chain,1), data.xdata), 2);
rd = (Ym - data.ydata) ./ data.ydata;
histogram(rd)

%%
param_mean = mean(chain, 1);
chain_end = chain(end,:);
hm_end = asrl(chain_end, input_masked);
disp(['Mean parameters: ', num2str(param_mean)])
disp(['Mean deviation: ', num2str(mean(hm_end-ho))])
disp(['RMSE: ', num2str(rmse(hm_end, ho))])

%%
figure
[~, ax] = plotmatrix(chain);
ylabel(ax(1,1), "h_{ref}")
ylabel(ax(2,1), "$\gamma$", 'Interpreter', 'latex', 'FontWeight','bold')
ylabel(ax(3,1), "WUE") 
ylabel(ax(4,1), "LAI")
ylabel(ax(5,1), "$\frac{\Delta h}{\Delta r_{c}}$", 'Interpreter','latex', 'FontWeight','bold')
ylabel(ax(6,1), "$\sigma$", 'Interpreter','latex' ,'FontWeight','bold')
xlabel(ax(6,1), "h_{ref}")
xlabel(ax(6,2), "$\gamma$", 'Interpreter', 'latex', 'FontWeight','bold')
xlabel(ax(6,3), "WUE") 
xlabel(ax(6,4), "LAI")
xlabel(ax(6,5), "$\frac{\Delta h}{\Delta r_{c}}$", 'Interpreter','latex', 'FontWeight','bold')
xlabel(ax(6,6), "$\sigma$", 'Interpreter','latex' ,'FontWeight','bold')

figure
ax = ecornerplot(chain, 'ks',true, 'color',[0.6,0.35,0.3]);
ylabel(ax(1,1), "h_{ref}")
ylabel(ax(2,1), "$\gamma$", 'Interpreter', 'latex', 'FontWeight','bold')
ylabel(ax(3,1), "WUE") 
ylabel(ax(4,1), "LAI")
ylabel(ax(5,1), "$\frac{\Delta h}{\Delta r_{c}}$", 'Interpreter','latex', 'FontWeight','bold')
ylabel(ax(6,1), "$\sigma$", 'Interpreter','latex' ,'FontWeight','bold')
xlabel(ax(6,1), "h_{ref}")
xlabel(ax(6,2), "$\gamma$", 'Interpreter', 'latex', 'FontWeight','bold')
xlabel(ax(6,3), "WUE") 
xlabel(ax(6,4), "LAI")
xlabel(ax(6,5), "$\frac{\Delta h}{\Delta r_{c}}$", 'Interpreter','latex', 'FontWeight','bold')
xlabel(ax(6,6), "$\sigma$", 'Interpreter','latex' ,'FontWeight','bold')

%% Functions
% Log prior for individual parameter
function logprior = uniformPrior(param, paramPriorMin, paramPriorMax)
    if param>=paramPriorMin && param<=paramPriorMax
        logprior = log(1/(paramPriorMax - paramPriorMin));
    else
        logprior = -Inf;
    end
end

function logprior = normalPrior(param, paramPriorMean, paramPriorSigma)
    logprior = log(normpdf(param, paramPriorMean, paramPriorSigma));
end

% Lod likelihood (univariate normal)
function loglike = uvnLike(paramV, data)
    % Extract necessary info
    errStd = paramV(6); % residual std at any pixel
    Yo = data.ydata; 
    % Derive modeled Y 
    Ym = asrl(paramV, data.xdata);

    if errStd < 0
        loglike = -Inf;
        return
    end
    loglike = sum(log(normpdf(Yo, mean(Ym,2), errStd)));
end

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

% Sample k 1-valued elements from a 0-1 matrix
function mask = randsp (matrix01, k)
    [nrow, ncol] = size(matrix01);
    % Find the indices of the 1-valued cells in the matrix
    [row_indices, col_indices] = find(matrix01 == 1);
    % Randomly sample specified number of indices
    random_indices = randperm(length(row_indices), k);
    % Extract the sampled indices
    sampled_row_indices = row_indices(random_indices);
    sampled_col_indices = col_indices(random_indices);
    % Convert row and column indices to 1-dimensional indices
    sampled_linear_indices = sub2ind(size(matrix01), sampled_row_indices, sampled_col_indices);
    % Create a new matrix with sampled cells set to 1 and others to 0
    mask = zeros(nrow, ncol);
    mask(sampled_linear_indices) = 1;
end

% Sigmoid function
function x = Sigmoid(num)
    x = 1./ (1 + exp(-num));
end