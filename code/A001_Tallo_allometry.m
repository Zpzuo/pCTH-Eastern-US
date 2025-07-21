addpath(genpath("/usr4/ge646/zpzuo/tools/matlab/"));

%% Load the Tallo database
Tallo = readtable("/projectnb/amazondr/data10/cliveg/zpzuo/Tallo_1.0.0/Tallo.csv", ...
    "TreatAsMissing", "NA");
[species_count, species_names] = groupcounts(table2array(Tallo(:, "species")));

%% Get species names
species = struct( ...
    "ftg_100", ["Pinus banksiana","Pinus resinosa","Pinus strobus","Tsuga canadensis"], ...
    "ftg_120", ["Abies balsamea","Picea glauca","Picea rubens","Picea mariana","Larix laricina","Thuja occidentalis","Abies fraseri"], ...
    "ftg_140", ["Pinus palustris","Pinus elliottii"], ...
    "ftg_160", ["Pinus taeda","Pinus echinata","Pinus virginiana","Pinus clausa","Pinus pungens","Pinus serotina","Pinus rigida","Pinus glabra"], ...
    "ftg_400", ["Pinus strobus","Quercus falcata","Quercus rubra","Fraxinus americana","Juniperus virginiana","Pinus palustris","Pinus echinata","Pinus virginiana","Pinus taeda","Pinus elliottii"], ...
    "ftg_500", ["Quercus stellata","Quercus marilandica","Quercus montana","Quercus alba","Quercus falcata","Quercus rubra","Carya aquatica","Carya cordiformis","Carya glabra","Carya laciniosa","Carya ovata","Carya texana","Carya tomentosa","Liriodendron tulipifera","Sassafras albidum","Diospyros virginiana","Quercus macrocarpa","Quercus coccinea","Juglans nigra","Robinia psuedoacacia","Quercus ilicifolia","Quercus velutina"], ...
    "ftg_600", ["Quercus michauxii","Quercus falcata","Liquidambar styraciflua","Quercus texana","Quercus phellos","Quercus lyrata","Carya aquatica","Chamaecyparis thyoides","Taxodium distichum","Nyssa aquatica","Magnolia virginiana","Nyssa biflora","Taxodium distichum"], ...
    "ftg_800", ["Acer saccharum","Fagus grandifolia","Betula alleghaniensis","Prunus serotina","Prunus pensylvanica","Fraxinus americana","Fraxinus nigra","Fraxinus pennsylvanica","Liriodendron tulipifera","Tilia americana","Ulmus americana","Ulmus rubra","Robinia pseudoacacia","Acer rubrum"], ...
    "ftg_900", ["Populus grandidentata","Populus tremuloides","Betula papyrifera","Betula populifolia","Populus balsamifera","Prunus pensylvanica"]);
%% Fit allometric relationship for one forest group
% Determine the species of interest
sp_cd = "ftg_900";
sp = getfield(species, sp_cd);
sp_idx = ismember(table2array(Tallo(:,"species")), sp); % table rows of the species
% Extract variables
ht_grp = table2array(Tallo(sp_idx, "height_m")); % measured heights in meters of the species individuals
rc_grp = table2array(Tallo(sp_idx, "crown_radius_m")); % measured crown radius in meters of the species individuals
ds_grp = table2array(Tallo(sp_idx, "stem_diameter_cm")) / 100; % measured stem diameter in meters of the species individuals

%%
allom_rc_ht = extractHtRcAllometryFun(ht_grp, rc_grp);
powRcHt_bar = allom_rc_ht(1);
coeRcHt_bar = exp(allom_rc_ht(2));
coe_rc_ht = rc_grp./ht_grp.^(powRcHt_bar);
figure
histogram(coe_rc_ht, 200, 'Normalization','pdf')
hold on
plot(-1:0.01:1, ...
     gampdf(-1:0.01:1, ...
            mean(coe_rc_ht,'all','omitnan')^2/std(coe_rc_ht,0,'all','omitnan')^2, ...
            std(coe_rc_ht,0,'all','omitnan')^2/mean(coe_rc_ht,'all','omitnan')), ...
     'Color','red'); % gamma distribution fit
plot(-1:0.01:1, ...
     normpdf(-1:0.01:1, ...
             mean(coe_rc_ht,'all','omitnan'), ...
             std(coe_rc_ht,0,'all','omitnan')), ...
    'Color','green'); % normal distribution fit
hold off
disp(['Forest group: ', num2str(sp_cd)])
disp(['mean: ', num2str(mean(coe_rc_ht,'all','omitnan'))])
disp(['std: ', num2str(std(coe_rc_ht,0,'all','omitnan'))])
disp(['a: ', num2str(mean(coe_rc_ht,'all','omitnan')^2/std(coe_rc_ht,0,'all','omitnan')^2)])
disp(['b: ', num2str(std(coe_rc_ht,0,'all','omitnan')^2/mean(coe_rc_ht,'all','omitnan'))])
fprintf('\n')
%%
allom_ds_ht = extractHtRsAllometryFun(ht_grp, ds_grp);
powRsHt_bar = allom_ds_ht(1);
coeRsHt_bar = exp(allom_ds_ht(2));
%% Function that handles the regression
function b = extractHtRcAllometryFun(ht, rc) 
    [nrow, ~] = size(ht);
    valid_idx = ht>0 & rc>0;
    sample_idx = rand([nrow, 1]) < 1;
    ht = ht(valid_idx & sample_idx);
    rc = rc(valid_idx & sample_idx);
    % Get the logrithmics
    ht_log = log(ht);
    rc_log = log(rc);
    % Filter the records by a kernel density threshold
    kd = ksdensity([ht_log,rc_log], [ht_log,rc_log]);
    kd_idx = kd >= quantile(kd, 0.0);
    ht_log = ht_log(kd_idx);
    rc_log = rc_log(kd_idx);
    % Run the regression: X -- ln(Ht); Y -- ln(Rc); y = b2 + b1*x
    b = [1.0, mean(rc_log-1.0*ht_log,'omitnan')];
    %b = polyfit(ht_log, rc_log, 1);
    rc_fit = polyval(b, ht_log);
    figure
    scatter_kde(ht_log, rc_log,'filled', 'MarkerSize', 10)
    hold on
    plot(ht_log, rc_fit, 'Color', 'r')
    xlabel('ln(tree height meters)')
    ylabel('ln(crown radius meters)')
    grid on
end

function b = extractHtRsAllometryFun(ht, ds) 
    [nrow, ~] = size(ht);
    valid_idx = ht>0 & ds>0;
    sample_idx = rand([nrow, 1]) < 1;
    ht = ht(valid_idx & sample_idx);
    ds = ds(valid_idx & sample_idx);
    % Get the logrithmics
    ht_log = log(ht);
    ds_log = log(ds);
    % Filter the records by a kernel density threshold
    kd = ksdensity([ht_log,ds_log], [ht_log,ds_log]);
    kd_idx = kd >= quantile(kd, 0.0);
    ht_log = ht_log(kd_idx);
    ds_log = ds_log(kd_idx);
    % Run the regression: X -- ln(Ht); Y -- ln(Ds); y = b2 + b1*x
    b = [3/2, mean(ds_log-3/2*ht_log,'omitnan')];
    %b = polyfit(ht_log, ds_log, 1); 
    ds_fit = polyval(b, 0:0.1:max(ht_log));
    figure
    scatter_kde(ht_log, ds_log,'filled', 'MarkerSize', 10)
    hold on
    plot(0:0.1:max(ht_log), ds_fit, 'Color', 'r')
    xlabel('ln(tree height meters)')
    ylabel('ln(stem diameter meters)')
    hold off
    grid on
    figure 
    scatter_kde(exp(ht_log), exp(ds_log),'filled', 'MarkerSize', 10)
    hold on
    plot(exp(0:0.1:max(ht_log)), exp(ds_fit), 'Color', 'r')
    xlabel('tree height meters')
    ylabel('stem diameter meters')
    hold off
    grid on
end
