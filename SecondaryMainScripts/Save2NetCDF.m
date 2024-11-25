clear all

folder = '/net/kryo/work/ursho/Damiano_Presence_data/presence_absence_tables_ensemble_averages/Group_specific_background_approach';
addpath(genpath(folder))

cd(folder)
%% load data

load('/net/kryo/work/ursho/Damiano_Presence_data/presence_absence_tables_ensemble_averages/Group_specific_background_approach/Data/HelpVariables2.mat', 'latchl')
load('/net/kryo/work/ursho/Damiano_Presence_data/presence_absence_tables_ensemble_averages/Group_specific_background_approach/Data/HelpVariables2.mat', 'lonchl')

lats = latchl;
lons = lonchl;

load('/net/kryo/work/ursho/Damiano_Presence_data/presence_absence_tables_ensemble_averages/Group_specific_background_approach/Data/07Biomes/Seasonally_corrected_original.mat')
load('/net/kryo/work/ursho/Damiano_Presence_data/presence_absence_tables_ensemble_averages/Group_specific_background_approach/Data/07Biomes/Biomes_annual_monthly.mat')

I = [5     1     9     6     2     8     7     3     4];

annual_tmp = ones(1,180,360).*NaN;
season_tmp = corr_season_smooth;
month_tmp = ones(12,180,360).*NaN;


for i = 1:9
    
    annual_tmp(smooth_annual_map == I(i)) = i;
 
    month_tmp(smooth_map == I(i)) = i;
end

plotSOM(annual_tmp,1,8)
hold on;
cmap = morgenstemning(12);
[cmap_tmp] = shuffle_colormap(cmap);
[cmap_tmp1] = shuffle_colormap(cmap_tmp);


cmap = ametrine(12);
[cmap_tmp] = shuffle_colormap(cmap);
[cmap_tmp2] = shuffle_colormap(cmap_tmp);

cmap = isolum(12);
[cmap_tmp] = shuffle_colormap(cmap);
[cmap_tmp3] = shuffle_colormap(cmap_tmp);

comb_cmap = [cmap_tmp2([2,1,3],:);cmap_tmp1(7,:);cmap_tmp3(5,:);cmap_tmp1(8,:);cmap_tmp1(9,:);cmap_tmp1(4,:)];

colormap(comb_cmap(1:8,:))
c = colorbar;
set( c, 'YDir', 'reverse' );
[positions] = get_ticks_centered(8);
c.YTick = positions;
c.TickLabels = {'(1) TRP','(2) HIL','(3) WIS','(4) SUS','(5) HIT ', '(6) MTR',...
    '(7) PEU', '(8) SMN'};
set(findall(gcf,'-property','FontSize'),'FontSize',30)
set(findall(gcf,'-property','LineWidth'),'LineWidth',3)

annual_map = annual_tmp;
seasonal_map = season_tmp;
monthly_map = month_tmp;

% 
annual_map = ones(360,180,1).*NaN;
seasonal_map = ones(360,180,4).*NaN;
monthly_map = ones(360,180,12).*NaN;

annual_map(:,:,1) = squeeze(annual_tmp(1,:,:)).';

for t = 1:4
    seasonal_map(:,:,t) = squeeze(season_tmp(t,:,:)).';
end

for t = 1:12
    monthly_map(:,:,t) = squeeze(month_tmp(t,:,:)).';
end

unique(annual_map(~isnan(annual_map)))

unique(seasonal_map(~isnan(seasonal_map)))
unique(monthly_map(~isnan(monthly_map)))
%% Save the masks in a NetCDF file

cd('/net/kryo/work/ursho/Damiano_Presence_data/presence_absence_tables_ensemble_averages/Group_specific_background_approach/Data/07Biomes')

disp('Masks calculated, now saving netcdf...');

today=date; % get today's date
Season = ["Spring", "1"; "Summer", "2"; "Fall", "3"; "Winter",'4'];
Month = ["January",'1'; "February", '2'; "March", '3'; "April", '4';
    "May", "5"; "June", "6"; "July", "7"; "August", "8";
    "September", "9"; "October", "10"; "November", "11";
    "December", "12"];

% --------------------------
%     WRITE netcdf file!
% --------------------------
NetCDF_filename = 'j.pocean.2021.102530_Data_file.nc';

nccreate(NetCDF_filename,'lats','Dimensions',{'lat',180},'Format','netcdf4');
nccreate(NetCDF_filename,'lons','Dimensions',{'lon',360},'Format','netcdf4');

nccreate(NetCDF_filename,'annual_biomes','Dimensions',{'time',1,'lat',180,'lon',360},'Datatype','double','Format','netcdf4','FillValue',NaN);
nccreate(NetCDF_filename,'seasonal_biomes','Dimensions',{'season', 4,'lat',180,'lon',360},'Datatype','double','Format','netcdf4','FillValue',NaN);
nccreate(NetCDF_filename,'monthly_biomes','Dimensions',{'month', 12,'lat',180,'lon',360},'Datatype','double','Format','netcdf4','FillValue',NaN);

ncwrite(NetCDF_filename,'lons',lons);
ncwriteatt(NetCDF_filename, 'lons', 'standard_name', 'longitude');
ncwriteatt(NetCDF_filename, 'lons', 'units', 'degrees East');
ncwriteatt(NetCDF_filename, 'lons', 'axis', 'X');

ncwrite(NetCDF_filename,'lats',lats);
ncwriteatt(NetCDF_filename, 'lats', 'standard_name', 'latitude');
ncwriteatt(NetCDF_filename, 'lats', 'units', 'degrees North');
ncwriteatt(NetCDF_filename, 'lats', 'axis', 'Y');

ncwrite(NetCDF_filename,'annual_biomes',annual_map);
ncwriteatt(NetCDF_filename, 'annual_biomes', 'description', 'Annual biomes. The number indicates the label of the corresponding biome for each grid cell');
ncwriteatt(NetCDF_filename, 'annual_biomes', 'resolution', '1°latitude x 1°longitude');
ncwriteatt(NetCDF_filename, 'annual_biomes', 'units', '1');
ncwriteatt(NetCDF_filename, 'annual_biomes', 'range', '1 - 7');
ncwriteatt(NetCDF_filename, 'annual_biomes', 'value', '1: TRP, 2: HIL, 3: WIS, 4: SUS, 5: HIT, 6: MTR, 7: PEU');

ncwrite(NetCDF_filename,'seasonal_biomes',seasonal_map);
ncwriteatt(NetCDF_filename, 'seasonal_biomes', 'description', 'Seasonal biomes for global spring, summer, fall, and winter. The number indicates the label of the corresponding biome for each grid cell');
ncwriteatt(NetCDF_filename, 'seasonal_biomes', 'resolution', '1°latitude x 1°longitude');
ncwriteatt(NetCDF_filename, 'seasonal_biomes', 'units', '1');
ncwriteatt(NetCDF_filename, 'seasonal_biomes', 'range', '1 - 8');
ncwriteatt(NetCDF_filename, 'seasonal_biomes', 'value', '1: TRP, 2: HIL, 3: WIS, 4: SUS, 5: HIT, 6: MTR, 7: PEU, 8: SMN');

ncwrite(NetCDF_filename,'monthly_biomes',monthly_map);
ncwriteatt(NetCDF_filename, 'monthly_biomes', 'description', 'Monthly biomes. The number indicates the label of the corresponding biome for each grid cell');
ncwriteatt(NetCDF_filename, 'monthly_biomes', 'resolution', '1°latitude x 1°longitude');
ncwriteatt(NetCDF_filename, 'monthly_biomes', 'units', '1');
ncwriteatt(NetCDF_filename, 'monthly_biomes', 'range', '1 - 9');
ncwriteatt(NetCDF_filename, 'monthly_biomes', 'value', '1: TRP, 2: HIL, 3: WIS, 4: SUS, 5: HIT, 6: MTR, 7: PEU, 8: SMN, 9: -');

% Create Global Attributes
ncwriteatt(NetCDF_filename,'/','Title of dataset',...
    ['Phytoplankton based biome identity at sea surface at 1° latitude x 1° longitude resolution, for the annual, seasonal, and monthly scale' ]);
ncwriteatt(NetCDF_filename,'/','Title of associated publication',...
    ['Biome partitioning of the global ocean based on phytoplankton biogeography. DOI:10.1016/j.pocean.2021.102530' ]);
ncwriteatt(NetCDF_filename,'/','Authors',...
    ['Urs Hofmann Elizondo, Damiano Righetti, Fabio Benedetti, Meike Vogt' ]);
ncwriteatt(NetCDF_filename,'/','Journal',...
    ['Progress in Oceanography. Volume 194' ]);
ncwriteatt(NetCDF_filename,'/','Year',...
    ['2021' ]);
ncwriteatt(NetCDF_filename,'/','Date of creation',today);
ncwriteatt(NetCDF_filename,'/','Source of data','Urs Hofmann Elizondo');
ncwriteatt(NetCDF_filename,'/','Institution','ETH Zurich, Zurich Switzerland');
ncwriteatt(NetCDF_filename,'/','Department','Environmental Systems Science');
ncwriteatt(NetCDF_filename,'/','Contact','Urs Hofmann Elizondo (urs.hofmann@usys.ethz.ch)');
