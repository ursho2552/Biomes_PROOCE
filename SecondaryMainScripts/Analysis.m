%% Analysis

% =========================================================================
% After producing biomes on a monthly, seasonal, and annual scale, we now
% analyse their location and spatio-temporal evolution, their species
% composition, differences in between and whithin biome-composition,
% core, common or indicator speces found in biomes, the species networks,
% and the environmental conditions of each biome
% =========================================================================

% =========================================================================
% This script is a stand-alone script that uses the monthly, seasonal,
% annual biomes, and sometimes the raw presence/absence projections. Thus
% we first make sure that our workspace is empty, and everything is setup
% correctly.
% =========================================================================

%% Restart workspace

clear all
folder_main = '/net/kryo/work/ursho/PhD/Projects/Biomes/Scripts/Biomes_PROOCE';
addpath(genpath(folder_main))
cd(folder_main)

% =========================================================================
% Load all variables needed
% =========================================================================

%load trained SOM
cd('/net/kryo/work/ursho/PhD/Projects/Biomes/Scripts/Biomes_PROOCE/Data/01NeuronsError/')
load('Single_run_11.mat')

%load help variables
cd('/net/kryo/work/ursho/PhD/Projects/Biomes/Scripts/Biomes_PROOCE/Data/')
load('HelpVariables.mat')
load('Area_map.mat')

%construct or load simplified version of raw data
cd('/net/kryo/work/ursho/PhD/Projects/Biomes/Scripts/Biomes_PROOCE/Data/00Probabilities/')
load('Simple_sort_Data.mat')
load('Transformed_CompleteSuitePhyto.mat')
load('Seasonally_corrected_data.mat')
load('Names_species')

%load biomes
cd('/net/kryo/work/ursho/PhD/Projects/Biomes/Scripts/Biomes_PROOCE/Data/05Biomes/')
load('No_mean_PCA_biomes_9_v2_5_perc.mat')
load('No_mean_PCA_biomes_annual_9_v2_5_perc.mat')
load('No_mean_PCA_biomes_seasonal_9_v2_5_perc.mat')




[ID_maps] = prepare2plot( [No_nan_phyto_simple(:,2:4),No_nan_phyto_simple(:,1)]);
corrected_monthly_raw = NaN(12,180, 360);
corrected_monthly_ID = corrected_monthly_raw;
%shift Southern Hemisphere by 6 months
for i =1:12
    if(i<7)
        j = mod(i+6,13);
    else
        j = mod(i+6,13) +1;
    end
    %i and j are the indices of thmatlab e months that need to be combined
    corrected_monthly_raw(i,91:end,:) = raw_monthly_maps(i,91:end,:);
    corrected_monthly_raw(i,1:90,:) = raw_monthly_maps(j,1:90,:);

    corrected_monthly_ID(i,91:end,:) = ID_maps(i,91:end,:);
    corrected_monthly_ID(i,1:90,:) = ID_maps(j,1:90,:);

end



%% Sort biomes by size

% =========================================================================
% Change label numbering using descending mean area; UHE 09/08/2019
% =========================================================================

%Calculate area of biomes on a monthly basis
n_clusters = 9;
area_biome = NaN(n_clusters,12);
tmp_map = corrected_monthly_smooth;
for m = 1:12
    for i = 1:n_clusters
        tot_area = sum(area_map(~isnan(tmp_map(m,:,:))));
        area_biome(i,m) = sum(area_map(tmp_map(m,:,:) == i))/tot_area;
    end
end
area_biome(area_biome == 0) = NaN;

%Calculate mean area across all months and sort in descending order
[B, I] = sort(mean(area_biome,2,'omitnan'),'descend');

%print results
disp('Biomes sorted in descending order of mean area coverage:')
area_biome_sorted = area_biome(I,:)'

%allocate new matrices to store "sorted" biomes
corr_smooth_annual_map = smooth_annual_map;
corr_season_smooth = NaN(4,180,360);
corr_corrected_monthly_smooth = NaN(12,180,360);
corr_corrected_monthly_raw = NaN(12,180,360);
corr_monthly_smooth = NaN(12,180,360);
for i = 1:n_clusters
    corr_monthly_smooth(smooth_map == I(i)) = i;
    corr_smooth_annual_map(smooth_annual_map == I(i)) = i;
    corr_season_smooth(season_map_smooth == I(i)) = i;
    corr_corrected_monthly_smooth(corrected_monthly_smooth == I(i)) = i;
    corr_corrected_monthly_raw(corrected_monthly_raw == I(i)) = i;
end


%safe biomes with corrected labels
cd('/net/kryo/work/ursho/PhD/Projects/Biomes/Scripts/Biomes_PROOCE/Data/05Biomes/')

%save the dataset
if isfile('Seasonally_corrected_original.mat')
    disp('File already exists!')
else
    save('Seasonally_corrected_original','corr_smooth_annual_map','corr_season_smooth','corr_corrected_monthly_smooth','corr_monthly_smooth')
end

cd(folder_main)

%% Data for Table 3 in manuscript

% =========================================================================
% Calculate min, max and mean biome area from monthly area
% =========================================================================

area_biome(isnan(area_biome)) = 0
max_area = max(area_biome,[],2)
min_area = min(area_biome,[],2)

%print max, min, and mean area
disp('Max, min, and mean area of biomes:')
round(1000*max_area(I))/10
round(1000*min_area(I))/10
round(1000*B)/10

%% Figure 2 Annual biomes
% =========================================================================
% Plot annual biomes; UHE 09/08/2019
% =========================================================================

plotSOM(corr_smooth_annual_map,1,8)

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

%% Figure 3 Seasonal biomes

% =========================================================================
% Plot seasonal biomes; UHE 09/08/2019
% =========================================================================

%print unique labels per season to make sure we know which biomes are
%present
for i = 1:4
    unique(corr_season_smooth(i,~isnan(corr_season_smooth(i,:,:))))
end

%plot each season
for s = 1:4
    plotSOM(corr_season_smooth,s,8)
    hold on;
    colormap(comb_cmap)
    set(findall(gcf,'-property','FontSize'),'FontSize',30)
    set(findall(gcf,'-property','LineWidth'),'LineWidth',3)
    hold off
end

%% Table Annual and seasonal area (not used in final version)

% =========================================================================
% Get area of annual and seasonal biomes; UHE 09/08/2019
% =========================================================================

area_seasonal = NaN(9,4);
for s = 1:5
    for i = 1:n_clusters
        if(s<5)
            tot_area = sum(area_map(~isnan(corr_season_smooth(s,:,:))));
            area_seasonal(i,s) = sum(area_map(corr_season_smooth(s,:,:) == i))/tot_area;
        else
            tot_area = sum(area_map(~isnan(corr_smooth_annual_map(1,:,:))));
            area_seasonal(i,s) = sum(area_map(corr_smooth_annual_map(1,:,:) == i))/tot_area;

        end
    end
end

disp('Seasonal area coverage by biomes:')
area_seasonal  = round(1000*area_seasonal)/10

%% Figure A.15 Plot differences between seasonal biome location and annual location


diff_maps = NaN(4,180,360);
area_diff = NaN(4,1);
for s = 1:4

    for i = 1:180
        for j = 1:360
            if(corr_smooth_annual_map(1,i,j) == corr_season_smooth(s,i,j))
                diff_maps(s,i,j) = corr_smooth_annual_map(1,i,j);
            else
                diff_maps(s,i,j) = NaN;
            end
        end
    end

    diff_maps(s,isnan(corr_season_smooth(s,:,:))) = NaN;
    area_diff(s) = sum(area_map(isnan(diff_maps(s,:,:))))/sum(area_map(~isnan(corr_season_smooth(s,:,:))));

end


for i = 1:4

    plotSOM(diff_maps,i,8)
    cmap = parula(9);
    [cmap_tmp] = shuffle_colormap(cmap);
    [cmap_tmp] = shuffle_colormap(cmap_tmp);

    colormap(comb_cmap)

    set(findall(gcf,'-property','FontSize'),'FontSize',30)
    set(findall(gcf,'-property','LineWidth'),'LineWidth',3)

end


%% Find average latitude of biomes

% =========================================================================
% The analysis of the area showed that one biome is absent on seasonal
% scales and is overall the smallest one. Thus, we do not consider this
% biome in our analysis. This means that we only use 8 biomes
% =========================================================================


% =========================================================================
% Now we try to find the median location of each biome, i.e. their centroid
% =========================================================================

n_clusters = 8;
abs_lats = abs(lats);
median_lat = NaN(12,n_clusters);
for m = 1:12
    for i = 1:n_clusters
        [r, c] = find(squeeze(corr_corrected_monthly_smooth(m,:,:)) == i);
        if(~isempty(r))
            %get for each c the median of r, then the median of r
            c_un = unique(c);
            r_all = c_un.*NaN;
            tmp = NaN;
            for j = 1:length(c_un)
                r_all(j) = median(abs_lats(r(c == c_un(j))),'omitnan');
            end
             median_lat(m,i) = median(r_all,'omitnan');

        end
    end
end

%plot median latitude of biomes
figure
hold on
plot(1:8,median(median_lat,1,'omitnan'),'+-')
xticklabels({'(1) TRP','(2) HIL','(3) WIS','(4) SUS','(5) HIT ', '(6) MTR',...
    '(7) PEU','(8) SMN'})
grid on

% get sequence of biomes
disp('The sequence in decreasing median latitude is:')
[sorted, sequence] = sort(median(median_lat,1,'omitnan'));
sequence = flip(sequence)


%% Figure 4a Plot dendrogram of clusters/biomes

%dendrogram of biome centroids
corr_new_weights = new_weights(I,:);
corr_new_weights(end,:) = [];

Z = linkage(corr_new_weights,'weighted','cityblock');

yvalues =  {'(1) TRP','(2) HIL','(3) WIS','(4) SUS','(5) HIT ', '(6) MTR',...
    '(7) PEU','(8) SMN'};

[h,nodes, orig] = dendrogram(Z,0,'Orientation','left');

dendro_clusters = gcf;

figure(dendro_clusters);
hold on;
yticklabels(yvalues(orig));
xlabel('Manhattan distance')
set(h,'LineWidth',3)


%% Figure 4b Non-metric multi-dimensional scaling

% =========================================================================
% For the manuscript we used "get_NMDS_annual" as it considers all months
% at the same time, instead of each month separately. Thus, being in line
% with our definition of biomes that persist throughout the months and that
% are not specific for only a few months
% =========================================================================


n_clusters = 8;
cd('/net/kryo/work/ursho/Damiano_Presence_data/presence_absence_tables_ensemble_averages/Group_specific_background_approach/Data/01NeuronsError/')
load('Single_run_11.mat')
cd(folder_main)

%NMDS for each month separately
[ Y_dis1,stress1,disparities1,spread,convex_hull,biome_points ] =...
    get_NMDS(net, corr_corrected_monthly_smooth,[No_nan_phyto_simple(:,2:4),classes],8,1:12);

%NMDS for all months at the same time (name of function might not be the
%best)
[Y_dis,stress,disparities,spread,convex_hull_points,points_biome,location_biome ] =...
    get_NMDS_annual(net,corr_corrected_monthly_smooth,[No_nan_phyto_simple(:,2:4),classes], n_clusters);

%% Table 4 Number of species per biome per month
% =========================================================================
%Get number of species per biome per month
% =========================================================================

Season_corr_classes = No_nan_phyto_simple(:,1:5).*NaN;
cc = 1;
for m = 1:12
    for lat = 180:-1:1
        for lon = 1:360
            tmp = corr_corrected_monthly_smooth(m,lat,lon);
            if(~isnan(tmp))
                IDs = corrected_monthly_ID(m,lat,lon);
                Season_corr_classes(cc,:) = [IDs,lon,lat,m,tmp];
                cc = cc + 1;
            end
        end
    end
end
Season_corr_classes(isnan(Season_corr_classes(:,1)),:) = [];


%plot mean richness
Season_obs = [Season_corr_classes(:,1:4), sum(No_nan_phyto_simple(Season_corr_classes(:,1),5:end-1),2)];
[monthly_richness_map] = prepare2plot(Season_obs(:,2:end));
plotSOM(monthly_richness_map,1,NaN)

[mean_richness_map] = nanmean(monthly_richness_map,1);
mean_richness_map(mean_richness_map==0) = NaN;
plotSOM(mean_richness_map,13,NaN)

num_species = NaN(12,8);
for m = 1:12
    for i = 1:n_clusters
        num_species(m,i) = nanmean(monthly_richness_map(m,corr_corrected_monthly_smooth(m,:,:) == i));
    end
end

disp('Median and IQR species richness:')
nanmedian(num_species,1)
prctile(num_species,75,1)-prctile(num_species,25,1)

%get sum for each latitude and reproduce Damianos curve (Fig 1 in Global pattern of phytoplankton diversity...)
sum_lats = NaN(180,1);
for i = 1:180
    sum_lats(i,1) = nanmean(mean_richness_map(1,i,:),3);
end

figure
plot(sum_lats,-90:89)

%Seems to be working properly!!
plotSOM(corrected_monthly_ID,2,NaN)

Season_obs = [Season_corr_classes(:,1:4), No_nan_phyto_simple(Season_corr_classes(:,1),5:end)];

num_species_month = NaN(12,8);
for m = 1:12
    for i = 1:n_clusters
        %get IDs from map
        IDs = corrected_monthly_ID(m,corr_corrected_monthly_smooth(m,:,:) == i);
        num_species_month(m,i) = sum(sum(No_nan_phyto_simple(IDs,5:end-1),1)> 0)
    end
end

% Sanity check
num_species_monthV2 = NaN(12,8);
for m = 1:12
    for i = 1:n_clusters
        %get IDs from map
        tmp = Season_obs(Season_corr_classes(:,end) == i & Season_corr_classes(:,4) == m,5:end-1);
        num_species_monthV2(m,i) = median(sum(tmp,2),'omitnan')
%         num_species_monthV2(m,i) = sum(sum(tmp,1)>0)
    end
end

% Sanity check 2
Season_obs = [Season_corr_classes(:,1:4), sum(No_nan_phyto_simple(Season_corr_classes(:,1),5:end-1),2)];
[mean_richness_map] = prepare2plot(Season_obs(:,2:end));
num_species_monthV3 = NaN(12,8);
for m = 1:12
    for i = 1:n_clusters
        %get IDs from map
        tmp = mean_richness_map(m,corr_corrected_monthly_smooth(m,:,:) == i);
        num_species_monthV3(m,i) = sum(tmp>0)
    end
end


% =========================================================================
% All 3 versions appear to be the same
% =========================================================================

yvalues =  {'(1) TRP','(2) HIL','(3) WIS','(4) SUS','(5) HIT ', '(6) MTR',...
    '(7) PEU','(8) SMN'};
%get rid of zeros
num_species_month(num_species_month == 0) = NaN;

figure
hold on
for i = 1:12
    hh(i) = plot(1:8,num_species_month(i,:),'+');
end

pp = errorbar(1:8,median(num_species_month,1,'omitnan'),...
    median(num_species_month,1,'omitnan')-prctile(num_species_month,25),...
    prctile(num_species_month,75)-median(num_species_month,1,'omitnan'),'r');
  xlim([0.5 8.5])
xticklabels(yvalues)
ylabel('Number of species')
xlabel('Biome')
legend([hh(1), pp],{'Monthly projection','Median \pm IQR'})
    iqr(num_species_month,1)

[~, I ] = sort(median(num_species_month,1,'omitnan'),'descend');

%print the names of biomes in descending number of species
disp('Biomes in descending species richness:')
yvalues(I)


%% Turnover from monthly
turnover_map = NaN(1,180,360);
for i = 1:180
    for j = 1:360
        tmp_all = NaN;
        for m = 1:12
            tmp_all = [tmp_all;corrected_monthly_ID(m,i,j)];
        end
        tmp_all(isnan(tmp_all)) = [];
        if(~isnan(tmp_all))
            %calculate turnover
            dat = No_nan_phyto_simple(tmp_all,5:end-1);
            for m = 2:size(dat,1)
                S1 = sum(dat(m-1,:));
                S2 = sum(dat(m,:));
                [~, c] = find(dat(m-1,:) == 1 & dat(m,:) == 1);
                c = length(c);
                turnover_map(1,i,j) = (S1-c)+(S2-c);
            end
        end
    end
end

plotSOM(turnover_map,1,NaN)
turnover_mean_annual = NaN(n_clusters,2);
for n = 1:n_clusters
    turnover_mean_annual(n,1) = median(turnover_map(corr_smooth_annual_map == n),'omitnan');
    turnover_mean_annual(n,2) = iqr(turnover_map(corr_smooth_annual_map == n));
end
%% Turnover on a seasonal scale
turnover_map = NaN(4,180,360);

seasons = [3 4 5;6 7 8;9 10 11;12,1,2];
for s = 1:4
    for i = 1:180
        for j = 1:360
            tmp_all = NaN;
            for mm = 1:3
                m = seasons(s,mm);
                tmp_all = [tmp_all;corrected_monthly_ID(m,i,j)];
            end
            tmp_all(isnan(tmp_all)) = [];
            if(~isnan(tmp_all))
                %calculate turnover
                dat = No_nan_phyto_simple(tmp_all,5:end-1);
                for m = 2:size(dat,1)
                    S1 = sum(dat(m-1,:));
                    S2 = sum(dat(m,:));
                    [~, c] = find(dat(m-1,:) == 1 & dat(m,:) == 1);
                    c = length(c);
                    turnover_map(s,i,j) = (S1-c)+(S2-c);
                end
            end
        end
    end

end
turnover_mean = NaN(n_clusters,5);

for n = 1:n_clusters
    turnover_mean(n,1) = prctile(turnover_map(corr_season_smooth == n),25);
    turnover_mean(n,2) = median(turnover_map(corr_season_smooth == n),'omitnan');
    turnover_mean(n,3) = prctile(turnover_map(corr_season_smooth == n),75);
    turnover_mean(n,4) = iqr(turnover_map(corr_season_smooth == n));
    turnover_mean(n,5) = mean(turnover_map(corr_season_smooth == n),'omitnan');
end
disp('Mean species turnover:')
turnover_mean


%% Fraction of species shared between biomes

species_biomes = NaN(12,8,536);
for m = 1:12
    for n = 1:n_clusters
        tmp = corrected_monthly_ID(m,corr_corrected_monthly_smooth(m,:,:) == n);
        min_tmp = 80*length(tmp)/100;
        tmp2 = sum(No_nan_phyto_simple(tmp,5:end-1),1,'omitnan');
        tmp2(tmp2 < min_tmp) = 0;
        tmp2(tmp2 > 0) = 1;
        species_biomes(m,n,:) = tmp2;
    end
end

%for each month get the number of species that is equal between two biomes
matrix_shared_species = NaN(12,n_clusters,n_clusters);
for m = 1:12
    for i = 1:n_clusters
        for j = 1:n_clusters
            [r, ~] = find(species_biomes(m,i,:) == 1 & species_biomes(m,j,:) == 1);
            matrix_shared_species(m,i,j) = length(r)/sum(species_biomes(m,i,:));
        end
    end
end

disp('Median fraction of species shared between biomes:')
mean_mat_shared = squeeze(median(matrix_shared_species,1,'omitnan'))

%% Species richness separated by the top largest groups

% =========================================================================
% Calculate the median number of species for total, dino, bacillario, and
% hapto
% =========================================================================

%find rows for dino, bacillario and haptogit

[r_din,c_din] = find(name_genus_phylum(:,3) == 'dinoflagellata');
[r_bac,c_bac] = find(name_genus_phylum(:,3) == 'bacillariophyceae');
[r_hap,c_hap] = find(name_genus_phylum(:,3) == 'haptophyta');


num_species_month_din = NaN(12,8);
for m = 1:12
    for i = 1:n_clusters
        %get IDs from map
        tmp = Season_obs(Season_corr_classes(:,end) == i & Season_corr_classes(:,4) == m,r_din+4);
        num_species_month_din(m,i) = sum(sum(tmp,1)>0);
    end
end
num_species_month_din(num_species_month_din == 0) = NaN;

num_species_month_bac = NaN(12,8);
for m = 1:12
    for i = 1:n_clusters
        %get IDs from map
        tmp = Season_obs(Season_corr_classes(:,end) == i & Season_corr_classes(:,4) == m,r_bac+4);
        num_species_month_bac(m,i) = sum(sum(tmp,1)>0);
    end
end
num_species_month_bac(num_species_month_bac == 0) = NaN;

num_species_month_hap = NaN(12,8);
for m = 1:12
    for i = 1:n_clusters
        %get IDs from map
        tmp = Season_obs(Season_corr_classes(:,end) == i & Season_corr_classes(:,4) == m,r_hap+4);
        num_species_month_hap(m,i) = sum(sum(tmp,1)>0);
    end
end
num_species_month_hap(num_species_month_hap == 0) = NaN;

disp('Median and IQR of total and group-specific species richness:')
total_species = [median(num_species_month,1,'omitnan');iqr(num_species_month,1);...
                median(num_species_month_din,1,'omitnan');iqr(num_species_month_din,1);...
                median(num_species_month_bac,1,'omitnan');iqr(num_species_month_bac,1);...
                median(num_species_month_hap,1,'omitnan');iqr(num_species_month_hap,1)]

disp('Mean percentage of Dinoflagellata and Bacillariophyceae:')
mean(num_species_month_din./num_species_month,1,'omitnan')
mean(num_species_month_bac./num_species_month,1,'omitnan')

%%
% =========================================================================
% Plot median and IQR to see overlap better for each category m. m = 1 is
% the total, m = 3 is for Dinoflagellata, m = 5 for Bacillariophyceae, m =
% 7 for Haptophyta. See definition of total_species in l.583-586
% =========================================================================
category = [1,3,5,7];

figure
hold on
for m_idx = 1:length(category)
    m = category(m_idx);
    switch m
        case 1
            str = 'Total species';
        case 3
            str = 'Dinoflagellata';
        case 5
            str = 'Bacillariophyceae';
        case 7
            str = 'Haptophyta';
    end

    for i = 1:n_clusters
        plot([total_species(m,i)-total_species(m+1,i),...
            total_species(m,i),total_species(m,i)+total_species(m+1,i)],[i i i])
        plot(total_species(m,i),i,'k*')
    end

    grid on
    ylim([0.5 8.5])
    title(str)
    set (gca,'YDir','reverse')
end

%% Figure 5 Core-Satellite hypothesis
% =========================================================================
% Calculate the area coverage of each species and determine whether they
% are a core or a satellite species
% =========================================================================

tmp_area_map = area_map;
[coverage_month, coverage] = get_occupancy( tmp_area_map, Season_obs,8,corr_corrected_monthly_smooth);

%get average coverage of all species that are present in at least 1 month
for i = 1:n_clusters
    mean(coverage(i,coverage(i,:)>0),'omitnan')
end
% get top 10 species for each biome
top_10 = NaN(n_clusters,10);
for n = 1:n_clusters
    [~, I] = sort(coverage(n,:),'descend');
    top_10(n,:) = I(1:10);
end

disp('Top 10 species per biome:')
coverage(n,top_10(n,:))

% =========================================================================
% Calculate sequence of species based on monthly global area coverage
% =========================================================================

% get heatmaps
yvalues =  {'(1) TRP','(2) HIL','(3) WIS','(4) SUS','(5) HIT ', '(6) MTR',...
    '(7) PEU','(8) SMN'};
[top_occ,I] = get_heatmaps(coverage, yvalues, n_clusters,name_genus_phylum);

%save coverage for further analysis
cd('/net/kryo/work/ursho/PhD/Projects/Biomes/Scripts/Biomes_PROOCE/Data/07Analysis/IndicatorSpecies')
if isfile('Coverage.mat')
    disp('File already exists!')
else
    save('Coverage','coverage_month','coverage','yvalues','name_genus_phylum','top_occ','I')
end

cd(folder_main)
%%
% =========================================================================
% General description
% =========================================================================
%correct wrong entries
name_genus_phylum(237,3) = "dictyochophyceae";
name_genus_phylum(239,3) = "dictyochophyceae";

%% Table A.15 (not used in final version)

%load data
cd('/net/kryo/work/ursho/PhD/Projects/Biomes/Scripts/Biomes_PROOCE/Data/07Analysis/IndicatorSpecies')
load('Coverage.mat')

% correlation between barcodes of biomes
R_species = corr(top_occ','Type','Spearman') ;
R_species(R_species == 1) = NaN;
mean(R_species,2,'omitnan')

%% Figure A.16
% =========================================================================
% Construct maps of coverage of species 1 179 358 and 536
% =========================================================================

my_species = No_nan_phyto_simple(:,[2 3 4 I(1)+4 I(179)+4 I(358)+4 I(536)+4]);
% my_species = No_nan_phyto_simple(:,[2 3 4 I+4]);

disp('The names of the species are:')
my_names = name_genus_phylum([I(1) I(179) I(358) I(536)],1)

for i = 4:7
    my_map = prepare2plot(my_species(:,[1 2 3 i]));
    my_map = sum(my_map,1,'omitnan');
    my_map(my_map == 0) = NaN;
    plotSOM(my_map,1,NaN)
    cmap = ametrine;
    colormap(cmap);
    cc = colorbar;
    caxis([1 12]);
    cc.YTick = [1:12];
    cc.TickLabels = {'1','2','3','4','5', '6','7','8','9','10','11','12'};
    ylabel(cc, 'Number of months')
     set(findall(gcf,'-property','FontSize'),'FontSize',30)
    set(findall(gcf,'-property','LineWidth'),'LineWidth',3)
%     title(my_names(i-3))
    jj = input('next? (Press Enter to continue)');
end

name_genus_phylum([I(1) I(179) I(358) I(536)],1)
unique(my_map(~isnan(my_map)))

%% Overall coverage of tripos extensus globally

my_map = prepare2plot(No_nan_phyto_simple(:,[2 3 4 I(1)+4]));
for i = 1:12
    tmp = squeeze(my_map(i,:,:).*area_map);
    monthly_area = sum(sum(squeeze(area_map(~isnan(my_map(i,:,:))))));
    tmp_sum(i) = sum(sum(tmp,'omitnan'),'omitnan')/monthly_area;
end
mean(tmp_sum)
plotSOM(my_map,1,NaN)
%% Get IDs of top 100 species coverage for each biome
%{
% %mark species in each biome that have at least one month with more than 50%
% %of area coverage
%
% presence_over_thershold = coverage_month.*0;
% presence_over_thershold(coverage_month > 0.5) = 1;
% presence_over_thershold(isnan(coverage_month)) = NaN;
% presence_over_thershold = squeeze(sum(presence_over_thershold,1,'omitnan'));
%
% %get only species that are present in a biome (with the above threshold)
% %for the absolute majority of months
% threshold_presence = squeeze(sum(~isnan(coverage_month),1,'omitnan'));
% present_species = presence_over_thershold.*0;
% present_species(presence_over_thershold >= threshold_presence) = 1;
%
% % sum_presence_species = sum(present_species,2)
% sum_presence_species = sum(present_species,1);
% sum_presence_species(sum_presence_species > 1) = 1;
% sum(sum_presence_species)
% %get the individual species and coverage of species in present_species
%
% for n = 1:length(sum_presence_species)
%     tmp = coverage(n,present_species(n,:) == 1)
% end
%
% for n = 1:8
% [B,tmp] = sort(coverage(n,:),'descend');
% [~,r]=sort(tmp)
% ranks(:,n) = r'
% end
%}
cd('/net/kryo/work/ursho/PhD/Projects/Biomes/Scripts/Biomes_PROOCE/Data/07Analysis/IndicatorSpecies')
load('Coverage.mat')

% Create table of top100 species for each biome
n_clusters = 8;
n_features = size(coverage,2);

ranks = NaN(n_features,n_clusters);
for n = 1:n_clusters
    [~,tmp] = sort(coverage(n,:),'descend');
    [~,r] = sort(tmp);
    ranks(:,n) = r';
end

[rr, ~] = find(ranks <= 100);

top100_spec = name_genus_phylum(unique(rr),[3 1]);
top100_num = coverage(:,unique(rr))';
top100_ranks = ranks(unique(rr),:);

top100_num(top100_ranks > 100) = NaN;

cd(folder_main)

%% Indicator species analysis

% =========================================================================
% Categorize species into satellite (1), subordinate (2) and satellite (3)
% =========================================================================

%load species coverage data
cd('/net/kryo/work/ursho/PhD/Projects/Biomes/Scripts/Biomes_PROOCE/Data/07Analysis/IndicatorSpecies')
load('Coverage.mat')
% get sorted list of species

sorted_species_list = name_genus_phylum(I,:);

sorted_species_coverage = coverage(:,I)';

[~, c] = find(coverage(2,:) == max(coverage(2,:)));
name_genus_phylum(c,:)
mean(coverage,2)
std(coverage,[],2)
max(coverage,[],2)
min(coverage,[],2)
tmp_coverage = coverage_month;

% Categorize into satellite (1), subordinate (2) and satellite (3)
core_species = coverage_month.*0;
satellite_species = coverage_month.*0;
subordinate_species = coverage_month.*0;
for m = 1:12
    tmp_tmpcoverage = squeeze(tmp_coverage(m,:,:));
    for i = 1:size(tmp_tmpcoverage,1)
        [r1, c1] = find(tmp_tmpcoverage(i,:) > 0.9);
        if(~isempty(c1))
            core_species(m,i,c1) = 1;
        end
        [r1, c2] = find(tmp_tmpcoverage(i,:) <= 0.9 & tmp_tmpcoverage(i,:) > 0.1);
        if(~isempty(c2))
            subordinate_species(m,i,c2) = 1;
        end
        [r1, c3] = find(tmp_tmpcoverage(i,:) <= 0.1);
        if(~isempty(c3))
            satellite_species(m,i,c3) = 1;
        end
    end
end

%get core-sub-sat species across all months

for i = 1:12
    [sum(squeeze(core_species(i,:,:)),2,'omitnan'),...
        sum(squeeze(subordinate_species(i,:,:)),2,'omitnan'),...
        sum(squeeze(satellite_species(i,:,:)),2,'omitnan')]
    jj = input('next? (Press Enter to continue)')
end

%get core-sub-sat species (aggregated) over all months
tmp_core_species = core_species;
tmp_subordinate_species = subordinate_species;
tmp_satellite_species = satellite_species;

tmp_core_species(isnan(tmp_core_species)) = 0;
tmp_subordinate_species(isnan(tmp_subordinate_species)) = 0;
tmp_satellite_species(isnan(tmp_satellite_species)) = 0;

tmp_core_species = squeeze(sum(tmp_core_species,1));
tmp_subordinate_species = squeeze(sum(tmp_subordinate_species,1));
tmp_satellite_species = squeeze(sum(tmp_satellite_species,1));

tmp_core_species(tmp_core_species > 0) = 1;
tmp_subordinate_species(tmp_subordinate_species > 0) = 1;
tmp_satellite_species(tmp_satellite_species > 0) = 1;

[sum(tmp_core_species,2,'omitnan'),...
        sum(tmp_subordinate_species,2,'omitnan'),...
        sum(tmp_satellite_species,2,'omitnan')]


cd(folder_main)

%% Get indicator species from monthly data
satellite_species(isnan(satellite_species)) = 1;
for m = 1:12

%     tmp = corr_corrected_monthly_smooth;
%     tmp(corr_corrected_monthly_smooth==9) = NaN;
%     num_biomes = length(unique(tmp(m,~isnan(tmp(m,:,:)))));
    num_biomes = length(unique(corr_corrected_monthly_smooth(m,~isnan(corr_corrected_monthly_smooth(m,:,:)))));


    sum_core = sum(squeeze(core_species(m,:,:)),1,'omitnan');
    sum_core(sum_core > 1) = 0;
    [r, c] = find(sum_core ~= 0);

    sum_sat = sum(squeeze(satellite_species(m,:,:)),1,'omitnan');
    [r, c] = find(sum_core ==1 & sum_sat >=num_biomes-1);
    ind_species = c

end
% =========================================================================
% There are no indicator species on monthly or on annually averaged biomes
% =========================================================================

% find at which level one would find indicator species
ind_species = NaN(8,2);
for i = 1:8
    tmp_ind = NaN;
    tmp_bio = NaN;
    for m = 1:12

        sum_core = sum(squeeze(core_species(m,:,:)),1,'omitnan');
        sum_core(sum_core > 1) = 0;

        sum_sat = sum(squeeze(satellite_species(m,:,:)),1,'omitnan');

        [r_2, c] = find(sum_core ==1 & sum_sat >=i-1)
        if(~isempty(c))
            tmp_ind = [tmp_ind,c]
            tmp = squeeze(core_species(m,:,:))
            [r, c_1] = find(tmp(:,c) == 1)
            tmp_bio = [tmp_bio;r]
        end

    end
    ind_species(i,1) = length(unique(tmp_bio(~isnan(tmp_bio))))
    ind_species(i,2) = length(unique(tmp_ind(~isnan(tmp_ind))))
end



figure
hold on
grid on
plot(0:7,ind_species(:,2))
xlabel('Number of biomes where an indicator species is a satellite species')
ylabel('Number of indicator species')



set(findall(gcf,'-property','FontSize'),'FontSize',30)
set(findall(gcf,'-property','LineWidth'),'LineWidth',3)




%% find indicator species, i.e. core species in only one biome
sum_core = squeeze(nansum(core_species(1:8,:),1))
sum_core(sum_core > 1) = 0;
[r c] = find(sum_core ~= 0)
% =========================================================================
% There are 75 core species that are specific for a biome
% =========================================================================

%make table with these 75 species, and mark them as core, subordinate or
%satellite

%choose threshold for number of satellite as at least 4 satellite
sum_sat = squeeze(sum(satellite_species(1:8,:),1))
tabulate(sum_sat)

num_Satspecies = ones(1,8)
% Indicator species defined as those that are satellite species in most
% other biomes
for i = 1:8


    [r c] = find(sum_core ~=0 & sum_sat >=i-1)
    num_Satspecies(i) = length(r)
end
figure
hold on
plot(0:7,num_Satspecies)
 %[q b] = findchangepts(num_Satspecies,'MaxNumChanges',1,'Statistic','linear')
% three as optimal --> 49 indicator species if we amend the definition, but
% then they are not rather "fuzzy indicator species"

%% Figure A.13
%Indicator species Table
num_biome_ind = (1:8).*NaN;
for i = 1:8
    sum_sat = squeeze(nansum(satellite_species,1));


    [r, c] = find(sum_core ~=0 & sum_sat >=i-1);
    ind_species = c;
    %find number of biomes with indicator species under different scenarios

    mat_cov =coverage(:,ind_species)';
    [rr, cc] = find(mat_cov >0.9);

    num_biome_ind(i) = length(unique(cc));
    [i-1 length(unique(cc))]

end

figure
hold on
grid on
plot(0:7,num_biome_ind)
xlabel('Number of biomes where an indicator species is a satellite species')
ylabel('Number of biomes with indicator species')
ylim([0 8])
yticks([0:8])

vline = line([2 2], [0 8])
vline.Color = 'k';
vline.LineStyle = '-.'

set(findall(gcf,'-property','FontSize'),'FontSize',30)
set(findall(gcf,'-property','LineWidth'),'LineWidth',3)

% =========================================================================
%at i = 3, we get the most biomes with indicator species and reduce the
%number of satellite species. UHE 14 AUG 2019
% =========================================================================

i = 3
sum_sat = squeeze(nansum(satellite_species,1));
[r c] = find(sum_core ~=0 & sum_sat >=i-1);
ind_species = c;


%get the names
labs_core = name_genus_phylum(ind_species,:)
mat_cov =coverage(:,ind_species)';

% =========================================================================
% ============== THERE ARE NO INDICATOR SPECIES! ==========================
% =========================================================================


%% Network analysis

% =========================================================================
% Here we perform our species network analysis
% =========================================================================

%load data needed
cd('/net/kryo/work/ursho/PhD/Projects/Biomes/Scripts/Biomes_PROOCE/Data/05Biomes/')
load('Seasonally_corrected_original')
cd(folder_main)


tic
Scores = monthly_Dunning(corr_corrected_monthly_smooth,Season_obs,n_clusters);
toc

cd('/net/kryo/work/ursho/PhD/Projects/Biomes/Scripts/Biomes_PROOCE/Data/07Analysis/Networks/')
if isfile('monthlyScores04Oct2019.mat')
    disp('File already exists!')
else
    save('monthlyScores04Oct2019','Scores','Season_obs')
end

cd(folder_main)
%% get monthly species pairs

cd('/net/kryo/work/ursho/PhD/Projects/Biomes/Scripts/Biomes_PROOCE/Data/07Analysis/IndicatorSpecies/')
load('Coverage.mat')
cd('/net/kryo/work/ursho/PhD/Projects/Biomes/Scripts/Biomes_PROOCE/Data/07Analysis/Networks/')
load('monthlyScores04Oct2019.mat')

filter_Scores = Scores;

all_biome_pairs = [NaN NaN NaN NaN];
for m = 1:12
    tmp_Scores = squeeze(filter_Scores(:,m,:,:));
    tmp_coverage = squeeze(coverage_month(m,:,:));

    [biome_pairs] = get_species_pairs(tmp_Scores,tmp_coverage);
    biome_pairs
    all_biome_pairs = [all_biome_pairs;biome_pairs];
    m

end
all_biome_pairs(1,:) = []



[B, I] = sort(all_biome_pairs(:,1),'ascend');
sorted_all_biome_pairs = all_biome_pairs(I,:);

% get the network species
final_pairs = [NaN NaN NaN];
for i = 1:8
    tmp = sorted_all_biome_pairs(sorted_all_biome_pairs(:,1) == i,:);
    if(~isempty(tmp))
        connected_edges =  tmp(:,2:3);
        tmp(:,4)
        [connected_edges,num_pairings] = find_most_connected([tmp(:,2), tmp(:,3)],tmp(:,4));
        num_pairings
        final_pairs = [final_pairs; [connected_edges, connected_edges(:,1).*0+i]];

    end

end

final_pairs(1,:) = []
size(final_pairs)


network_species = unique([final_pairs(:,[1,3]);final_pairs(:,[2,3])],'rows')
size(unique(network_species(:,1)))


%save data
cd('/net/kryo/work/ursho/PhD/Projects/Biomes/Scripts/Biomes_PROOCE/Data/07Analysis/Networks/')
if isfile('Network_species_11Oct2019.mat')
    disp('File already exists!')
else
    save('Network_species_11Oct2019','network_species','final_pairs')
end

cd(folder_main)

%% get monthly network species


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

copy_combination = final_pairs.*NaN;
correspondence_labels = [NaN NaN];
new_ID = 1;
flag = 0;
for i = 1:length(network_species)
    r1 = find(final_pairs(:,1) == network_species(i,1))
    r2 = find(final_pairs(:,2) == network_species(i,1))


    for j = 1:length(r1)
        if(isnan(copy_combination(r1(j),1)))
            copy_combination(r1,1) = new_ID;
            flag = 1;
        end
    end

    for j = 1:length(r2)
        if(isnan(copy_combination(r2(j),2)))
            copy_combination(r2,2) = new_ID;
            flag = 1;
        end
    end
    if(flag == 1)
        correspondence_labels = [correspondence_labels;[network_species(i,1),new_ID]];
     new_ID = new_ID + 1;
     flag = 0
    end

end
correspondence_labels(1,:) = [];
copy_combination(:,3) = final_pairs(:,3);


%for plotting purposes get rid of duplicates
[plotting_comb, ia, ic] = unique(copy_combination(:,1:2),'rows','stable');

classes_plotting = final_pairs(ia,3);
plotting_comb = [plotting_comb,classes_plotting];
plotting_comb(isnan(plotting_comb(:,1)),:) = [];
    final_pairs(ia,:)

G = graph(plotting_comb(:,1)', plotting_comb(:,2)');

Connected_comp = conncomp(G);

colormap(comb_cmap);

figure()
h = plot(G)
hold on;
n_clusters = 8
% highlight the different clusters with different colors
for i = 1:n_clusters

    tmp = plotting_comb(plotting_comb(:,3) == i,1:2);
    tmp = unique(reshape(tmp,1,[]));
    highlight(h,[tmp],'NodeColor',comb_cmap(i,:));
    hold off
    jj = input('next? (Press Enter to continue)')
end
h.NodeLabel = correspondence_labels(:,1);




%% Get tables

%get rank of observations

PotScores = Scores;
for m = 1:12
    coverage_data = squeeze(coverage_month(m,:,:));
    upper_area = prctile(coverage_data,90,2);
    lower_area = prctile(coverage_data,10,2);
    for n = 1:size(PotScores,1)

        for i = 1:size(PotScores,3)

            for j = i+1:size(PotScores,4)
                %get the value for each pairing
                if(min(coverage_data(n,i),coverage_data(n,j)) <= lower_area(n) || Scores(n,m,i,j) <= 0 )
                    PotScores(n,m,i,j) = NaN;
                    PotScores(n,m,j,i) = NaN;
                end
            end
        end
    end
end







PotScores = Scores;
for m = 1:12
    coverage_data = squeeze(coverage_month(m,:,:));
    upper_area = prctile(coverage_data,90,2);
    lower_area = prctile(coverage_data,10,2);
    for n = 1:size(PotScores,1)

        for i = 1:size(PotScores,3)

            for j = i+1:size(PotScores,4)
                %get the value for each pairing
                if(min(coverage_data(n,i),coverage_data(n,j)) <= lower_area(n) || Scores(n,m,i,j) <= 0 )
                    PotScores(n,m,i,j) = NaN;
                    PotScores(n,m,j,i) = NaN;
                end
            end
        end
    end
end


%of species pair, and values in each biome (SCORE)
final_pairs
%get for each pair the area coverage
% Scores(2,18,72)
mat_pair = NaN(length(final_pairs),8);
sorted_all_biome_pairs
months_found = NaN
for i = 1:size(final_pairs,1)

    first = min(final_pairs(i,1:2));
    second = max(final_pairs(i,1:2));
    [r c] = find(sorted_all_biome_pairs(:,1) == final_pairs(i,3) &...
        sorted_all_biome_pairs(:,2) == first &...
        sorted_all_biome_pairs(:,3) == second)
    tmp =    max(sorted_all_biome_pairs(r,4))
    %find tmp in Scores(
    tmp2 = Scores(final_pairs(i,3),:,first,second)
    [r m] = find(tmp2 == tmp)%c is  the month

    tmp3 = Scores(:,m,first,second)

    months_found= [months_found,m]
    mat_pair(i,:) = tmp3(1:end)';
end
months_found(1) = [];
months_found = months_found'
%coverage of species in biome

mat_pair_area_sp1 = NaN(length(final_pairs),8);
mat_pair_area_sp2 = mat_pair_area_sp1;
for i = 1:size(final_pairs,1)

    mat_pair_area_sp1(i,:) = coverage_month(months_found(i),:,final_pairs(i,1))';
    mat_pair_area_sp2(i,:) = coverage_month(months_found(i),:,final_pairs(i,2))';
%     mat_pair_area(i,2) = coverage(:,final_pairs(i,2));
end
name_genus_phylum(237,3) = "dictyochophyceae";
name_genus_phylum(239,3) = "dictyochophyceae";
name_genus_phylum(unique(network_species(:,1)),[1,3])
name_genus_phylum(final_pairs(:,1),[3,1])
name_genus_phylum(final_pairs(:,2),[3,1])

%% Create a matrix with a flag for those that are significant

flag_significant = mat_pair.*0;
for i = 1:size(final_pairs,1)
    flag_significant(i,final_pairs(i,3)) = 1;
end

%% schemaball

% For plotting purposes, mark values that are significant for each biome
% and map to interval 0.5 to 0.99
%all other values are mapped to interval -0.4 to 0

%this is needed, since the algorithm only takes values between -1 and 1,
%and we are mostly interested in showing the differences between biomes in
%their species networks. The value of the interaction is then shown in the
%table

original_mat_pair = mat_pair;
mat_pair = original_mat_pair;
mat_pair(mat_pair==0) = NaN;
mat_pair_tmp = mat_pair.*NaN;

for i = 1:8
    [rp, cp] = find(flag_significant(:,i) == 1);
    if(~isempty(rp))
       %scale to 0.1 to 1
       if(length(rp) == 1)
           tmp = 0.99;
       else
           tmp = mapminmax(mat_pair(rp,i)',0.5,0.99);
       end
       for j = 1:length(rp)
           mat_pair_tmp(rp(j),i) = tmp(j);
       end
    end
    [rn, cn] = find(flag_significant(:,i) == 0);
    if(~isempty(rn))
       %scale to 0.1 to 1
       if(length(rn) == 1)
           tmp = -0.4;
       else
           tmp = mapminmax(mat_pair(rn,i)',0,-0.4);
       end
       for j = 1:length(rn)
           mat_pair_tmp(rn(j),i) = tmp(j);
       end
    end
end

final_pairs(:,1)
final_pairs(:,2)

%transform mat_pair to a NxN matrix

all_species_net = unique([final_pairs(:,1);final_pairs(:,2)])
name_genus_phylum(all_species_net,[3,1])


for i = 1:8
    %for each biome i
    matrix_schemaball = NaN(length(all_species_net));
    for j = 1:size(mat_pair_tmp,1)
        if mat_pair_tmp(j,i) ~= 0
            [r1 c1] = find(all_species_net == final_pairs(j,1));
            [r2 c2] = find(all_species_net == final_pairs(j,2));

            matrix_schemaball(r1,r2) = mat_pair_tmp(j,i);
            matrix_schemaball(r2,r1) = mat_pair_tmp(j,i);
        end

    end
    %plot
    matrix_schemaball(matrix_schemaball==0) = NaN;
    if i == 2
        schemaball(matrix_schemaball,cellstr(num2str(all_species_net)),[0.5 0.5 0.5; 1 0 0],[0, 0, 0],[1 0 0])
    else
        schemaball(matrix_schemaball,cellstr(num2str(all_species_net)),[0.5 0.5 0.5; 1 0 0],[0, 0, 0],[0 0 0])
    end

end

%%

f = figure;
hold on;
p = uipanel('Parent',f,'BorderType','none');
%p.Title = 'Normalized Occurrence of PFTs';
p.TitlePosition = 'centertop';
p.FontSize = 12;
p.FontWeight = 'bold';
p.BackgroundColor = [1 1 1];
j = 1;
for i= Spec_seq%1:size(species_per_PFT,2)
    s = subplot(3,3,j,'Parent',p)
    s.Visible = 'off'
    tmp = [new_yearly_classes(:,1:2), species_per_PFT(:,i)];
    species_map = prepare2plot(tmp,4,LatLon,0,0);
    species_map(boundaries ==1) = NaN;
    h = plotSOM(species_map,13,latchl,lonchl,104);
    fig = gcf;
    figure(fig)
    hold on;
    title(unique_pft(i))
    colorbar
    colormap(jet)
    h1 = get(gca,'children')
    copyobj(h1,s)
    hold off;
    j = j+1;
end



aa = num2cell(all_species_net)

aa{1}


%% Species networks SOM

coded_phyto = No_nan_phyto_simple(:,[1:4,unique(network_species(:,1))'+4,end]);

cd(folder_main)

%% Train SOM with network species

% =========================================================================
% Run the following snippet on a cluster not your local machine!
% =========================================================================

cd('/net/kryo/work/ursho/PhD/Projects/Biomes/Scripts/Biomes_PROOCE/Data/07Analysis/Networks/')
tic
d1 = 31;
d2 = d1;
optimal_epoch = 200;
[netw_classes, netw_net] = My_SOM( coded_phyto, d1,d2, optimal_epoch,'mandist' );
toc

%save data
cd('/net/kryo/work/ursho/PhD/Projects/Biomes/Scripts/Biomes_PROOCE/Data/07Analysis/Networks/')
if isfile('Coded_values_SOM_11Oct2019_V3_species.mat')
    disp('File already exists!')
else
    save('Coded_values_SOM_11Oct2019_V3_species','netw_net','netw_classes','coded_phyto','network_species','-v7.3')
end

cd(folder_main)




%% Construct network biomes

tic
original_weights = netw_net.IW{1};
%original_weights = bsxfun(@minus,original_weights,mean(original_weights));
[coeff,score,latent,tsquared,explained] = pca(original_weights);

[r,c] = find(latent > 1);
clearvars latent
%[coeff,score,latent,tsquared,explained]
[coeff,score,~,~,explained] = pca(original_weights,'NumComponents',r(end));
toc

[raw_monthly_maps_network, new_weights_network,~,~] = Calculate_biomes(netw_net, netw_classes,...
    coded_phyto, 9,coeff);

plotSOM(raw_monthly_maps_network,1,9)

[smooth_map_network] = Clean_up_biomes( raw_monthly_maps_network,new_weights_network,...
area_map,0.5,4,0);

plotSOM(smooth_map_network,1,9)
 [annual_map_smooth_netw, annual_map_netw] = aggregate_months(smooth_map_network,new_weights_network,area_map,coded_phyto,0.5,0,4);

plotSOM(annual_map_smooth_netw,1,9)

% =========================================================================
% Compare the monthly biomes; UHE 16.08.2019
% =========================================================================
n_choice = 9
%     load('CompleSuitePhyto.mat')
%     load('Simple_sort_Data.mat')
cd('/net/kryo/work/ursho/PhD/Projects/Biomes/Scripts/Biomes_PROOCE/Data/05Biomes/')
load('No_mean_PCA_biomes_9_v2_5_perc.mat')
cd(folder_main)

ref_weights = new_weights;
ref_map = smooth_map;

new_map = smooth_map_network;
removed_weights = new_weights_network;

tmp_ref = ref_weights(:,network_species);

% correspondence part
[metric_1,metric_2,metric_3,overlap_pairs1,overlap_pairs2]  = compare_overlap(ref_map,new_map,area_map,tmp_ref,removed_weights);

overlap_maps = NaN(12,180,360);
for i = 1:length(overlap_pairs1)
    overlap_maps(smooth_map == overlap_pairs1(i,1) & smooth_map_network == overlap_pairs2(i,2)) = 1;
end

num_months = sum(~isnan(smooth_map),1);
plotSOM(num_months,1,NaN)
sum_overlap = sum(overlap_maps,1,'omitnan');
plotSOM(sum_overlap./num_months,1,NaN)

for m = 1:12
    plotSOM(smooth_map_network,m,9)
end

%% Environmental analysis

% Use the respective scripts (Niches_from_model, Get_environmental_data,
% and scripts in folder functions/Environmental

%% Get niches from model
%plot distribution of the cooccurrences for each biome
Niches_from_model
