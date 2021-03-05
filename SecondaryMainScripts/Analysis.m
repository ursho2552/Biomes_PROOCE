clear all
folder_main = '/net/kryo/work/ursho/Damiano_Presence_data/presence_absence_tables_ensemble_averages/Group_specific_background_approach'
addpath(genpath(folder_main))
cd(folder_main)

% =========================================================================
% ================== Analysis =============================================
% =========================================================================

load('/net/kryo/work/ursho/Damiano_Presence_data/presence_absence_tables_ensemble_averages/Group_specific_background_approach/Data/01Neurons_error/Single_run_11.mat')
%load help variables
load('/net/kryo/work/ursho/Damiano_Presence_data/presence_absence_tables_ensemble_averages/Group_specific_background_approach/Data/HelpVariables2.mat', 'LatLon')
load('/net/kryo/work/ursho/Damiano_Presence_data/presence_absence_tables_ensemble_averages/Group_specific_background_approach/Data/HelpVariables2.mat', 'latchl')
load('/net/kryo/work/ursho/Damiano_Presence_data/presence_absence_tables_ensemble_averages/Group_specific_background_approach/Data/HelpVariables2.mat', 'lonchl')

%construct or load simplified version of raw data
load('/net/kryo/work/ursho/Damiano_Presence_data/presence_absence_tables_ensemble_averages/Group_specific_background_approach/Data/00Probabilities/Transformed_CompleteSuitePhyto.mat')
load('/net/kryo/work/ursho/Damiano_Presence_data/presence_absence_tables_ensemble_averages/Group_specific_background_approach/Data/00Probabilities/Simple_sort_Data.mat')

load('/net/kryo/work/ursho/Damiano_Presence_data/presence_absence_tables_ensemble_averages/Group_specific_background_approach/Data/06Monthly_biomes/PCA_no_mean_v2/Latent_threshold/No_mean_PCA_biomes_9_v2_5_perc.mat')
load('/net/kryo/work/ursho/Damiano_Presence_data/presence_absence_tables_ensemble_averages/Group_specific_background_approach/Data/06Monthly_biomes/PCA_no_mean_v2/Latent_threshold/No_mean_PCA_biomes_annual_9_v2_5_perc.mat')
load('/net/kryo/work/ursho/Damiano_Presence_data/presence_absence_tables_ensemble_averages/Group_specific_background_approach/Data/06Monthly_biomes/PCA_no_mean_v2/Latent_threshold/No_mean_PCA_biomes_seasonal_9_v2_5_perc.mat')
load('/net/kryo/work/ursho/Damiano_Presence_data/presence_absence_tables_ensemble_averages/Group_specific_background_approach/Data/00Probabilities/Names_species.mat')
load('/net/kryo/work/ursho/Damiano_Presence_data/presence_absence_tables_ensemble_averages/Group_specific_background_approach/Data/06Monthly_biomes/PCA_no_mean_v2/Latent_threshold/Seasonally_corrected_data.mat')

load('/net/kryo/work/ursho/Damiano_Presence_data/presence_absence_tables_ensemble_averages/Group_specific_background_approach/Data/Area_map.mat')



%%
% =========================================================================
% Change label numbering using descending mean area; UHE 09/08/2019
% =========================================================================


n_clusters = 9
area_biome = ones(n_clusters,12).*NaN;
tmp_map = corrected_monthly_smooth;
for m = 1:12
    for i = 1:n_clusters
        tot_area = sum(area_map(~isnan(tmp_map(m,:,:))));
        area_biome(i,m) = sum(area_map(tmp_map(m,:,:) == i))/tot_area;
    end
end
 area_biome(area_biome == 0) = NaN;
[B I] = sort(mean(area_biome,2,'omitnan'),'descend')
aa = area_biome(I,:)'
corr_smooth_annual_map = smooth_annual_map;
corr_season_smooth = ones(4,180,360).*NaN;%season_map_smooth;
corr_corrected_monthly_smooth = ones(12,180,360).*NaN;
corr_corrected_monthly_raw = ones(12,180,360).*NaN;
for i = 1:n_clusters
    
    corr_smooth_annual_map(smooth_annual_map == I(i)) = i;
    corr_season_smooth(season_map_smooth == I(i)) = i;
    corr_corrected_monthly_smooth(corrected_monthly_smooth == I(i)) = i;
    corr_corrected_monthly_raw(corrected_monthly_raw == I(i)) = i;
end

%safe biomes with corrected labels
cd('/net/kryo/work/ursho/Damiano_Presence_data/presence_absence_tables_ensemble_averages/Group_specific_background_approach/Data/05Biomes/')

save('Seasonally_corrected_original','corr_smooth_annual_map','corr_season_smooth','corr_corrected_monthly_smooth')


cd(folder_main)
% get min and max monthly area
%% Data for Table 3
max_area = max(area_biome,[],2);
round(1000*max_area(I))/10
area_biome(isnan(area_biome)) = 0
min_area = min(area_biome,[],2);
round(1000*min_area(I))/10

round(1000*B)/10

%% Figure 2
% =========================================================================
% Plot annual biomes; UHE 09/08/2019
% =========================================================================
plotSOM(corr_smooth_annual_map,1,8)
hold on;
cmap = morgenstemning(12)
[cmap_tmp] = shuffle_colormap(cmap);
[cmap_tmp1] = shuffle_colormap(cmap_tmp);


cmap = ametrine(12)
[cmap_tmp] = shuffle_colormap(cmap);
[cmap_tmp2] = shuffle_colormap(cmap_tmp);

cmap = isolum(12)
[cmap_tmp] = shuffle_colormap(cmap);
[cmap_tmp3] = shuffle_colormap(cmap_tmp);

comb_cmap = [cmap_tmp2([2,1,3],:);cmap_tmp1(7,:);cmap_tmp3(5,:);cmap_tmp1(8,:);cmap_tmp1(9,:);cmap_tmp1(4,:)]

colormap(comb_cmap(1:8,:))
c = colorbar
set( c, 'YDir', 'reverse' );
[positions] = get_ticks_centered(8);
c.YTick = positions;
c.TickLabels = {'(1) TRP','(2) HIL','(3) WIS','(4) SUS','(5) HIT ', '(6) MTR',...
    '(7) PEU', '(8) SMN'}
set(findall(gcf,'-property','FontSize'),'FontSize',30)
set(findall(gcf,'-property','LineWidth'),'LineWidth',3)

%% Figure 3
% =========================================================================
% Plot seasonal biomes; UHE 09/08/2019
% =========================================================================

for i = 1:4
    unique(corr_season_smooth(i,~isnan(corr_season_smooth(i,:,:))))
end


for s = 1:4
    plotSOM(corr_season_smooth,s,8)
    hold on;
    colormap(comb_cmap)
    set(findall(gcf,'-property','FontSize'),'FontSize',30)
    set(findall(gcf,'-property','LineWidth'),'LineWidth',3)
    hold off
end
%% Table A.4
% =========================================================================
% Get area of annual and seasonal biomes; UHE 09/08/2019
% =========================================================================
area_seasonal = ones(9,4).*NaN;
for s = 1:5
    for i = 1:n_clusters
        if(s<5)
            tot_area = sum(area_map(~isnan(corr_season_smooth(s,:,:))));
            area_biomes(i,s) = sum(area_map(corr_season_smooth(s,:,:) == i))/tot_area;
        else
            tot_area = sum(area_map(~isnan(corr_smooth_annual_map(1,:,:))));
            area_biomes(i,s) = sum(area_map(corr_smooth_annual_map(1,:,:) == i))/tot_area;
  
        end
    end
end
round(1000*area_biomes)/10
area_biomes  = area_biomes*100     

%% Figure A.11

diff_maps = ones(4,180,360).*NaN;
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


area_diff;

%% Find average latitude of biomes
n_clusters = 8
abs_latchl = abs(latchl);
median_lat = ones(12,n_clusters).*NaN;
for m = 1:12
    for i = 1:n_clusters
        [r c] = find(squeeze(corr_corrected_monthly_smooth(m,:,:)) == i);
        if(~isempty(r))
            %get for each c the median of r, then the median of r
            c_un = unique(c);
            r_all = c_un.*NaN;
            tmp = NaN;
            for j = 1:length(c_un)
                r_all(j) = median(abs_latchl(r(c == c_un(j))),'omitnan');
            end
             median_lat(m,i) = median(r_all,'omitnan');   
            
            
            
            
        end
    end
end

figure
hold on
plot(1:8,median(median_lat,1,'omitnan'),'+-')
xticklabels({'(1) TRP','(2) HIL','(3) WIS','(4) SUS','(5) HIT ', '(6) MTR',...
    '(7) PEU','(8) SMN'})   
grid on


%%
area_biomes
%% Figure 4a
%dendrogram of weights
corr_new_weights = new_weights(I,:);
corr_new_weights(end,:) = []
Z = linkage(corr_new_weights,'weighted','cityblock');

yvalues =  {'(1) TRP','(2) HIL','(3) WIS','(4) SUS','(5) HIT ', '(6) MTR',...
    '(7) PEU','(8) SMN'}
[h,nodes, orig] = dendrogram(Z,0,'Orientation','left');

dendro_clusters = gcf;

figure(dendro_clusters);
hold on;
yticklabels(yvalues(orig));
xlabel('Manhattan distance')
set(h,'LineWidth',3)


%% Figure 4b 
%NMDS
n_clusters = 8

load('/net/kryo/work/ursho/Damiano_Presence_data/presence_absence_tables_ensemble_averages/Group_specific_background_approach/Data/01Neurons_error/Single_run_11.mat')
%corr_corrected_monthly_raw
[ Y_dis1,stress1,disparities1,spread,convex_hull,biome_points ] =...
    get_NMDS(net, corr_corrected_monthly_smooth,[No_nan_phyto_simple(:,2:4),classes],8,1:12);
stress1
round(10.*spread)/10

maps = prepare2plotV2([No_nan_phyto_simple(:,2:4),classes]);



[Y_dis,stress,disparities,spread,convex_hull_points,points_biome,location_biome ] =...
    get_NMDS_annual(net,corr_corrected_monthly_smooth,[No_nan_phyto_simple(:,2:4),classes], n_clusters);

%% Table 4
% =========================================================================
%Get number of species per biome per month
% =========================================================================


%Calculate Season_obs from corr_corrected_monthly_smooth
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
    m
end
Season_corr_classes(isnan(Season_corr_classes(:,1)),:) = [];


%plot mean richness
Season_obs = [Season_corr_classes(:,1:4), sum(No_nan_phyto_simple(Season_corr_classes(:,1),5:end-1),2)];
[monthly_richness_map] = prepare2plotV2(Season_obs(:,2:end));
plotSOM(monthly_richness_map,1,NaN)

[mean_richness_map] = nanmean(monthly_richness_map,1);
mean_richness_map(mean_richness_map==0) = NaN;
plotSOM(mean_richness_map,13,NaN)

num_species = ones(12,8).*NaN;
for m = 1:12
    for i = 1:n_clusters
        num_species(m,i) = nanmean(monthly_richness_map(m,corr_corrected_monthly_smooth(m,:,:) == i));
    end
end
    
nanmedian(num_species,1)
    
prctile(num_species,75,1)-prctile(num_species,25,1)


%get sum for each latitude and reproduce Damianos curve (Fig 1 in Global pattern of phytoplankton diversity...) 
sum_lats = ones(180,1).*NaN
for i = 1:180
    sum_lats(i,1) = nanmean(mean_richness_map(1,i,:),3);
end

figure
plot(sum_lats,-90:89)

%Seems to be working properly!!
plotSOM(corrected_monthly_ID,2,NaN)


Season_obs = [Season_corr_classes(:,1:4), No_nan_phyto_simple(Season_corr_classes(:,1),5:end)];

num_species_month = ones(12,8).*NaN;
for m = 1:12
    for i = 1:n_clusters
        %get IDs from map
        IDs = corrected_monthly_ID(m,corr_corrected_monthly_smooth(m,:,:) == i);
        num_species_month(m,i) = sum(sum(No_nan_phyto_simple(IDs,5:end-1),1)> 0)
    end
end

% Safety check
num_species_monthV2 = ones(12,8).*NaN;
for m = 1:12
    for i = 1:n_clusters
        %get IDs from map
        tmp = Season_obs(Season_corr_classes(:,end) == i & Season_corr_classes(:,4) == m,5:end-1);
        num_species_monthV2(m,i) = median(sum(tmp,2),'omitnan')
%         num_species_monthV2(m,i) = sum(sum(tmp,1)>0)
    end
end

% Safety check 2
Season_obs = [Season_corr_classes(:,1:4), sum(No_nan_phyto_simple(Season_corr_classes(:,1),5:end-1),2)];
[mean_richness_map] = prepare2plotV2(Season_obs(:,2:end));
num_species_monthV2 = ones(12,8).*NaN;
for m = 1:12
    for i = 1:n_clusters
        %get IDs from map
        tmp = mean_richness_map(m,corr_corrected_monthly_smooth(m,:,:) == i)
        num_species_monthV2(m,i) = sum(tmp>0)
    end
end
%% get turnover
turnover_map = ones(1,180,360).*NaN;
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
                [r c] = find(dat(m-1,:) == 1 & dat(m,:) == 1);
                c = length(c);
                turnover_map(1,i,j) = (S1-c)+(S2-c);
            end
        end
    end
end
%%
plotSOM(turnover_map,1,NaN)
turnover_mean_annual = ones(n_clusters,2);
for n = 1:n_clusters
    turnover_mean_annual(n,1) = median(turnover_map(corr_smooth_annual_map == n),'omitnan');
    turnover_mean_annual(n,2) = iqr(turnover_map(corr_smooth_annual_map == n));
end
%% Turnover on a seasonal scale
turnover_map = ones(4,180,360).*NaN;

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
                    [r c] = find(dat(m-1,:) == 1 & dat(m,:) == 1);
                    c = length(c);
                    turnover_map(s,i,j) = (S1-c)+(S2-c);
                end
            end
        end
    end
    s
end
turnover_mean = ones(n_clusters,5).*NaN;

for n = 1:n_clusters
    turnover_mean(n,1) = prctile(turnover_map(corr_season_smooth == n),25);
    turnover_mean(n,2) = median(turnover_map(corr_season_smooth == n),'omitnan');
    turnover_mean(n,3) = prctile(turnover_map(corr_season_smooth == n),75);
    turnover_mean(n,4) = iqr(turnover_map(corr_season_smooth == n));
    turnover_mean(n,5) = mean(turnover_map(corr_season_smooth == n),'omitnan');
end
turnover_mean


%% get fraction of species shared between biomes
species_biomes = ones(12,8,536).*NaN;
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
matrix_shared_species = ones(12,n_clusters,n_clusters).*NaN;
for m = 1:12
    for i = 1:n_clusters
        for j = 1:n_clusters
            [r c] = find(species_biomes(m,i,:) == 1 & species_biomes(m,j,:) == 1);
            matrix_shared_species(m,i,j) = length(r)/sum(species_biomes(m,i,:));
        end
    end
end
mean_mat_shared = squeeze(median(matrix_shared_species,1,'omitnan'))        

%%


% num_species_monthV3 = ones(12,8).*NaN;
% for m = 1:12
%     for i = 1:n_clusters
% 
%         if(length(r) >= 0)
%             [rr cc] = find(coverage_month(m,i,:) > 0 & ~isnan(coverage_month(m,i,:)));
%             num_species_monthV3(m,i) = length(rr)
% 
%         end
%     end
% end

% =========================================================================
% All versions of num_species_month are equal
% =========================================================================
% get names of biomes
yvalues =  {'(1) TRP','(2) HIL','(3) WIS','(4) SUS','(5) HIT ', '(6) MTR',...
    '(7) PEU','(8) SMN'}
%get rid of zeros
num_species_month(num_species_month == 0) = NaN;

figure
hold on
for i = 1:12
    hh(i) = plot(1:8,num_species_month(i,:),'+')
end

pp = errorbar(1:8,median(num_species_month,1,'omitnan'),...
    median(num_species_month,1,'omitnan')-prctile(num_species_month,25),...
    prctile(num_species_month,75)-median(num_species_month,1,'omitnan'),'r')
  xlim([0.5 8.5])  
xticklabels(yvalues)
ylabel('Number of species') 
xlabel('Biome')
legend([hh(1), pp],{'Monthly projection','Median \pm IQR'})
    iqr(num_species_month,1)
    
[B I ] = sort(median(num_species_month,1,'omitnan'),'descend')
yvalues(I)



% =========================================================================
% Calculate the median number of species for total, dino, bacillario, and
% hapto
% =========================================================================

%find rows for dino, bacillario and haptogit

[r_din,c_din] = find(name_genus_phylum(:,3) == 'dinoflagellata')
[r_bac,c_bac] = find(name_genus_phylum(:,3) == 'bacillariophyceae')
[r_hap,c_hap] = find(name_genus_phylum(:,3) == 'haptophyta')


num_species_month_din = ones(12,8).*NaN;
for m = 1:12
    for i = 1:n_clusters
        %get IDs from map
        tmp = Season_obs(Season_corr_classes(:,end) == i & Season_corr_classes(:,4) == m,r_din+4);
        num_species_month_din(m,i) = sum(sum(tmp,1)>0)
    end
end
num_species_month_din(num_species_month_din == 0) = NaN;

num_species_month_bac = ones(12,8).*NaN;
for m = 1:12
    for i = 1:n_clusters
        %get IDs from map
        tmp = Season_obs(Season_corr_classes(:,end) == i & Season_corr_classes(:,4) == m,r_bac+4);
        num_species_month_bac(m,i) = sum(sum(tmp,1)>0)
    end
end
num_species_month_bac(num_species_month_bac == 0) = NaN;

num_species_month_hap = ones(12,8).*NaN;
for m = 1:12
    for i = 1:n_clusters
        %get IDs from map
        tmp = Season_obs(Season_corr_classes(:,end) == i & Season_corr_classes(:,4) == m,r_hap+4);
        num_species_month_hap(m,i) = sum(sum(tmp,1)>0)
    end
end
num_species_month_hap(num_species_month_hap == 0) = NaN;

total_species = [median(num_species_month,1,'omitnan');iqr(num_species_month,1);...
                median(num_species_month_din,1,'omitnan');iqr(num_species_month_din,1);...
                median(num_species_month_bac,1,'omitnan');iqr(num_species_month_bac,1);...
                median(num_species_month_hap,1,'omitnan');iqr(num_species_month_hap,1)]


mean(num_species_month_din./num_species_month,1,'omitnan')
mean(num_species_month_bac./num_species_month,1,'omitnan')

%%
% =========================================================================
% Plot median and IQR to see overlap better
% =========================================================================
m = 7

figure
hold on
for i = 1:n_clusters
    plot([total_species(m,i)-total_species(m+1,i),...
        total_species(m,i),total_species(m,i)+total_species(m+1,i)],[i i i])
    plot(total_species(m,i),i,'k*')
end
grid on        
ylim([0.5 8.5])
set (gca,'YDir','reverse')
%% Figure 5
% =========================================================================
% Core-Satellite hypothesis
% =========================================================================
tmp_area_map = area_map;%.*0 + 1;
[coverage_month, coverage] = get_occupancy( tmp_area_map, Season_obs,8,corr_corrected_monthly_smooth);

%get average coverage of all species that are present in at least 1 month
for i = 1:n_clusters
    mean(coverage(i,coverage(i,:)>0),'omitnan')
end
% get top 10 species for each biome
top_10 = ones(n_clusters,10).*NaN;
for n = 1:n_clusters
    [B I] = sort(coverage(n,:),'descend');
    top_10(n,:) = I(1:10)
end

aa = coverage(n,top_10(n,:))

% =========================================================================
% Calculate sequence of species based on monthly global area coverage
% =========================================================================


% get heatmaps
yvalues =  {'(1) TRP','(2) HIL','(3) WIS','(4) SUS','(5) HIT ', '(6) MTR',...
    '(7) PEU','(8) SMN'}
[top_occ,I] = get_heatmaps(coverage, yvalues, n_clusters,name_genus_phylum);%,mean_global_area_coverage);

cd('/net/kryo/work/ursho/Damiano_Presence_data/presence_absence_tables_ensemble_averages/Group_specific_background_approach/Data/11Indicator_species')

save('CoverageV2','coverage_month','coverage','yvalues','name_genus_phylum','top_occ','I')
load('CoverageV2.mat')

cd(folder_main)
%%
% =========================================================================
% General description
% =========================================================================
%correct wrong entries
name_genus_phylum(237,3) = "dictyochophyceae";
name_genus_phylum(239,3) = "dictyochophyceae";

%% Table A.15
% correlation between barcodes of biomes
R_species = corr(top_occ','Type','Spearman') 
R_species(R_species == 1) = NaN;
mean(R_species,2,'omitnan')

%% Figure A.12
% =========================================================================
% Construct maps of coverage of species 1 179 358 and 536
% =========================================================================


my_species = No_nan_phyto_simple(:,[2 3 4 I(1)+4 I(179)+4 I(358)+4 I(536)+4]);
my_species = No_nan_phyto_simple(:,[2 3 4 I+4]);

my_names =name_genus_phylum([I(1) I(179) I(358) I(536)],1)
for i = 4:539
    my_map = prepare2plotV2(my_species(:,[1 2 3 i]));
    my_map = sum(my_map,1,'omitnan')
    my_map(my_map == 0) = NaN;
    plotSOM(my_map,1,NaN)
    cmap = ametrine;
    colormap(cmap);
    cc = colorbar;
    caxis([1 12]);
    cc.YTick = [1:12];
    cc.TickLabels = {'1','2','3','4','5', '6','7','8','9','10','11','12'}
    ylabel(cc, 'Number of months')
     set(findall(gcf,'-property','FontSize'),'FontSize',30)
    set(findall(gcf,'-property','LineWidth'),'LineWidth',3)
%     title(my_names(i-3))
    jj = input('dfd')
end
name_genus_phylum([I(1) I(179) I(358) I(536)],1)
unique(my_map(~isnan(my_map)))

%% Overall coverage of tripos extensus globally

my_map = prepare2plotV2(No_nan_phyto_simple(:,[2 3 4 I(1)+4]));
for i = 1:12
    tmp = squeeze(my_map(i,:,:).*area_map);
    monthly_area = sum(sum(squeeze(area_map(~isnan(my_map(i,:,:))))));
    tmp_sum(i) = sum(sum(tmp,'omitnan'),'omitnan')/monthly_area
end
mean(tmp_sum)
plotSOM(my_map,1,NaN)
%%

%mark species in each biome that have at least one month with more than 50%
%of area coverage

presence_over_thershold = coverage_month.*0;
presence_over_thershold(coverage_month > 0.5) = 1;
presence_over_thershold(isnan(coverage_month)) = NaN;
presence_over_thershold = squeeze(sum(presence_over_thershold,1,'omitnan'));

%get only species that are present in a biome (with the above threshold)
%for the absolute majority of months
threshold_presence = squeeze(sum(~isnan(coverage_month),1,'omitnan'));
present_species = presence_over_thershold.*0;
present_species(presence_over_thershold >= threshold_presence) = 1;
aa = sum(present_species,2)
aa = sum(present_species,1)
aa(aa > 1) = 1
sum(aa)
%get the individual species and coverage of species in present_species

for n = 1:length(aa)
    tmp = coverage(n,present_species(n,:) == 1)
end
    
for n = 1:8
[B,tmp] = sort(coverage(n,:),'descend');
[~,r]=sort(tmp)
ranks(:,n) = r'
end

ranks = ones(536,n_clusters).*NaN;
for n = 1:n_clusters
    [B,tmp] = sort(coverage(n,:),'descend');
    [~,r]=sort(tmp)
    ranks(:,n) = r'
end

[rr c] = find(ranks <= 100)

top10_spec = name_genus_phylum(unique(rr),[3 1])
top10_num = coverage(:,unique(rr))'
top10_ranks = ranks(unique(rr),:)


top10_num(top10_ranks > 100) = NaN;

tmp = coverage';

sorted_cov = coverage(n,tmp)';
sorted_cov_name = name_genus_phylum(tmp,[3 1]);

%% Indicator species analysis

% get sorted list of species

sorted_species_list = name_genus_phylum(I,:);

sorted_species_coverage = coverage(:,I)';

[r c] = find(coverage(2,:) == max(coverage(2,:)))
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
        [r1 c1] = find(tmp_tmpcoverage(i,:) > 0.9);
        if(~isempty(c1))
            core_species(m,i,c1) = 1;
        end
        [r1 c2] = find(tmp_tmpcoverage(i,:) <= 0.9 & tmp_tmpcoverage(i,:) > 0.1);
        if(~isempty(c2))
            subordinate_species(m,i,c2) = 1;
        end
        [r1 c3] = find(tmp_tmpcoverage(i,:) <= 0.1);
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
    jj = input('dfd')
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

%% Additional analysis on carbon and  biovolume

% Add carbon and biovolume to name_genus_phylum
% read data from Sal et al (2013) Table 2 manually!!!
save('Sal_et_al_biovolume','Biovolume_names','Biovolume_values')


Values_per_species = ones(length(name_genus_phylum),2).*NaN;
for i = 1:length(name_genus_phylum)
    [r c] = find(Biovolume_names == name_genus_phylum{i,1})
    if(~isempty(r))
        tmp = mean(Biovolume_values(r,:),1,'omitnan');
        Values_per_species(i,:) = tmp;
    end
end
        
%get list of number of species per genera 

[r_bac c] = find(name_genus_phylum(:,3) == 'bacillariophyceae')
bacillario = name_genus_phylum(r_bac,:); 
bac_gen = unique(bacillario(:,2))
 
for i = 1:length(bac_gen)
    [r c] = find(bacillario(:,2) == bac_gen{i})
    bac_gen_num(i,1) = length(r)
end

[r_din c] = find(name_genus_phylum(:,3) == 'dinoflagellata')
dino = name_genus_phylum(r_din,:); 
dino_gen = unique(dino(:,2))
 dino_gen_num = ones(length(dino_gen),1).*NaN;
for i = 1:length(dino_gen)
    [r c] = find(dino(:,2) == dino_gen{i})
    dino_gen_num(i,1) = length(r)
end

[r_hap c] = find(name_genus_phylum(:,3) == 'haptophyta')
hapto = name_genus_phylum(r_hap,:); 
hapto_gen = unique(hapto(:,2))
  hapto_gen_num = ones(length(hapto_gen),1).*NaN;

for i = 1:length(hapto_gen)
    [r c] = find(hapto(:,2) == hapto_gen{i})
    hapto_gen_num(i,1) = length(r)
end

%get carbon and biovolume for each biome
carbon_biovolume_all = ones(12,8,2).*NaN;
for m = 1:12
    for i = 1:n_clusters
        carbon_biovolume_all(m,i,1) = sum(num_species_monthV3(m,:,i)'.*Values_per_species(:,1),1,'omitnan');
        carbon_biovolume_all(m,i,2) = sum(num_species_monthV3(m,:,i)'.*Values_per_species(:,2),1,'omitnan');
    end
end
carbon_biovolume_all(carbon_biovolume_all == 0) = NaN;

carbon_biovolume_bac = ones(12,8,2).*NaN;
for m = 1:12
    for i = 1:n_clusters
        carbon_biovolume_bac(m,i,1) = sum(num_species_monthV3(m,r_bac,i)'.*Values_per_species(r_bac,1),1,'omitnan');
        carbon_biovolume_bac(m,i,2) = sum(num_species_monthV3(m,r_bac,i)'.*Values_per_species(r_bac,2),1,'omitnan');
    end
end
carbon_biovolume_bac(carbon_biovolume_bac == 0) = NaN;

carbon_biovolume_din = ones(12,8,2).*NaN;
for m = 1:12
    for i = 1:n_clusters
        carbon_biovolume_din(m,i,1) = sum(num_species_monthV3(m,r_din,i)'.*Values_per_species(r_din,1),1,'omitnan');
        carbon_biovolume_din(m,i,2) = sum(num_species_monthV3(m,r_din,i)'.*Values_per_species(r_din,2),1,'omitnan');
    end
end
carbon_biovolume_din(carbon_biovolume_din == 0) = NaN;

mean(carbon_biovolume_all,1,'omitnan')
mean(carbon_biovolume_bac,1,'omitnan')
mean(carbon_biovolume_din,1,'omitnan')


%% get all core species
i = 8
[r c] = find(core_species == 1);

all_core_species = unique(c)
labs_core = name_genus_phylum(all_core_species,:)
mat_cov = round(100.*coverage(:,all_core_species)')./100;
Values_per_species(all_core_species,:)

%% Get indicator species from monthly data

for m = 1:12
    
    num_biomes = length(unique(corr_corrected_monthly_smooth(m,~isnan(corr_corrected_monthly_smooth(m,:,:)))))
    sum_core = sum(squeeze(core_species(m,:,:)),1,'omitnan');
    sum_core(sum_core > 1) = 0;
    [r c] = find(sum_core ~= 0);
    
    sum_sat = sum(squeeze(satellite_species(m,:,:)),1,'omitnan');
    [r c] = find(sum_core ==1 & sum_sat >=num_biomes-1);
    ind_species = c
    
end
% =========================================================================
% There are no indicator species on monthly or on annually averaged biomes
% =========================================================================

% find at which level one would find indicator species
ind_species = ones(8,2).*NaN;
for i = 1:8
    tmp_ind = NaN;
    tmp_bio = NaN;
    for m = 1:12
%
        sum_core = sum(squeeze(core_species(m,:,:)),1,'omitnan');
        sum_core(sum_core > 1) = 0;
%         [r_1 c_core] = find(sum_core ~= 0);

        sum_sat = sum(squeeze(satellite_species(m,:,:)),1,'omitnan');

        [r_2 c] = find(sum_core ==1 & sum_sat >=i-1)
        if(~isempty(c))
            tmp_ind = [tmp_ind,c]
            tmp = squeeze(core_species(m,:,:))
            [r c_1] = find(tmp(:,c) == 1)
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
% three as optimal --> 49 indicator species

%% Figure A.13
%Indicator species Table
num_biome_ind = (1:8).*NaN;
for i = 1:8
    sum_sat = squeeze(nansum(satellite_species,1));
    

    [r c] = find(sum_core ~=0 & sum_sat >=i-1);
    ind_species = c;
    %find number of biomes with indicator species under different scenarios
    
    mat_cov =coverage(:,ind_species)';
    [rr cc] = find(mat_cov >0.9);
    
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
Values_per_species(ind_species,:)
%% SOM indicator (NOTE: at level 0 there are no indicator species!!!)

ind_phyto = No_nan_phyto_simple(:,[1:4,(ind_species+4), end]);
%get SOM only using indicator species
tic
d1 = 31;
d2 = d1;
optimal_epoch = 200;
[ind_classes, ind_net] = My_SOM( ind_phyto, d1,d2, optimal_epoch,'mandist' );
toc
cd('/net/kryo/work/ursho/Damiano_Presence_data/presence_absence_tables_ensemble_averages/Group_specific_background_approach/Data/11Indicator_species')
save('Indicator_SOM_48','ind_phyto','ind_classes','ind_net' ,'ind_species')

cd(folder_main)


%% Construct biomes
tic
original_weights = ind_net.IW{1};
%original_weights = bsxfun(@minus,original_weights,mean(original_weights));
[coeff,score,latent,tsquared,explained] = pca(original_weights);

[r,c] = find(latent > 1)
if(isempty(r))
    clearvars latent
    [r,c] = find(explained > 10)
end
%[coeff,score,latent,tsquared,explained]
[coeff,score,~,~,explained] = pca(original_weights,'NumComponents',r(end));
toc

[raw_monthly_maps_indicator, new_weights_indicator,~,~] = Calculate_biomesV2(ind_net, ind_classes,...
    No_nan_phyto_simple, 9,coeff);

[smooth_map_indicator] = Clean_up_biomesV3( raw_monthly_maps_indicator,new_weights_indicator,...
area_map,0.5,4,0);

[annual_map_smooth_ind, annual_map_ind] = aggregate_months(smooth_map_indicator,new_weights_indicator,area_map,ind_phyto,0.5,0,4);


%get seasonal


% =========================================================================
% Compare the monthly biomes; UHE 15.08.2019
% =========================================================================
n_choice = 9
%     load('CompleSuitePhyto.mat')
%     load('Simple_sort_Data.mat')
load(horzcat('/net/kryo/work/ursho/Damiano_Presence_data/presence_absence_tables_ensemble_averages/Group_specific_background_approach/Data/06Monthly_biomes/PCA_no_mean_v2/Latent_threshold/No_mean_PCA_biomes_',int2str(n_choice),'_v2_5_perc.mat'))
ref_weights = new_weights;
ref_map = smooth_map;%raw_monthly_maps;


new_map = smooth_map_indicator;%raw_monthly_maps_indicator;
removed_weights = new_weights_indicator;

tmp_ref = ref_weights(:,ind_species');

% correspondence part
[metric_1,metric_2,metric_3,overlap_pairs1,overlap_pairs2]  = compare_overlapV2(ref_map,new_map,area_map,tmp_ref,removed_weights);

overlap_maps = ones(12,180,360).*NaN;
for i = 1:length(overlap_pairs1)
    overlap_maps(smooth_map == overlap_pairs1(i,1) & smooth_map_indicator == overlap_pairs2(i,2)) = 1;
end

num_months = sum(~isnan(smooth_map),1);
plotSOM(num_months,1,NaN)
sum_overlap = sum(overlap_maps,1,'omitnan');
plotSOM(sum_overlap./num_months,1,NaN)


%% 
% =========================================================================
% Use smooth_map to create a vector containing the respective classes for
% each pixel; UHE 26 Jul 2019
% =========================================================================
% 
% cd('/net/kryo/work/ursho/Damiano_Presence_data/presence_absence_tables_ensemble_averages/Group_specific_background_approach/Data/11Indicator_species')
% classes_smooth = reverse_prepare2plotV2(No_nan_phyto_simple,smooth_map);
% 
% Smooth_data = No_nan_phyto_simple(~isnan(classes_smooth),5:end-1);
% classes_smooth(isnan(classes_smooth)) = [];
% 
% months = No_nan_phyto_simple(~isnan(classes_smooth),4);
% 
% writematrix(Smooth_data)
% 
% writematrix(classes_smooth)
% writematrix(months)
% %
% % Run R file!!!
% 
% 
% % Read in indval files
% indval=csvread('indval_all.csv',1,1)
% 
% 
% tmp_indval = indval;
% 
% for i = 1:9
%     tmp = indval(:,i);
%     upper = prctile(indval(:,i),99,1)
%     tmp_indval(indval(:,i) < upper,i) = NaN;
% %     jj = input('fdf')
% end
% 
% sum(~isnan(tmp_indval),1)
% 
% a = name_genus_phylum(:,1)
% a = a'

%% Network analysis
%Species networks Table         
% Scores = monthly_Dunning(tmp_raw_monthly_maps,No_nan_phyto_simple,n_clusters);
tic
Scores = monthly_Dunning(corr_corrected_monthly_smooth,Season_obs,n_clusters);
toc

cd('/net/kryo/work/ursho/Damiano_Presence_data/presence_absence_tables_ensemble_averages/Group_specific_background_approach/Data/14Networks')
% save('monthlyScores04Oct2019','Scores','Season_obs');
% save('AcrossmonthScores04Oct2019','Scores','Season_obs');
% save('overallScores04Oct2019','Scores','Season_obs');


% Use the area coverage or core species to make sure the pairs are significant
% for i = 1:8
%     Scores_biome = squeeze(Scores(i,:,:));
%     Scores_biome(Scores_biome <= 0| isnan(Scores_biome)) = 0;
%     writematrix(Scores_biome,horzcat('Scores_biome_',int2str(i),'.csv'))
% end
cd(folder_main)
%% get monthly species pairs


load('CoverageV2.mat')
load('monthlyScores04Oct2019.mat')

size(Scores)

filter_Scores = Scores;

all_biome_pairs = [NaN NaN NaN NaN];
for m = 1:12
tmp_Scores = squeeze(filter_Scores(:,m,:,:));
% tmp_Scores_coverage = squeeze(Scores_coverage(:,m,:,:));
tmp_coverage = squeeze(coverage_month(m,:,:));
% num_biomes = unique(corr_corrected_monthly_smooth(m,~isnan(corr_corrected_monthly_smooth(m,:,:))));
% num_biomes(num_biomes==9) = [];
[biome_pairs] = get_species_pairs(tmp_Scores,tmp_coverage);
biome_pairs
all_biome_pairs = [all_biome_pairs;biome_pairs];
m

end
all_biome_pairs(1,:) = []



[B I] = sort(all_biome_pairs(:,1),'ascend');
sorted_all_biome_pairs = all_biome_pairs(I,:);

% get the network species
final_pairs = [NaN NaN NaN];
for i = 1:8
    tmp = sorted_all_biome_pairs(sorted_all_biome_pairs(:,1) == i,:);
    if(~isempty(tmp))
        connected_edges =  tmp(:,2:3);   
        tmp(:,4)
        [connected_edges,num_pairings] = find_most_connectedV2([tmp(:,2), tmp(:,3)],tmp(:,4));
        num_pairings
        final_pairs = [final_pairs; [connected_edges, connected_edges(:,1).*0+i]];

    end
    
end

final_pairs(1,:) = []
size(final_pairs)


network_species = unique([final_pairs(:,[1,3]);final_pairs(:,[2,3])],'rows')
size(unique(network_species(:,1)))


cd('/net/kryo/work/ursho/Damiano_Presence_data/presence_absence_tables_ensemble_averages/Group_specific_background_approach/Data/14Networks')
% save('Network_species_11Oct2019','network_species','final_pairs')


load('Network_species_11Oct2019.mat')
%% get monthly network species
% load('monthlyScores22Aug2019.mat')
% all_final_pairs = [NaN NaN NaN];
% for m = 1:12
% tmp_Scores = squeeze(Scores(:,m,:,:));
% 
% [final_pairs,network_species,upper,lower] = get_network_speciesV5(tmp_Scores);
% close all
% all_final_pairs = [all_final_pairs;final_pairs];
% % final_pairs
% m
% 
% end
% 
% all_final_pairs(1,:) = [];
% 
% aa = unique(all_final_pairs,'rows')
% 
% [B I] = sort(aa(:,3),'ascend')
% aa(I,:)
% network_species = unique([aa(:,1);aa(:,2)]);
% final_pairs = aa(I,:);
% save('Network_species_22Aug2019','network_species','final_pairs')

cmap = morgenstemning(12)%parula(9);
[cmap_tmp] = shuffle_colormap(cmap);
[cmap_tmp1] = shuffle_colormap(cmap_tmp);


cmap = ametrine(12)%parula(9);
[cmap_tmp] = shuffle_colormap(cmap);
[cmap_tmp2] = shuffle_colormap(cmap_tmp);

cmap = isolum(12)%parula(9);
[cmap_tmp] = shuffle_colormap(cmap);
[cmap_tmp3] = shuffle_colormap(cmap_tmp);

comb_cmap = [cmap_tmp2([2,1,3],:);cmap_tmp1(7,:);cmap_tmp3(5,:);cmap_tmp1(8,:);cmap_tmp1(9,:);cmap_tmp1(4,:)];



copy_combination = final_pairs.*NaN;
correspondence_labels = [NaN NaN];
new_ID = 1;
flag = 0
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
       j = input('fdf')
end
h.NodeLabel = correspondence_labels(:,1);

%% Get niches from model
%plot distribution of the cooccurrences for each biome
Niches_from_model


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
mat_pair = ones(length(final_pairs),8).*NaN;
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

mat_pair_area_sp1 = ones(length(final_pairs),8).*NaN;
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

% place_global = final_pairs(:,1:2).*NaN;
% % find number in list I for each species
% for i = 1:length(final_pairs)
%     [r place_global(i,1)] = find(I == final_pairs(i,1));    
%     [r place_global(i,2)] = find(I == final_pairs(i,2)); 
% end
% 


% 
% 
% Values_per_species(final_pairs(:,1),:)
% Values_per_species(final_pairs(:,2),:)
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
mat_pair_tmp =mat_pair.*NaN;

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

%plot sst for all biomes
% f = figure;
% hold on;
% p = uipanel('Parent',f,'BorderType','none'); 
% %p.Title = 'Normalized Occurrence of PFTs'; 
% p.TitlePosition = 'centertop'; 
% p.FontSize = 12;
% p.FontWeight = 'bold';
% p.BackgroundColor = [0 0 0];
for i = 1:8
    %for each biome i
    matrix_schemaball = ones(length(all_species_net)).*NaN;
    for j = 1:size(mat_pair_tmp,1)
        if mat_pair_tmp(j,i) ~= 0
            [r1 c1] = find(all_species_net == final_pairs(j,1));
            [r2 c2] = find(all_species_net == final_pairs(j,2));
            
            matrix_schemaball(r1,r2) = mat_pair_tmp(j,i);
            matrix_schemaball(r2,r1) = mat_pair_tmp(j,i);
        end
            
    end
    %plot
%     s = subplot(3,3,i,'Parent',p)
%     s.Visible = 'off'
    matrix_schemaball(matrix_schemaball==0) = NaN;
    if i == 2
        schemaballV2(matrix_schemaball,cellstr(num2str(all_species_net)),[0.5 0.5 0.5; 1 0 0],[0, 0, 0],[1 0 0])
    else  
        schemaballV2(matrix_schemaball,cellstr(num2str(all_species_net)),[0.5 0.5 0.5; 1 0 0],[0, 0, 0],[0 0 0])
    end
%     h1 = get(gca,'children')
%     copyobj(h1,s)
end

%%

morgenstemning(2)
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
% name_genus_phylum(final_pairs(:,1),[3,1])
% name_genus_phylum(final_pairs(:,2),[3,1])

%% Species networks SOM



%code stands for type of interaction
% tmp_com = final_pairs(~isnan(final_pairs(:,1)),:);
% coded_phyto = ones(size(No_nan_phyto_simple,1),length(tmp_com))*NaN;
% 
% for i = 1:size(coded_phyto,1)
%     for j = 1:size(coded_phyto,2)
%         %get indeces of species
%         ind1 = tmp_com(j,1)+4;
%         ind2 = tmp_com(j,2)+4;
%         if(No_nan_phyto_simple(i,ind1) == 1 && No_nan_phyto_simple(i,ind2) == 1)
%             coded_phyto(i,j) = 1;
%         elseif(No_nan_phyto_simple(i,ind1) == 0 && No_nan_phyto_simple(i,ind2) == 1)
%             coded_phyto(i,j) = 0.5;
%             
%         elseif(No_nan_phyto_simple(i,ind1) == 1 && No_nan_phyto_simple(i,ind2) == 0)
%             coded_phyto(i,j) = 0.5;
%             
%         elseif(No_nan_phyto_simple(i,ind1) == 0 && No_nan_phyto_simple(i,ind2) == 0)
%             coded_phyto(i,j) = 0;
%         end
%     end
% end
% 
% coded_phyto = [No_nan_phyto_simple(:,1:4),coded_phyto,No_nan_phyto_simple(:,end)];

coded_phyto = No_nan_phyto_simple(:,[1:4,unique(network_species(:,1))'+4,end]);


size(coded_phyto)

% save('Coded_values_21Aug2019_specV2','coded_phyto','-v7.3')
%% Train SOM with coded values 
%load('Coded_values21Nov2018_min5_99.mat')
tic
d1 = 31;
d2 = d1;
optimal_epoch = 200;
[netw_classes, netw_net] = My_SOM2( coded_phyto, d1,d2, optimal_epoch,'mandist' );
toc
cd('/net/kryo/work/ursho/Damiano_Presence_data/presence_absence_tables_ensemble_averages/Group_specific_background_approach/Data/14Networks')
save('Coded_values_SOM_11Oct2019_V3_species','netw_net','netw_classes','coded_phyto','network_species','-v7.3')
cd(folder_main)

%% Construct network biomes

tic
original_weights = netw_net.IW{1};
%original_weights = bsxfun(@minus,original_weights,mean(original_weights));
[coeff,score,latent,tsquared,explained] = pca(original_weights);

[r,c] = find(latent > 1)
clearvars latent
%[coeff,score,latent,tsquared,explained]
[coeff,score,~,~,explained] = pca(original_weights,'NumComponents',r(end));
toc

[raw_monthly_maps_network, new_weights_network,~,~] = Calculate_biomesV2(netw_net, netw_classes,...
    coded_phyto, 9,coeff);

plotSOM(raw_monthly_maps_network,1,9)

[smooth_map_network] = Clean_up_biomesV3( raw_monthly_maps_network,new_weights_network,...
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
load(horzcat('/net/kryo/work/ursho/Damiano_Presence_data/presence_absence_tables_ensemble_averages/Group_specific_background_approach/Data/06Monthly_biomes/PCA_no_mean_v2/Latent_threshold/No_mean_PCA_biomes_',int2str(n_choice),'_v2_5_perc.mat'))
ref_weights = new_weights;
ref_map = smooth_map;%raw_monthly_maps;


new_map = smooth_map_network;%raw_monthly_maps_network;
removed_weights = new_weights_network;

tmp_ref = ref_weights(:,network_species);

% correspondence part
[metric_1,metric_2,metric_3,overlap_pairs1,overlap_pairs2]  = compare_overlapV2(ref_map,new_map,area_map,tmp_ref,removed_weights);

overlap_maps = ones(12,180,360).*NaN;
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



%% Construct seasonal and annual biomes using indicator SOM and network SOM
%compare the overlap using the respective neurons
%compare the overlap only using the spatial extent




Compare_ind_net_to_original





%% Environmental analysis

% Use the respective scripts (Niches_from_model, Get_environmental_data,
% and scripts in folder functions/Environmental









