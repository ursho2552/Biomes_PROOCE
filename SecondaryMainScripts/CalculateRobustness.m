%% Robustness test of end product (seasonally corrected smooth months)

% =========================================================================
% Here we test the robustness of our biomes to information loss or feature 
% loss. This is a stand-alone script which is called from the main script
%  UHE 05/12/2019
% =========================================================================

%load area map
cd('/net/kryo/work/ursho/Damiano_Presence_data/presence_absence_tables_ensemble_averages/Group_specific_background_approach/Data/')
load('Area_map.mat')

%get monthly biomes with full data
cd('/net/kryo/work/ursho/Damiano_Presence_data/presence_absence_tables_ensemble_averages/Group_specific_background_approach/Data/05Biomes/')
load('No_mean_PCA_biomes_9_v2_5_perc.mat')

if isfile('Seasonally_corrected_monthly_biomes.mat')
    load('Seasonally_corrected_monthly_biomes.mat')
else
    %get seasonally corrected version
    corrected_monthly_raw = ones(12,180, 360)*NaN;
    for i =1:12
        if(i<7)
            j = mod(i+6,13);
        else
            j = mod(i+6,13) +1;
        end
        %i and j are the indices of thmatlab e months that need to be combined
        corrected_monthly_raw(i,91:end,:) = raw_monthly_maps(i,91:end,:);
        corrected_monthly_raw(i,1:90,:) = raw_monthly_maps(j,1:90,:); 
    end

    %Smooth seasonally corrected biomes
    [corrected_monthly_smooth] = Clean_up_biomes( corrected_monthly_raw,new_weights,...
    area_map,0.5,4,0);
    %save seasonally corrected biomes
    save('Seasonally_corrected_monthly_biomes','corrected_monthly_smooth')
end


%Construct seasonally corrected biomes for spatial/temporal loss and
%feature loss
cd('/net/kryo/work/ursho/Damiano_Presence_data/presence_absence_tables_ensemble_averages/Group_specific_background_approach/Data/06Robustness/SpatialTemporalLoss/Biomes/')
fr = [1, 5, 10, 20, 30];
for ii = 1:length(fr)
    
    load(horzcat('Leaky_biomes_9_fr_',int2str(fr(ii)),'_5_perc.mat'))
    corrected_monthly_raw = ones(12,180, 360)*NaN;
    for i =1:12
        if(i<7)
            j = mod(i+6,13);
        else
            j = mod(i+6,13) +1;
        end
        %i and j are the indices of thmatlab e months that need to be combined
        corrected_monthly_raw(i,91:end,:) = raw_monthly_maps(i,91:end,:);
        corrected_monthly_raw(i,1:90,:) = raw_monthly_maps(j,1:90,:); 
    end

    %get smooth version
    [Lcorrected_monthly_smooth] = Clean_up_biomesV3( corrected_monthly_raw,new_weights,...
    area_map,0.5,4,0);
    save(horzcat('Leaky_corrected_month_fr_',int2str(fr(ii))),'Lcorrected_monthly_smooth')
end
    
    
%Now for feature loss expriments

cd('/net/kryo/work/ursho/Damiano_Presence_data/presence_absence_tables_ensemble_averages/Group_specific_background_approach/Data/06Robustness/FeatureLoss/Biomes/')
fr = [1, 5, 10, 20, 30];

for ii = 1:length(fr)
    
    load(horzcat('Removed_biomes_9_',int2str(fr(ii)),'_5_perc.mat'))
    corrected_monthly_raw = ones(12,180, 360)*NaN;
    for i =1:12
        if(i<7)
            j = mod(i+6,13);
        else
            j = mod(i+6,13) +1;
        end
        %i and j are the indices of thmatlab e months that need to be combined
        corrected_monthly_raw(i,91:end,:) = raw_monthly_maps(i,91:end,:);%raw_monthly_maps(i,91:end,:);%smooth_map(i,91:end,:);
        corrected_monthly_raw(i,1:90,:) = raw_monthly_maps(j,1:90,:); % raw_monthly_maps(j,1:90,:); %smooth_map(j,1:90,:); % 
    end

    %get smooth version
    [Rcorrected_monthly_smooth] = Clean_up_biomesV3( corrected_monthly_raw,new_weights,...
    area_map,0.5,4,0);
    save(horzcat('Removed_corrected_month_fr_',int2str(fr(ii))),'Rcorrected_monthly_smooth')
end

%% Kappa analysis 

% =========================================================================
% Here we use the area-weighted Kappa index to assess the robustness of our
% biomes to spatio-temporal information loss and feature loss
% =========================================================================
cd('/net/kryo/work/ursho/Damiano_Presence_data/presence_absence_tables_ensemble_averages/Group_specific_background_approach/Data/05Biomes/')
load('Seasonally_corrected_monthly_biomes.mat')

cd('/net/kryo/work/ursho/Damiano_Presence_data/presence_absence_tables_ensemble_averages/Group_specific_background_approach/Data/06Robustness/FeatureLoss/Biomes/')
load('Removed_corrected_month_fr_1.mat')
R1 = Rcorrected_monthly_smooth;
load('Removed_corrected_month_fr_5.mat')
R5 = Rcorrected_monthly_smooth;
load('Removed_corrected_month_fr_10.mat')
R10 = Rcorrected_monthly_smooth;
load('Removed_corrected_month_fr_20.mat')
R20 = Rcorrected_monthly_smooth;
load('Removed_corrected_month_fr_30.mat')
R30 = Rcorrected_monthly_smooth;

cd('/net/kryo/work/ursho/Damiano_Presence_data/presence_absence_tables_ensemble_averages/Group_specific_background_approach/Data/06Robustness/SpatialTemporalLoss/Biomes/')
load('Leaky_corrected_month_fr_1.mat')
L1 = Lcorrected_monthly_smooth;
load('Leaky_corrected_month_fr_5.mat')
L5 = Lcorrected_monthly_smooth;
load('Leaky_corrected_month_fr_10.mat')
L10 = Lcorrected_monthly_smooth;
load('Leaky_corrected_month_fr_20.mat')
L20 = Lcorrected_monthly_smooth;
load('Leaky_corrected_month_fr_30.mat')
L30 = Lcorrected_monthly_smooth;
%%
%kappa analysis for spatio-temporal information loss
for i = 1:5
    kappas = ones(1,12).*NaN;
    switch i
        case 1
            data = L1;
        case 2
            data = L5;
        case 3
            data = L10;
        case 4
            data = L20;
        case 5
            data = L30;
    end


    [correspondence, maps] = get_correspondence_partitioningsV2( corrected_monthly_smooth,data,area_map);

    changed_data = data.*NaN;
    for m = 1:size(correspondence,1)
        changed_data(data == correspondence(m,2)) = correspondence(m,1);
    end


    for m = 1:12

        kappas(m) = cohensKappa(corrected_monthly_smooth(m,:,:),changed_data(m,:,:));

    end
    disp('Res:')
    mean(kappas)
    std(kappas)
    jj = input('next?')
end

%kappa analysis for feature loss
for i = 1:5
    kappas = ones(1,12).*NaN;
    switch i
        case 1
            data = R1;
        case 2
            data = R5;
        case 3
            data = R10;
        case 4
            data = R20;
        case 5
            data = R30;
    end


    [correspondence, maps] = get_correspondence_partitioningsV2( corrected_monthly_smooth,data,area_map);

    changed_data = data.*NaN;
    for m = 1:size(correspondence,1)
        changed_data(data == correspondence(m,2)) = correspondence(m,1);
    end


    for m = 1:12

        kappas(m) = cohensKappa(corrected_monthly_smooth(m,:,:),changed_data(m,:,:));

    end
    disp('Res:')
    mean(kappas)
    std(kappas)
    jj = input('next?')
end
