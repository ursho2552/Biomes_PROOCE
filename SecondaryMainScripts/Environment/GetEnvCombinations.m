%% Preprocessing
%in this script the MLD average values for the environmental data is
%calculated. Then we construct matrices of complete data suites, and look 
%at the similarity between all variables
cd('/net/kryo/work/ursho/PhD/Projects/Biomes/Scripts/Biomes_PROOCE/Data/07Analysis/Environment/')
load('EnvData.mat')

%% MLD average
%need to calculate for nitrate, phosphate, salinity, silicate, and
%temperature

%calculate excess nitrate and excess phosphate
Pstar = phosphate - nitrate./16;
Nstar = nitrate - 16.*phosphate;

%also get only the surface values
nitrate_surf = permute(squeeze(nitrate(:,:,1,:)),[3,2,1]);
phosphate_surf = permute(squeeze(phosphate(:,:,1,:)),[3,2,1]);
silicate_surf = permute(squeeze(silicate(:,:,1,:)),[3,2,1]);
salinity_surf = permute(squeeze(salinity(:,:,1,:)),[3,2,1]);
temperature_surf = permute(squeeze(temperature(:,:,1,:)),[3,2,1]);
Pstar_surf = permute(squeeze(Pstar(:,:,1,:)),[3,2,1]);
Nstar_surf = permute(squeeze(Nstar(:,:,1,:)),[3,2,1]);

MLD_surf = permute(MLD,[3,2,1]);
NPP_surf = permute(NPP,[3,2,1]);
PAR_surf = permute(PAR,[3,2,1]);
chl_surf = permute(chl,[3,2,1]);
pCO2_surf = permute(pCO2,[3,2,1]);
wind_surf = permute(wind,[3,2,1]);


phosphate = permute(phosphate, [3,4,2,1]);
nitrate = permute(nitrate, [3,4,2,1]);


[r100 c100] = find(phosphate_depth == 100);
[r200 c200] = find(phosphate_depth == 200);
P_deep = squeeze(mean(phosphate(c100(1)+1:c200(1),:,:,:),1,'omitnan'));
N_deep = squeeze(mean(nitrate(c100(1)+1:c200(1),:,:,:),1,'omitnan'));
P_surf = squeeze(mean(phosphate(1:c100(1),:,:,:),1,'omitnan'));
N_surf = squeeze(mean(nitrate(1:c100(1),:,:,:),1,'omitnan'));

Delta_P = (P_deep-P_surf);
Delta_N = (N_deep-N_surf);


%% create matrices of complete data suites
%format is: ID, Lon, Lat, month, data1, data2, ..., dataN

%number of env variables (nitrate, phosphate, silicate, N*, P*,
%temperature, salinity, MLD, NPP, PAR, chl, pCO2, wind)
Env_data_month = ones(12,360*180,18)*NaN;
for i =1:12
    Env_data_month(i,:,1:3) = LatLon(:,1:3);
end

%copy each observation to each row
for month = 1:12
    for obs = 1:size(LatLon,1)
        %get indices
        lon = LatLon(obs,4);
        lat = LatLon(obs,5);
       
        Env_data_month(month,obs,4:end) = [month, nitrate_surf(month,lat,lon),...
            phosphate_surf(month,lat,lon), silicate_surf(month,lat,lon), Pstar_surf(month,lat,lon),...
            Nstar_surf(month,lat,lon), salinity_surf(month,lat,lon),temperature_surf(month,lat,lon),...
            MLD_surf(month,lat,lon), NPP_surf(month,lat,lon), PAR_surf(month,lat,lon),...
            chl_surf(month,lat,lon), pCO2_surf(month,lat,lon), wind_surf(month,lat,lon),Delta_P(month,lat,lon)];
    end
end

%look at each month, and for each variable get the rows that are NaN,
%delete them
Complete_env =  ones(1,18)*NaN;
for i = 1:12
    tmp = squeeze(Env_data_month(i,:,:));
    for j = 5:size(Env_data_month,3)
        tmp(isnan(tmp(:,j)),:) = [];
    end
    Complete_env = [Complete_env;tmp];
end
Complete_env(1,:) = [];

cd('/net/kryo/work/ursho/PhD/Projects/Biomes/Scripts/Biomes_PROOCE/Data/07Analysis/Environment/')

if isfile('MonthlyEnv_Oct2019.mat')
    disp('File already exists!')
else
    save('MonthlyEnv_Oct2019','Complete_env','Env_data_month')
end

cd(folder_main)

%% Construct Dendrogram

%construct Dendrogram to see which environmental parameters to use
t = 0.3;

labels_surf = { 'N_{surf}', 'P_{surf}', 'Si_{surf}', 'P*_{surf}', 'N*_{surf}',...
    'SSS','SST', 'MLD', 'log(NPP)', 'PAR', 'log(Chl)', 'pCO_{2}', 'Wind','\Delta P'}


%use only data within the 95%confidence interval of Complete_env

Complete95 = Complete_env;
for i = 5:size(Complete95,2)
    q2_5 = prctile(Complete_env(:,i),2.5,1);
    q97_5 = prctile(Complete_env(:,i),97.5,1);
    Complete95(Complete95(:,i) > q97_5,:) = [];
    Complete95(Complete95(:,i) < q2_5,:) = [];
end
%normalize data to unit variance

Normalized_env = Complete_env;
for i = 5:size(Normalized_env,2)
    mean_tmp = mean(Complete_env(:,i));
    std_tmp = std(Complete_env(:,i));
    for j = 1:length(Normalized_env)
        Normalized_env(j,i) = (Complete_env(j,i)-mean_tmp)/std_tmp;
    end
end

Standardized_env = Complete_env;

%only keep one, mld or surf
Surf_env = Standardized_env;

[ R_spearman_annual_env, orig_annual_env,T_annual_env, EVA1, EVA2, EVA3 ] =...
    get_dendrogram('spearman',Surf_env(:,5:end)',labels_surf,2:12,1,'weighted');

