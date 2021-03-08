%% Setup folder
clear all;

folder_main = '/net/kryo/work/ursho/Damiano_Presence_data/presence_absence_tables_ensemble_averages/Group_specific_background_approach';
addpath(genpath(folder_main))

cd(folder_main)


%% Read in data

% =========================================================================
% This part takes a while. Use the scripts RData2CSV and
% CSV2MAT.py instead. These are found under "./Scripts/functions/Read_Data/"
% Alternatively you can create a script to import the files in Matlab by
% double clicking the csv file, marking the section you need, and selecting
% "create script" option
% =========================================================================

names = {'Jan_GAM','Feb_GAM','Mar_GAM','Apr_GAM','May_GAM','Jun_GAM','Jul_GAM',...
    'Aug_GAM','Sep_GAM','Oct_GAM','Nov_GAM','Dec_GAM'};
dir = '/net/kryo/work/ursho/Damiano_Presence_data/presence_absence_tables_ensemble_averages/Group_specific_background_approach/Data/00Probabilities/GAM/';
for i=1:12
    %define the name of the file
    str = horzcat('dfs_month_',int2str(i),'_trsp_gam_pa_gridded_gr_bg(sit_ove).csv');
    %define variable name
    varname = genvarname(names(i));
    %read data, importfile1() was created by matlab
    data_read = importfile1(str);
    %assign data to varname
    eval([varname{i} '=data_read;'])
    %save file to the directory, change if not the same as the data
    save(horzcat(dir,varname{i}),varname{i})

end

% =========================================================================
% Import header, i.e. the column names. Here is the same as for the data
% itself, either use python or with a matlab script, or manually
% =========================================================================

labels_phyto = importfile2('dfs_month_1_trsp_gam_pa_gridded_gr_bg(sit_ove).csv');

% =========================================================================
% load data saved above
% =========================================================================

load('/net/kryo/work/ursho/Damiano_Presence_data/presence_absence_tables_ensemble_averages/Group_specific_background_approach/Data/00Probabilities/GAM/Jan_GAM.mat')
load('/net/kryo/work/ursho/Damiano_Presence_data/presence_absence_tables_ensemble_averages/Group_specific_background_approach/Data/00Probabilities/GAM/Feb_GAM.mat')
load('/net/kryo/work/ursho/Damiano_Presence_data/presence_absence_tables_ensemble_averages/Group_specific_background_approach/Data/00Probabilities/GAM/Mar_GAM.mat')
load('/net/kryo/work/ursho/Damiano_Presence_data/presence_absence_tables_ensemble_averages/Group_specific_background_approach/Data/00Probabilities/GAM/Apr_GAM.mat')
load('/net/kryo/work/ursho/Damiano_Presence_data/presence_absence_tables_ensemble_averages/Group_specific_background_approach/Data/00Probabilities/GAM/May_GAM.mat')
load('/net/kryo/work/ursho/Damiano_Presence_data/presence_absence_tables_ensemble_averages/Group_specific_background_approach/Data/00Probabilities/GAM/Jun_GAM.mat')
load('/net/kryo/work/ursho/Damiano_Presence_data/presence_absence_tables_ensemble_averages/Group_specific_background_approach/Data/00Probabilities/GAM/Jul_GAM.mat')
load('/net/kryo/work/ursho/Damiano_Presence_data/presence_absence_tables_ensemble_averages/Group_specific_background_approach/Data/00Probabilities/GAM/Aug_GAM.mat')
load('/net/kryo/work/ursho/Damiano_Presence_data/presence_absence_tables_ensemble_averages/Group_specific_background_approach/Data/00Probabilities/GAM/Sep_GAM.mat')
load('/net/kryo/work/ursho/Damiano_Presence_data/presence_absence_tables_ensemble_averages/Group_specific_background_approach/Data/00Probabilities/GAM/Oct_GAM.mat')
load('/net/kryo/work/ursho/Damiano_Presence_data/presence_absence_tables_ensemble_averages/Group_specific_background_approach/Data/00Probabilities/GAM/Nov_GAM.mat')
load('/net/kryo/work/ursho/Damiano_Presence_data/presence_absence_tables_ensemble_averages/Group_specific_background_approach/Data/00Probabilities/GAM/Dec_GAM.mat')

% =========================================================================
% Remove header/column descriptor. Check prior to deleting!
% =========================================================================

Jan = Jan_GAM(2:end,:);
Feb = Feb_GAM(2:end,:);
Mar = Mar_GAM(2:end,:);
Apr = Apr_GAM(2:end,:);
May = May_GAM(2:end,:);
Jun = Jun_GAM(2:end,:);
Jul = Jul_GAM(2:end,:);
Aug = Aug_GAM(2:end,:);
Sep = Sep_GAM(2:end,:);
Oct = Oct_GAM(2:end,:);
Nov = Nov_GAM(2:end,:);
Dec = Dec_GAM(2:end,:);

% =========================================================================
% Stack months in one matrix
% =========================================================================

All_phyto = [Jan;Feb;Mar;Apr;May;Jun;Jul;Aug;Sep;Oct;Nov;Dec];

% =========================================================================
% Delete first column to define an ID. Check if this is still needed!!
% =========================================================================

All_phyto(:,1) = [];
All_phyto(:,1) = 1:size(All_phyto,1);
labels_phyto(1) = [];

% =========================================================================
% Only consider complete columns, i.e. go through all columns, if a whole 
% column is NaN, delete it
% =========================================================================

n = size(All_phyto,1);
r = NaN;
for i = 1:size(All_phyto,2)
    tmp = sum(isnan(All_phyto(:,i)));
    
    if(tmp == n)
        r = [r;i];
    end
    
end
r(1) = [];
All_phyto(:,r) = [];
labels_phyto(r) = [];

% =========================================================================
% Only consider complete rows, i.e. delete rows where we find a NaN for any
% species (i.e. only complete communities allowed)
% =========================================================================

for i =1:size(All_phyto,2)
    tmp = isnan(All_phyto(:,i));
    All_phyto(tmp == 1,:) = [];
end

% =========================================================================
% All_phyto now has no NaNs and the following structure
% ID|Lon|Lat|month|species1....n|buffer. Currently, the buffer space is 
% not used.
% =========================================================================

No_nan_phyto = [All_phyto(:,[1 3 4]),All_phyto(:,2),...
    All_phyto(:,5:end), (1:length(All_phyto))'*NaN];

% =========================================================================
% Save your complete, merged data (commented out to avoid overwritting)
% =========================================================================

%got to data folder
cd('/net/kryo/work/ursho/Damiano_Presence_data/presence_absence_tables_ensemble_averages/Group_specific_background_approach/Data/00Probabilities')
%save file
% save('CompleSuitePhyto','No_nan_phyto','labels_phyto')
%return to working folder
cd(folder_main)


%% Transform to binary (1/0) using threshold 0.5 (strictly above 0.5 is 1 and below or equal to 0.5 is 0)

% =========================================================================
% The values in No_nan_phyto are ensembles of presences. If more than
% half of the number of ensemble members project a presence we mark it as a
% presence in the ensemble.
% =========================================================================

%load the data
cd('/net/kryo/work/ursho/Damiano_Presence_data/presence_absence_tables_ensemble_averages/Group_specific_background_approach/Data/00Probabilities')
load('CompleSuitePhyto.mat')
%allocate new matrix
Transformed_phyto = No_nan_phyto;
%change the values according to your threshold
Transformed_phyto(No_nan_phyto <= 0.5) = 0;
Transformed_phyto(No_nan_phyto > 0.5) = 1;
%since you changed the whole matrix (including ID, lon, lat, month) you
%need to overwrite these four fields
Transformed_phyto(:,1:4) = No_nan_phyto(:,1:4);
%save the data
if isfile('Transformed_CompleteSuitePhyto.mat')
    disp('File already exists!')
else
    save('Transformed_CompleteSuitePhyto','Transformed_phyto','labels_phyto')
end

cd(folder_main)

%% Construct helper matrix for plotting

% =========================================================================
% Here you define a helper matrix that speeds up plotting and other
% sorting problems, so it is worth spending a few seconds on this to save
% minutes later
% =========================================================================

% =========================================================================
% Create matrix "LatLon", with the structure ID|LON|LAT|X-coord|Y-coord.
% LON starts at -179.5 and goes to 179.5, LAT starts as -89.5 and goes to
% 89.5, X- and Y-coord refer to the index (1 based) in a 2D matrix
% =========================================================================

%allocate matrix
LatLon = NaN(180*360,5);
%define ID
LatLon(:,1) = 1:size(LatLon,1);
lons = -179.5:179.5;
lats = -89.5:89.5;
%c is a dummy counter, which increases with each iteration
c = 1;
for i = 1:length(lons)
    for j = 1:length(lats)
        LatLon(c,2) = lons(i);
        LatLon(c,3) = lats(j);
        LatLon(c,4) = i;
        LatLon(c,5) = j;
        c = c + 1;
    end
end

% =========================================================================
% Save LatLon and lons, lats for possible use later
% =========================================================================

cd('/net/kryo/work/ursho/Damiano_Presence_data/presence_absence_tables_ensemble_averages/Group_specific_background_approach/Data/')
if isfile('HelpVariables.mat')
    disp('File already exists!')
else
    save('HelpVariables','LatLon','lons','lats','labels_phyto')
end

% =========================================================================
% Use LatLon to substitute lat lon for X- and Y-indeces in the merged,
% complete dataset
% =========================================================================

%allocate data
No_nan_phyto_simple = Transformed_phyto;
%loop over each latlon combination, and search for all instances in the
%dataset
for i = 1:length(LatLon)
    [r, ~] = find(Transformed_phyto(:,2) == LatLon(i,2) & Transformed_phyto(:,3) == LatLon(i,3));
    for j =1:length(r)
        No_nan_phyto_simple(r(j),2:3) = LatLon(i,4:5);
    end
end
%save the dataset
cd('/net/kryo/work/ursho/Damiano_Presence_data/presence_absence_tables_ensemble_averages/Group_specific_background_approach/Data/00Probabilities')
if isfile('Simple_sort_Data.mat')
    disp('File already exists!')
else
    save('Simple_sort_Data','No_nan_phyto_simple')
end
cd(folder_main)


%% Get species, genus, and phylum of species available in the complete data

% =========================================================================
% Make sure you have all species, genera that you actually intend to have
% =========================================================================

%load data
cd('/net/kryo/work/ursho/Damiano_Presence_data/presence_absence_tables_ensemble_averages/Group_specific_background_approach/Data/00Probabilities')
load('Transformed_CompleteSuitePhyto.mat')

%read in table with the following information: Species name|Species group 
%(e.g. diatoms) as columns, and rows for the modeled species
%You can read in and save the data manually by double-clicking it and
%getting in inside matlab, or again using a similar script as CSV2MAT.py
%In this case we read in gam.gr.bg.sit.ove.models_overview.csv (second and third column as
%string array). This matrix contains the names and genus of the
%phytoplankton in the cleaned data
load('models_overview_gam_gr_bg.mat')

%labels_phyto contains the names of columns in your data
%get only the name of species, and get rid of ID, lon, lat, month
tmp_labels_phyto = labels_phyto(5:end);
%search each entry and compare it to the table with species name and genera
for i = 1:length(tmp_labels_phyto)
    %split the name to get genus and species
    genus_species = split(tmp_labels_phyto(i));
    %in the table search for full species name
    [r, ~] = find(gam(:,1) == tmp_labels_phyto(i));
    %store the information in name_genus_phylum as full name of
    %species|genus|group
    name_genus_phylum(i,:) = [tmp_labels_phyto(i) genus_species(1) gam(r,2)];
end

%save the dataset
if isfile('Names_species.mat')
    disp('File already exists!')
else
    save('Names_species','name_genus_phylum')
end

cd(folder_main)
%% Table 1 in manuscript

% =========================================================================
% There are two errors according to algaebase.org. (1) Dictyocha fibula is of the class
% Dictyochophyceae (not Chrysophyceae). (2) Octactis speculum is of the
% class Dictyochophyceae (not Chrysophyceae)
% =========================================================================

%load data
cd('/net/kryo/work/ursho/Damiano_Presence_data/presence_absence_tables_ensemble_averages/Group_specific_background_approach/Data/00Probabilities')
load('Names_species.mat');

%get unique phyla modeled, and print it
unique_phyl = unique(name_genus_phylum(:,3))

for i = 1:length(unique_phyl)
    tmp = name_genus_phylum;
    tmp(name_genus_phylum(:,3) ~= unique_phyl(i),:) = []
    unique_phyl(i)
    num_species = size(tmp,1)
    num_genera = length(unique(tmp(:,2)))
    %this is used to run the loop step-wise with user input. If the console
    %prints "next?" simply press enter
    jj = input('next?')
end

cd(folder_main)

%% Run the Self-Organizing Map

% =========================================================================
% Run the "SOM" file. This file is a stand-alone file that can be called 
% independetly from the main file, i.e. sort of a "secondary main file"
% Running SOM as is will take a very long time. Thus, we instead use the
% bash script "Euler_run_SOM.sh" to run every choice in parallel (Note that
% some of the parameters in "Euler_run_SOM.sh" need to be changed prior to
% runing the script
% =========================================================================

SOM

%% Figure A.9a)


% =========================================================================
% Plot the total error as a function of increasing number of neurons
% =========================================================================

cd('/net/kryo/work/ursho/Damiano_Presence_data/presence_absence_tables_ensemble_averages/Group_specific_background_approach/Data/01Neurons_error')
str1 = 'Single_run_';
str2 = '_error.mat';
nn = 15;
error_neurons = NaN(1,nn)*NaN;
for i = 1:nn
    
        load(horzcat(str1,int2str(i),str2),'total_error')
        error_neurons(i) = total_error;
   
end
    
     
dim1= strcat(int2str(num_neurons(1:nn,1)),' x ',int2str(num_neurons(1:nn,2)));
label = dim1;
n4 = size(num_neurons,1);
h = num_neurons(1:nn,1).*num_neurons(1:nn,2);
x = 1:n4;
y = error_neurons(1:end);
%calculate first derivative
dydx = y.*NaN;

ASF = y*NaN;
%The change should be at the place of the added neurons, not before adding
%them --> easier to interpret!
for i = 2:length(y)-1
        dydx(i) = (y(i+1) - y(i-1))/(h(i+1)-h(i-1));
end

% 5% of the first decrease by increasing the number of neurons
first = dydx(2)*0.05;


figure()
hold on;
xlabel('Number of Neurons')
yyaxis left
plot(h,y)
ylabel('Total error')
yyaxis right
plot(h,-dydx)
ylabel('Total error change')
set(gca,'XTick', h)
set(gca,'XTickLabel',[label])
xtickangle(90)
hline = refline(0,-first);
hline.Color = 'r';
grid on
hold off;

% --> 31 by 31 neurons

%% Figure A.9b)


% =========================================================================
% Plot the total error as a function of increasing number of epochs
% =========================================================================

str1 = 'Single_run_epoch';
str2 = '_error.mat';
nn = 11;
error_neurons_ep = NaN(1,nn)*NaN;
for i = 1:nn
    if(i <= 8)
    
        load(horzcat(str1,int2str(i),str2),'total_error_ep')
    else
        load(horzcat(str1,'_',int2str(i),'.mat'),'total_error_ep')
        
    end
    error_neurons_ep(i) = total_error_ep;
   
end
epoch = [1, 5, 10, 20, 50, 100, 200, 300, 400, 500, 700, 1000];

x = epoch(1:nn);
y = error_neurons_ep(1:nn);

figure()
hold on;
plot(x,y)
xlabel('Epochs')
ylabel('Total error')
grid on
vline = line([200 200], [40 max(error_neurons_ep)]);
vline.Color = 'r';
hold off;


%% Cross-Validation test

% =========================================================================
% We performed an extended Cross-Validation (CV) test of increasing
% leave-out data fraction. We trained SOMS with 90%, 80%, 66%, and 50% of
% the data, while using the remaining 10%, 20%, 33% and 50% for validation.
% This was done such that each data partition was used once for validation
% and 9, 4, 2, 1 for training. 
% =========================================================================

% =========================================================================
% The first step is to create the CV folds. For this we use the stnad-alone
% script "Create_CV_folds
% =========================================================================

CreateCVfolds

% =========================================================================
% After determining the folds, we train SOMs using 90%, 80%, 66%, and 50% 
% of the data. For this we use the bash script "Euler_run_CV.sh" which 
% calls the function CV_SOM with specified parameters to run everything in
% parallel on a cluster
% =========================================================================
% =========================================================================
% =========================================================================
% ===== CALLING THE SCRIPT "Euler_run_CV.sh"IS DONE MANUALLY !!!! =========
% =========================================================================
% =========================================================================
% =========================================================================


% =========================================================================
% After training the SOMs, we validate against the leave-out data using the
% following stand-alone script. The script "Validation" performs the
% dimensionality reduction, clustering, and validation. In the end it plots
% the error of the clustering as a function of increasing number of
% clusters. This has to be adapted to each application as the threshold
% lines in the plots are hardcoded. This script was used to plot Figure 1
% in the manuscript
% =========================================================================

Validation


cd(folder_main)

%% Determine parameters for the definition of biomes

% =========================================================================
% After determining the optimal number of clusters (9), we now cluster our
% optimal setup SOM (trained with the entire dataset), and perform a
% dimensionality reduction of the neurons, cluster them into 9 clusters,
% and map them
% =========================================================================

%load data
cd('/net/kryo/work/ursho/Damiano_Presence_data/presence_absence_tables_ensemble_averages/Group_specific_background_approach/Data/01NeuronsError/')
load('Single_run_11.mat')

%load help variables
cd('/net/kryo/work/ursho/Damiano_Presence_data/presence_absence_tables_ensemble_averages/Group_specific_background_approach/Data/')
load('HelpVariables.mat')

if isfile('Area_map.mat')
    load('Area_map.mat')
else
    %calculate an area map for area weighted calculations
    area = NaN(size(LatLon,1),3);
    area(:,1:2) = LatLon(:,2:3);
    %LatLon (2 is Lon and 3 is Lat), for every combination get the area in km
    earthellipsoid = referenceSphere('earth','km');
    for i = 1:size(LatLon,1)
        lat_tmp = LatLon(i,3);
        lon_tmp = LatLon(i,2);
        lat = [lat_tmp+0.5, lat_tmp-0.5];
        lon = [lon_tmp+0.5, lon_tmp-0.5];

        area(i,3) = areaquad(lat(1),lon(1),lat(2),lon(2),earthellipsoid);
    end
    area_map = prepare2plot([LatLon(:,4:5),LatLon(:,4)*0+1,area(:,3)]);
    %save the map of pixel area for the future
    save('Area_map','area_map')
end

% =========================================================================
% Now we  test our definition of centroid, i.e. how are the
% centroids best defined, (1) weighted by the occurrence frequency of each
% memeber or (2) simply using the unique members. For this we use the
% stand-alone script CompareCentroidMethod, which plots Figure A.10 in the
% manuscript.
% =========================================================================

CompareCentroidMethod

%The once that consider the occurrence frequency perform better

% =========================================================================
% Now we test how small biomes can be without loosing a lot of information
% =========================================================================

%perform dimensionality reduction on trained neurons
tic
original_weights = net.IW{1};
%original_weights = bsxfun(@minus,original_weights,mean(original_weights));
[coeff,score,latent,tsquared,explained] = pca(original_weights);


%Use Kaiser's rule, i.e. dismiss principal components that explain less
%variance than a single original variable
[r,c] = find(latent > 1)

%print variance explained
round(explained(r(1:end))*100)/100

%print cumulative variance explained
sum(explained(r(1:end)))
clearvars latent
%[coeff,score,latent,tsquared,explained]
[coeff,score,~,~,explained] = pca(original_weights,'NumComponents',r(end));
toc


% =========================================================================
% REVIEW 23/11/2020 loadings START
% =========================================================================
max(coeff,[],1)
min(coeff,[],1)
% =========================================================================
% REVIEW 23/11/2020 loadings END
% =========================================================================

cd('/net/kryo/work/ursho/Damiano_Presence_data/presence_absence_tables_ensemble_averages/Group_specific_background_approach/Data/05Biomes/PotentialBiomes/')
%calculate from 2 to 20 biomes
n_clusters = [2:20];
%define how small biomes should be, 0.1% of the global ocean area up to 10%
fracs = [0.1:0.1:2, 3:10];
for j = 1:length(fracs)
    frac = fracs(j);
    for i = 1:length(n_clusters)

        [raw_monthly_maps, new_weights,~,~] = Calculate_biomes(net, classes,...
            No_nan_phyto_simple, n_clusters(i),coeff);
        n = n_clusters(i);


        [smooth_map] = Clean_up_biomes( raw_monthly_maps,new_weights,...
        area_map,frac,4,0);


         save(horzcat('No_mean_PCA_biomes_',int2str(n_clusters(i)),'_v2_',int2str(frac*10),'_perc'),...
             'raw_monthly_maps','new_weights','smooth_map','n','frac')
    end
end

% =========================================================================
% Check the number of pixels changed from raw
% =========================================================================

[neuron_maps] = prepare2plot( [No_nan_phyto_simple(:,2:4),classes]);
orig_neurons = net.IW{1};
fracs = [0.1:0.1:2, 3:10];
frac_change = fracs.*NaN;

for ii = 1:length(fracs)
    tmp = NaN;
    for i = 1:length(n_clusters)
        load(horzcat('No_mean_PCA_biomes_',int2str(n_clusters(i)),'_v2_',int2str(fracs(ii)*10),'_perc.mat'));

        %for each label in raw and smooth get the difference to the
        %original neurons, then subtract the raw from smooth
        tmp_new_weights = [new_weights;new_weights(1,:).*0]; %add dummy weight
        raw_tmp = raw_monthly_maps;
        smooth_tmp = smooth_map;
        raw_tmp(isnan(raw_monthly_maps)) = 0;
        smooth_tmp(isnan(smooth_map)) = 0;

        diff_map = raw_tmp-smooth_tmp;
        diff_map(diff_map ~= 0) = 1;

        %calculate difference to neurons at places where diff_map is 1

        orig_raw = orig_neurons(neuron_maps(diff_map == 1),:);
        orig_smooth = orig_neurons(neuron_maps(diff_map == 1),:);

        tmp_new = orig_raw.*NaN;
        tmp_new(:,:) = new_weights(raw_monthly_maps(diff_map == 1),:);
        raw_dist =mean(sum(abs(orig_raw-tmp_new),2));

        tmp_new = orig_smooth.*NaN;
        smooth_tmp = smooth_map;
        smooth_tmp(isnan(smooth_map)) = size(tmp_new_weights,1);
        tmp_new(:,:) = tmp_new_weights(smooth_tmp(diff_map == 1),:);
        smooth_dist =mean(sum(abs(orig_smooth-tmp_new),2));

        tmp = [tmp;(smooth_dist-raw_dist)/smooth_dist];
        
    end
    frac_change(ii)  = mean(tmp(2:end));
end

figure
hold on
plot(fracs,frac_change)
grid on


figure
hold on
plot(fracs(1:20),frac_change(1:20))
grid on


%upper threshold at 2%, but 0.5% and 1.8% as critical changing points,
%using function findchangepts with 2 MaxNumChanges, since 1 does not find a
%changing point!!! --> lower likely better, since less added error

cd(folder_main)
%% Producing biomes

% =========================================================================
% We now produce biomes based on the choices tested above, i.e. 9 clusters
% are optimal, 0.5% of the global ocean area is the minimum area for a
% biome, and cluster centroids are weighted by the occurrence frequency of
% each cluster member.
% =========================================================================
cd('/net/kryo/work/ursho/Damiano_Presence_data/presence_absence_tables_ensemble_averages/Group_specific_background_approach/Data/05Biomes/PotentialBiomes/')

load('No_mean_PCA_biomes_9_v2_5_perc.mat')
[smooth_annual_map, annual_map] = aggregate_months(smooth_map,new_weights,area_map,No_nan_phyto_simple,0.5,0,4);
[uncertainty_smooth_annual_map, uncertainty_annual_map] = aggregate_months(smooth_map,new_weights,area_map,No_nan_phyto_simple,0.5,1,4);

cd('/net/kryo/work/ursho/Damiano_Presence_data/presence_absence_tables_ensemble_averages/Group_specific_background_approach/Data/05Biomes/')

save('No_mean_PCA_biomes_annual_9_v2_5_perc','smooth_annual_map','annual_map',...
    'uncertainty_smooth_annual_map','uncertainty_annual_map')

%plot the annual biome to see the structure
plotSOM(smooth_annual_map,1,9)


cd(folder_main)

%% Testing the robustness of our biomes

% =========================================================================
% In this section we test the sensitivity of our optimal SOM to information
% loss. First, we test the sensitivity of our the SOM is to loss of 
% spatial/temporal information (individual pixels). Second, we test the 
% sensitivity of our SOM to loss of features (species across the dataset).
% =========================================================================

%load helper variables
cd('/net/kryo/work/ursho/Damiano_Presence_data/presence_absence_tables_ensemble_averages/Group_specific_background_approach/Data/')
load('HelpVariables.mat')

%load simple sort data
cd('/net/kryo/work/ursho/Damiano_Presence_data/presence_absence_tables_ensemble_averages/Group_specific_background_approach/Data/00Probabilities/')
load('Simple_sort_Data.mat')


cd('/net/kryo/work/ursho/Damiano_Presence_data/presence_absence_tables_ensemble_averages/Group_specific_background_approach/Data/06Robustness/SpatialTemporalLoss/')
%create matrix indicating which obsrvations should be deleted
noisy_data_ind = No_nan_phyto_simple(:,1:30).*NaN;
nsamples = size(noisy_data_ind,1);

for fr = 1:30
    nsamples_fr = round(nsamples*fr/100);
    r_ind = randsample(nsamples,nsamples_fr);
    noisy_data_ind(r_ind,fr) = fr;
end
%save noisy data indeces
if isfile('Noisy_data_ind.mat')
    disp('File already exists!')
else
    save('Noisy_data_ind','noisy_data_ind')
end
% =========================================================================
% Now we again run a bash script "Euler_run_Robustness.sh to train 10
% different SOMs. This again has to be done manually!!!!
% =========================================================================

% =========================================================================
% ================= RUN Euler_run_Robustness.sh ===========================
% =========================================================================


% =========================================================================
% After running the scripts and having the trained SOMs, we now produce
% biomes for each robustness experiment. First we produce biomes for the
% spatio-temporal information loss experiment
% =========================================================================

%load area map
cd('/net/kryo/work/ursho/Damiano_Presence_data/presence_absence_tables_ensemble_averages/Group_specific_background_approach/Data/')
load('Area_map.mat')

fr = [1 5 10 20 30];
for i =1:5
    %load network
    cd('/net/kryo/work/ursho/Damiano_Presence_data/presence_absence_tables_ensemble_averages/Group_specific_background_approach/Data/06Robustness/SpatialTemporalLoss/SOMs/')
    str1 = horzcat('Leaky_SOM_frac_',int2str(fr(i)))
    load(horzcat(str1,'.mat'))

    tic
    [ classes_noise_all ] = Apply_SOM( No_nan_phyto_simple, noise_net );
    toc
    
    n_clusters = 9;
  
    tic
    original_weights = noise_net.IW{1};
    %original_weights = bsxfun(@minus,original_weights,mean(original_weights));
    [~,~,latent,~,~] = pca(original_weights);

    [r,~] = find(latent > 1);
    clearvars latent
    [coeff,~,~,~,~] = pca(original_weights,'NumComponents',r(end));
    toc 
  
    [raw_monthly_maps, new_weights,~,~] = Calculate_biomes(noise_net, classes_noise_all,...
        No_nan_phyto_simple, n_clusters,coeff);
    n = n_clusters;

    [smooth_map] = Clean_up_biomes( raw_monthly_maps,new_weights,...
    area_map,0.5,4,0);

    cd('/net/kryo/work/ursho/Damiano_Presence_data/presence_absence_tables_ensemble_averages/Group_specific_background_approach/Data/06Robustness/SpatialTemporalLoss/Biomes/')
    save(horzcat('Leaky_biomes_',int2str(n_clusters),'_fr_',int2str(fr(i)),'_5_perc'),...
    'raw_monthly_maps','new_weights','smooth_map','n','classes_noise_all')
      
end

% =========================================================================
% Now we produce biomes for the feature loss experiment
% =========================================================================

fr = [1 5 10 20 30];
for i =length(fr):-1:1
    %load network
    cd('/net/kryo/work/ursho/Damiano_Presence_data/presence_absence_tables_ensemble_averages/Group_specific_background_approach/Data/06Robustness/FeatureLoss/SOMs/')
    str1 = horzcat('Removed_SOM_frac_',int2str(fr(i)));
    load(horzcat(str1,'.mat'))
    
    n_clusters = 9;
    
    tic
    original_weights = removed_net.IW{1};
    [~,~,latent,~,~] = pca(original_weights);

    [r,~] = find(latent > 1);
    clearvars latent
    [coeff,~,~,~,~] = pca(original_weights,'NumComponents',r(end));
    toc
    
    [raw_monthly_maps, new_weights,~,~] = Calculate_biomes(removed_net, removed_classes,...
        tmp_phyto, n_clusters,coeff);
    n = n_clusters;

    [smooth_map] = Clean_up_biomes( raw_monthly_maps,new_weights,...
    area_map,0.5,4,0);

    cd('/net/kryo/work/ursho/Damiano_Presence_data/presence_absence_tables_ensemble_averages/Group_specific_background_approach/Data/06Robustness/FeatureLoss/Biomes/')
    save(horzcat('Removed_biomes_',int2str(n_clusters),'_',int2str(fr(i)),'_5_perc'),...
    'raw_monthly_maps','new_weights','smooth_map','n','tmp_phyto')
    
end


% =========================================================================
% Validation test yields 9 clusters, robustness shows around nine clusters,
% the change caused by adding more is not "significant", despite the area
% coverage being lower than for lower clusterings (we want the one that is
% stable), Testing for fraction to decide on smallest biome size is 0.5% as
% seen in the plot, on the lower tail, the changing point is at 0.5 and 1.8
% and we want to maintain most of the structures, thus 0.5% is our choice.
% =========================================================================

%% Compare the robustness of biomes (Table 2)

% =========================================================================
% Here we use the stand-alone script "CalculateRobustness" to first produce
% seasonally corrected monthly biomes for our trained SOM, and for the SOMs
% from the robustness experiments. Then we use the area-weighted Kappa
% index to calculate the overlap between the seasonally corrected biomes.
% This script was used to fill in Table 2 of the manuscript.
% =========================================================================

CalculateRobustness


%% Calculate seasonal biomes

% =========================================================================
% To finalize all preparations prior to the analysis, we now produce 
% seasonal biomes
% =========================================================================

%load monthly and annual biomes
cd('/net/kryo/work/ursho/Damiano_Presence_data/presence_absence_tables_ensemble_averages/Group_specific_background_approach/Data/05Biomes/')
load('No_mean_PCA_biomes_9_v2_5_perc.mat')
load('No_mean_PCA_biomes_annual_9_v2_5_perc.mat')


% =========================================================================
% In order to create seasonal biomes, we first need to seasonally correct
% the data set by season, i.e. the month column in 
% No_nan_phyto_simple dataset needs to be shifted by 6 months in the
% Southern Hemisphere. This way, when speaking of biomes in fall, we can
% show the location of biomes during austral and boreal fall in the same
% map.
% =========================================================================

%Create a map of the observation IDs, i.e. the rows

%load data
cd('/net/kryo/work/ursho/Damiano_Presence_data/presence_absence_tables_ensemble_averages/Group_specific_background_approach/Data/00Probabilities/')
load('Simple_sort_Data.mat')

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
    corrected_monthly_raw(i,91:end,:) = raw_monthly_maps(i,91:end,:)
    corrected_monthly_raw(i,1:90,:) = raw_monthly_maps(j,1:90,:); 
    
    corrected_monthly_ID(i,91:end,:) = ID_maps(i,91:end,:);
    corrected_monthly_ID(i,1:90,:) = ID_maps(j,1:90,:);

end

% =========================================================================
% The matrix "corrected_monthly_ID" now contains monthly maps of
% observation IDs in the seasonally corrected months. The same applies for
% "corrected_monthly_raw" but for the biome labels.
% =========================================================================

%Smooth the new seasonally corrected monthly biomes
[corrected_monthly_smooth] = Clean_up_biomes( corrected_monthly_raw,new_weights,...
area_map,0.5,4,0);

%reverse the map-data to matrix
tic
[Season_No_nan_phyto_simple] = reverse_prepare2plot(corrected_monthly_ID);
toc

cd('/net/kryo/work/ursho/Damiano_Presence_data/presence_absence_tables_ensemble_averages/Group_specific_background_approach/Data/00Probabilities/')
if isfile('Seasonally_corrected_data.mat')
    disp('File already exists!')
else
    save('Seasonally_corrected_data','Season_No_nan_phyto_simple')
end

% =========================================================================
% Using the seasonally corrected monthly maps, and seasonally corrected
% dataset (matrix) define seasonal maps
% =========================================================================

Season_obs = [Season_No_nan_phyto_simple(:,end), Season_No_nan_phyto_simple(:,1:3),...
    No_nan_phyto_simple(Season_No_nan_phyto_simple(:,end),5:end)];
seasons = [3 4 5;6 7 8;9 10 11;12,1,2];
season_map = NaN(4,180,360);
season_map_smooth = season_map;
uncertainty_smooth_seasonal_map = season_map;
uncertainty_season_map = season_map;


for s = 1:4
    tmp_obs = Season_obs(1,:).*NaN;
    for i = 1:3
        tmp_obs = [tmp_obs;Season_obs(Season_obs(:,4) == seasons(s,i),:)];
    end
    tmp_obs(1,:) = [];
    [season_map_smooth(s,:,:), season_map(s,:,:)] = aggregate_months(corrected_monthly_smooth(seasons(s,:),:,:),new_weights,area_map,tmp_obs,0.5,0,4);
    [uncertainty_smooth_seasonal_map(s,:,:), uncertainty_season_map(s,:,:)] = aggregate_months(corrected_monthly_smooth(seasons(s,:),:,:),new_weights,area_map,tmp_obs,0.5,1,4);
    s
end
 

cd('/net/kryo/work/ursho/Damiano_Presence_data/presence_absence_tables_ensemble_averages/Group_specific_background_approach/Data/05Biomes/')
if isfile('No_mean_PCA_biomes_seasonal_9_v2_5_perc.mat')
    disp('File already exists!')
else
    save('No_mean_PCA_biomes_seasonal_9_v2_5_perc','season_map_smooth','season_map',...
        'uncertainty_smooth_seasonal_map','uncertainty_season_map','Season_obs','corrected_monthly_smooth')
end


%% Analysis

% =========================================================================
% The analysis section is performed on the stand-alone analysis script
% "Analysis.m" and "AnalysisEnv.m"
% =========================================================================

Analysis

AnalysisEnv