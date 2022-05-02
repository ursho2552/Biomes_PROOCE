cd('/net/kryo/work/ursho/PhD/Projects/Biomes/Scripts/Biomes_PROOCE/Data/07Analysis/Environment/')
load('Simple_sort_Env_Data.mat')

%transform to a (13x12x180x360) map (first 13 is the variable

map_env = ones(14,12,180,360).*NaN;


% =========================================================================
% get values of each variable for each neuron
% =========================================================================
%match values in Surf_env_simple to No_nan_phyto_simple
%save the dataset
cd('/net/kryo/work/ursho/PhD/Projects/Biomes/Scripts/Biomes_PROOCE/Data/07Analysis/Environment/')
if isfile('Surf_env_simple_with_ID_Oct2019.mat')
    load('Surf_env_simple_with_ID_Oct2019.mat')
else
    for i = 1:length(Surf_env_simple)
        [r c] = find(No_nan_phyto_simple(:,2) == Surf_env_simple(i,2) &...
            No_nan_phyto_simple(:,3) == Surf_env_simple(i,3) &...
            No_nan_phyto_simple(:,4) == Surf_env_simple(i,4));
        if(~isempty(r))
            Surf_env_simple(i,1) = No_nan_phyto_simple(r,1);
        else
            Surf_env_simple(i,1) = NaN;
        end

    end
    save('Surf_env_simple_with_ID_Oct2019','Surf_env_simple')
end
sum(isnan(Surf_env_simple(:,1)))


%create "new weights" using the vector classes for the neurons, and the
%environmental values using a mean
Surf_env_simple(isnan(Surf_env_simple(:,1)),:) = [];

for i = 5:size(Surf_env_simple,2)
    Surf_env_simple(:,i) = rescale(Surf_env_simple(:,i),-1,1);
end


new_env_weights = NaN(961,13);
Surf_with_classes = [Surf_env_simple, classes(Surf_env_simple(:,1))];
for i = 1:961
    [r c] = find(Surf_with_classes(:,end) == i);
    
    if(~isempty(r))
        new_env_weights(i,:) = mean(Surf_with_classes(r,5:end-1),1,'omitnan');
    end
end
%SST, pco2, wind, Si, chl
data = new_env_weights;
data(:,5) = [];
%% NMDS part 1
% data(isnan(data(:,1)),:) =[];
dissimilarities = pdist(data,'cityblock');

[Y_dis,stress,disparities]  = mdscale(dissimilarities,2,'Start','random');
   
[neuron_maps] = prepare2plot( Surf_with_classes(:,[2:4,end]));
tmp_corrected_neuron = neuron_maps*NaN;
for i =1:12
    if(i<7)
        j = mod(i+6,13);
    else
        j = mod(i+6,13) +1;
    end
    %i and j are the indices of the months that need to be combined
    tmp_corrected_neuron(i,91:end,:) = neuron_maps(i,91:end,:);%raw_monthly_maps(i,91:end,:);%smooth_map(i,91:end,:);
    tmp_corrected_neuron(i,1:90,:) = neuron_maps(j,1:90,:); % raw_monthly_maps(j,1:90,:); %smooth_map(j,1:90,:); % 
end   
corrected_neuron = tmp_corrected_neuron;
monthly_maps = corr_corrected_monthly_smooth;

%% NMDS part 2
figure
hold on;
location_biome = [NaN NaN NaN];
available_labels = unique(monthly_maps(~isnan(monthly_maps)))
%Don't do cluster 9 as we do not analyze it!
available_labels(available_labels == 9) = [];


for ii = 1:length(available_labels)
    %get the first label
    i = available_labels(ii);
    %get all the neurons found in one biome and store as tmp, do it
    %individually for every month, is slower, but needed to match
    %location_tmp (I think)
    tmp = NaN;
    for mm = 1:12
        tmp = [tmp;corrected_neuron(mm,monthly_maps(mm,:,:) == i)'];
    end
    tmp(1) = [];
    tmp(isnan(tmp)) = [];

    
    %get all the 2D values of all neurons transformed, and consider the
    %frequency of occurrence of each neuron (i.e. not only unique neurons)
    Y = Y_dis(tmp,:);

    
    %separate into components x and y  (first dimension and second
    %dimension)
    x = Y(:,1);
    y = Y(:,2);
    XY = [x,y];

    %get some thresholds to define envelope
    p95 = prctile(XY,90,1); 
    p5 = prctile(XY,10,1); 

    XY(XY(:,1) > p95(1) | XY(:,2) > p95(2) | XY(:,1) < p5(1) | XY(:,2) < p5(2),:) = [];

    if(size(unique(XY,'rows'),1) > 2)
        k = convhull(XY(:,1),XY(:,2));
        h0(i) = plot(XY(k,1),XY(k,2),'LineWidth',4,'Color',cmap2(i,:));
    end 
    h1(i) = plot(x,y,'+','LineWidth',1,'Color',cmap2(i,:));
    h2(i) = plot(median(Y(:,1)),median(Y(:,2)),...
                'k*','LineWidth',8); 
            

end
xlabel('Dimension 1')
ylabel('Dimension 2')

%
legend_names = {'(1) TRP','(2) HIL','(3) WIS','(4) SUS','(5) HIT ', '(6) MTR',...
    '(7) PEU','(8) SMN','Cluster Mean'}
tmp = legend_names([available_labels;9])
h0
legend([h1(available_labels), h2(1)],tmp);




%% Alternative without NMDS just using SST and pCO2
% load('/net/kryo/work/ursho/Damiano_Presence_data/presence_absence_tables_ensemble_averages/Group_specific_background_approach/Data/15Environment/Simple_sort_Env_Data.mat')
load('/net/kryo/work/ursho/PhD/Projects/Biomes/Scripts/Biomes_PROOCE/Data/15Environment/Simple_sort_Env_Data_Oct2019.mat')

%transform to a (13x12x180x360) map (first 13 is the variable

map_env = ones(14,12,180,360).*NaN;

for i = 1:size(Surf_env_simple,2)-4
    map_env(i,:,:,:) = prepare2plot( Surf_env_simple(:,[2:4,i+4]));
end

corr_map_env = map_env;
%shift SH by six months
for m = 1:12
    if(m < 7)
        tmp_m = mod(m+6,13);
    else
        tmp_m = mod(m+6,13) + 1;
    end

    
    corr_map_env(:,m,91:end,:) = map_env(:,m,91:end,:);
    corr_map_env(:,m,1:90,:) = map_env(:,tmp_m,1:90,:);
end


map_env = corr_map_env;
map_env(5,:,:,:) = [];
new_data_annC = corr_corrected_monthly_smooth;



clear biome_env




n_clusters = 8
for i =1:n_clusters
    %get data from env_data for a biome
    tmp = map_env; %param x month x lat x lon

    tmp(:,corr_corrected_monthly_smooth ~= i) = NaN;
    tmp(~isnan(tmp));
    biome_env(i,:,:,:,:) = tmp;
end



orig_biome_env = biome_env;

boxtitle2  = { 'N_{surf}', 'P_{surf}', 'Si_{surf}', 'P*_{surf}',...
    'SSS','SST', 'MLD', 'log(NPP)', 'PAR', 'log(Chl)', 'pCO_{2}', 'Wind','\Delta P'}

%%
C = [6,11];% 6 10;10 11]

for mmm = 1:12 
f = figure;
p = uipanel('Parent',f,'BorderType','none'); 
%p.Title = 'Kruskal-Wallis Test'; 
p.TitlePosition = 'centertop'; 
p.FontSize = 12;
p.FontWeight = 'bold';
p.BackgroundColor = [1 1 1];
 
 
clear h1 h2 h3
for j = 1:size(C,1)
    subplot(size(C,1),1,j,'Parent',p)
    hold on;
    cmap2=colormap(parula(9));
 [cmap2] = shuffle_colormap(cmap2);
 [cmap2] = shuffle_colormap(cmap2)
    %figure
    %hold on;
    for i =1:n_clusters
        biome_env = orig_biome_env;
        %if(i ~= 5)
        x = squeeze(biome_env(i,C(j,1),mmm,:,:));
        x = x(~isnan(x));
        y = squeeze(biome_env(i,C(j,2),mmm,:,:));
        y = y(~isnan(y));
        length(y) 
        length(x)
        %get only pairwise complete data
       
        %get only data within xx% confidence interval
        if(length(y) > 2)
            XY = [x,y];
            p95 = prctile(XY,80,1); 
            p5 = prctile(XY,20,1); 
            XY(XY(:,1) > p95(1) | XY(:,2) > p95(2) | XY(:,1) < p5(1) | XY(:,2) < p5(2),:) = [];
            

        h2(i) = plot(median(XY(:,1)),median(XY(:,2)),'k*','LineWidth',8); 
        if(size(unique(XY,'rows'),1) > 2)
            k = convhull(XY(:,1),XY(:,2));
            h0(i) = plot(XY(k,1),XY(k,2),'LineWidth',4,'Color',cmap2(i,:));
        end 

        end


    end

xlabel(boxtitle2(C(j,1)))
ylabel(boxtitle2(C(j,2)))
set(gca, 'XDir','reverse')
box on
grid on
hold off
name_fig = horzcat('Orig_',boxtitle2{C(j,2)},'_',boxtitle2{C(j,1)},'.png')

end
end
legend([h0 h2(1)],{'(1) TRP','(2) HIL','(3) WIS','(4) SUS','(5) HIT ', '(6) MTR',...
    '(7) PEU','(8) SMN','Cluster Mean'})
