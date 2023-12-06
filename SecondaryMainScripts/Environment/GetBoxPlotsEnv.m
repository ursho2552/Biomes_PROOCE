%% Get maps for each environmental variable

cd('/net/kryo/work/ursho/PhD/Projects/Biomes/Scripts/Biomes_PROOCE/Data/07Analysis/Environment/')
load('MonhtlyEnv_Oct2019.mat')

Surf_env = Complete_env;
labels_surf = { 'N_{surf}', 'P_{surf}', 'Si_{surf}', 'P*_{surf}', 'N*_{surf}',...
    'SSS','SST', 'MLD', 'log(NPP)', 'PAR', 'log(Chl)', 'pCO_{2}', 'Wind','\Delta P'}

%create simple sort environmental data

%load help variables
cd('/net/kryo/work/ursho/PhD/Projects/Biomes/Scripts/Biomes_PROOCE/Data/')
load('HelpVariables.mat')
%construct or load simplified version of raw data

Surf_env_simple = Surf_env;
for i = 1:length(LatLon)
    [r c] = find(Surf_env(:,2) == LatLon(i,2) & Surf_env(:,3) == LatLon(i,3));
    for j =1:length(r)
        Surf_env_simple(r(j),2:3) = LatLon(i,4:5);
    end
    i
end

%save the dataset
cd('/net/kryo/work/ursho/PhD/Projects/Biomes/Scripts/Biomes_PROOCE/Data/07Analysis/Environment/')
if isfile('Simple_sort_Env_Data_Oct2019.mat')
    disp('File already exists!')
else
    save('Simple_sort_Env_Data_Oct2019','Surf_env_simple')
end
load('Simple_sort_Env_Data_Oct2019.mat')

cd(folder_main)



%transform to a (13x12x180x360) map (first 12 is the month

map_env = ones(13,12,180,360).*NaN;


for i = 1:size(Surf_env_simple,2)-4
    map_env(i,:,:,:) = prepare2plot( Surf_env_simple(:,[2:4,i+4]));
end






%% construct boxplots for each env variable throughout the months

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
new_data_annC = corr_corrected_monthly_smooth;
%omit ninth cluster
new_data_annC(new_data_annC == 9) = NaN;

%%
labels_surf = { 'N_{surf}', 'P_{surf}', 'Si_{surf}', 'P*_{surf}', 'N*_{surf}',...
    'SSS','SST', 'MLD', 'log(NPP)', 'PAR', 'log(Chl)', 'pCO_{2}', 'Wind'}
%get rid of N*
map_env(5,:,:,:) = [];


n_clusters = 8;
mat = ones(12,12,n_clusters+2,n_clusters+2);
f = figure;
p = uipanel('Parent',f,'BorderType','none'); 
%p.Title = 'Kruskal-Wallis Test'; 
p.TitlePosition = 'centertop'; 
p.FontSize = 12;
p.FontWeight = 'bold';
p.BackgroundColor = [1 1 1];
boxtitle2 = { 'N', 'P', 'Si', 'P*',...
    'SSS','SST', 'MLD', 'log(NPP)', 'PAR', 'log(chl)', 'pCO_{2}', 'Wind'}
median_env = ones(12,n_clusters).*NaN;

%The sequence was calculated in Analysis.m but it should be as follows:
%sequence = [2 5 4 3 8 1 6 7]
%This also applies to the xlabels which are hardcoded below

IQR = median_env;
for i = 1:size(map_env,1)
    re_anova_annual = ones(12*180*360,n_clusters).*NaN;

    for nn = 1:n_clusters
        %use sequence in latitude to draw box diagram
        n = sequence(nn)
        tmp = squeeze(map_env(i,:,:,:));
        tmp1 = tmp(new_data_annC == n);
        re_anova_annual(1:length(tmp1),nn) = tmp1;
    
    end
    median_env(i,:) = median(re_anova_annual,1,'omitnan');
    IQR(i,:) = prctile(re_anova_annual,75,1)-prctile(re_anova_annual,25,1)
    
    
        subplot(3,4,i,'Parent',p)
        boxplot(re_anova_annual)
        title(boxtitle2{i})
        if(i >= 9)
        xticklabels({'(2) HIL','(5) HIT','(4) SUS','(3) WIS', '(8) SMN','(1) TRP',...
            '(6) MTR','(7) PEU'})

        xtickangle(90)
        end
end


%% get ranking for each variable

for i = 1:size(median_env,1)
     [B I] = sort(median_env(i,:),'descend');
     r_all = NaN;
     for j = 1:8
         [r c] = find(I == j);
         r_all = [r_all,c];
         
     end
     r_all(1) = [];
     r_all
     
     jj = input('next? (Press Enter to continue)')
end

for i = 1:size(IQR,1)
     [B I] = sort(IQR(i,:),'descend');
     r_all = NaN;
     for j = 1:8
         [r c] = find(I == j);
         r_all = [r_all,c];
         
     end
     r_all(1) = [];
     r_all
     
     jj = input('next? (Press Enter to continue)')
end


%%  

cd('/net/kryo/work/ursho/PhD/Projects/Biomes/Scripts/Biomes_PROOCE/Data/07Analysis/Environment/')
load('Simple_sort_Env_Data.mat')



%transform to a (13x12x180x360) map (first 12 is the month

map_env = ones(13,12,180,360).*NaN;

for i = 1:size(Surf_env_simple,2)-4
    map_env(i,:,:,:) = prepare2plot( Surf_env_simple(:,[2:4,i+4]));
end




%% construct boxplots for each env variable throughout the months

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
new_data_annC = corr_corrected_monthly_smooth;
%omit ninth cluster
new_data_annC(new_data_annC == 9) = NaN;


labels_surf = { 'N_{surf}', 'P_{surf}', 'Si_{surf}', 'P*_{surf}',...
    'SSS','SST', 'MLD', 'log(NPP)', 'PAR', 'log(Chl)', 'pCO_{2}', 'Wind'}
%get rid of N*
map_env(5,:,:,:) = [];

not_signif = [NaN NaN NaN];%parameter|month|biome1|biome2
n_clusters = 8;
for i = 1:size(map_env,1)
    re_anova_annual = ones(12*180*360,n_clusters).*NaN;
%     for m = 1:12
        for n = 1:n_clusters

            tmp = squeeze(map_env(i,:,:,:));
            tmp1 = tmp(new_data_annC == n);
            re_anova_annual(1:length(tmp1),n) = tmp1;

        end
    
        [p,tbl{i},stats] = kruskalwallis(re_anova_annual,[],'off');

        %multicomparison test after kruskal
        [c c_m] = multcompare(stats,'Alpha',0.01,'Display','off');
        
        [r cc] = find(c(:,end) >0.01)
        if(~isempty(r))
            not_signif = [not_signif;[r.*0 + i,c(r,1:2)]];
        end
        
        
%     end
end
not_signif(1,:) = []; 
labels_surf(not_signif(:,1))
unique_pairs = unique(not_signif(:,2:3),'rows')
    
 
