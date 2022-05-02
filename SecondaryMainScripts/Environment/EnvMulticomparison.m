%% Do analysis of differences on monthly scale
cd('/net/kryo/work/ursho/PhD/Projects/Biomes/Scripts/Biomes_PROOCE/Data/07Analysis/Environment/')
load('Simple_sort_Env_Data.mat')


map_env = ones(13,12,180,360).*NaN;


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


n_clusters = 8;
%for twelve parameters and twelve months
mat = ones(12,12,n_clusters+2,n_clusters+2).*NaN;
f = figure;
p = uipanel('Parent',f,'BorderType','none'); 
p.TitlePosition = 'centertop'; 
p.FontSize = 12;
p.FontWeight = 'bold';
p.BackgroundColor = [1 1 1];
for i = 1:size(map_env,1)
    for m = 1:12
        %get all observations of a variable i in all biomes in a month m
        
        [mon_mat_anova] = get_observations_per_biome(map_env(i,m,:,:), new_data_annC(m,:,:),n_clusters); %do this, but for every month, append all months at the end
    
        %perform kruskal wallis test    
        [p,tbl{i},stats] = kruskalwallis(mon_mat_anova,[],'off');

        %multicomparison test after kruskal
        figure()
        [c c_m] = multcompare(stats,'Alpha',0.01,'Display','off');
        
        %get biome label for the values in c
        [r labels] = find(sum(~isnan(mon_mat_anova),1) ~= 0); 
        labels
        %change value in mat if comparison is bellow threshold
        for j = 1:length(c)

            if(c(j,6)<=0.01)  
                mat(i,m,labels(c(j,1))+1,1+labels(c(j,2))) = 0;
                mat(i,m,labels(c(j,2))+1,1+labels(c(j,1))) = 0;
            else
                mat(i,m,labels(c(j,1))+1,1+labels(c(j,2))) = 1;
                mat(i,m,labels(c(j,2))+1,1+labels(c(j,1))) = 1;
            end
        end
        close all
    end
    i
end

%mat is 0 for combinations that are significantly different from each other
%at p <= 1%, 1 for combinations that are not significantly different from
%each other, and NaN for non-existent comparisons
%%        
       
size(mat)
sum_mat = squeeze(sum(mat,2,'omitnan'));
size(sum_mat)
env_labels={ 'N', 'P', 'Si', 'P*',...
    'SSS','SST', 'MLD', 'log(NPP)', 'PAR', 'log(chl)', 'pCO_{2}', 'Wind'}
Get_env_combination

combined_mat = sum_mat;
for i =1:size(combined_mat,2)
    combined_mat(:,i,i) = 13;
end


for k = 2:2:size(combined_mat,1)
    for i = 1:size(combined_mat,2)
        for j = i+1:size(combined_mat,2)
            combined_mat(k-1,j,i) = sum_mat(k,i,j);
        end
    end
end
size(combined_mat)
%%
let_plot1 = {'(a)','(b)','(c)','(d)','(e)','(f)'};
let_plot= let_plot1
f = figure;
p = uipanel('Parent',f,'BorderType','none'); 
p.TitlePosition = 'centertop'; 
p.FontSize = 12;
p.FontWeight = 'bold';
p.BackgroundColor = [1 1 1];
j = 1;
env_labels={ 'N', 'P', 'Si', 'P*',...
    'SSS','SST', 'MLD', 'log(NPP)', 'PAR', 'log(chl)', 'pCO_{2}', 'Wind'}
upper = zeros(size(combined_mat,1),1);
lower = zeros(size(combined_mat,1),1);
for i = 1:2:size(combined_mat,1)
    subplot(3,2,j,'Parent',p)
%figure
    tmp_mat = combined_mat;
    tmp_mat(combined_mat > 3 & combined_mat ~= 13) = 0; %for not significantly differente
    tmp_mat(combined_mat <= 3) = 1; %for significantly different
    tmp_mat(combined_mat == 13) = 3; %diagonal elements
    pcolor(squeeze(tmp_mat(i,2:end,2:end)))
    set(gca,'Ydir','reverse')
    for ii1 = 1:size(combined_mat,2)-2
        for ii2 = 1:size(combined_mat,2)-2
            if(ii1 ~= ii2)
                text(ii2 + 0.2, ii1 + 0.5, num2str(combined_mat(i,ii1+1,ii2+1)),'FontSize',16)
            end
        end
    end

    
    xticks([1.5:1:10.5])
    yticks([1.5:1:10.5])
    xticklabels( {'(1) TRP','(2) HIL','(3) WIS','(4) SUS','(5) HIT ', '(6) MTR',...
    '(7) PEU', '(8) SMN'})
    xtickangle(90)
    yticklabels( {'(1) TRP','(2) HIL','(3) WIS','(4) SUS','(5) HIT ', '(6) MTR',...
    '(7) PEU', '(8) SMN'})
    colormap([[1 0 0];[1 1 1];[0 0 0]]);
    set(gca,'Ydir','reverse')
    descr1 = env_labels(i);
    if(i < 13)
    descr2 = env_labels(i+1);
    else
        descr2 = '';
    end
    
    set(findall(gcf,'-property','LineWidth'),'LineWidth',2)
    
    set(gca, 'XAxisLocation', 'top')
    text(9.5,4,descr1,'FontSize',25)
    h = text(4,10,descr2,'FontSize',25)
    text(-0.5,0,let_plot(j),'FontSize',25)
    j = j+1;
    axis square
    alpha(.5)

end
%%
upper(upper == 0) = [];
lower(lower == 0) = [];
[upper, lower]
all = NaN;
for i = 1:length(upper)
    all = [all;upper(i);lower(i)];
end
all(1) = [];
[B, I] = sort(all,'ascend')
env_labels{I}

%%
% =========================================================================
% Find out on which months the differentiation is not significant
% =========================================================================

size(mat) %second 12 is months
sum_mat = squeeze(sum(mat,2,'omitnan'));
size(sum_mat)
env_labels={ 'N', 'P', 'Si', 'P*',...
    'SSS','SST', 'MLD', 'log(NPP)', 'PAR', 'log(chl)', 'pCO_{2}', 'Wind'}
Get_env_combination

combined_mat = sum_mat;
for i =1:size(combined_mat,2)
    combined_mat(:,i,i) = 13;
end


for k = 2:2:size(combined_mat,1)
    for i = 1:size(combined_mat,2)
        for j = i+1:size(combined_mat,2)
            combined_mat(k-1,j,i) = sum_mat(k,i,j);
        end
    end
end
size(combined_mat)


%%
for mm = 1:12
    mat(isnan(mat)) = 0;
    combined_mat = squeeze(mat(:,mm,:,:));
    
    for i =1:size(combined_mat,2)
        combined_mat(:,i,i) = 13;
    end


    for k = 2:2:size(combined_mat,1)
        for i = 1:size(combined_mat,2)
            for j = i+1:size(combined_mat,2)
                combined_mat(k-1,j,i) = mat(k,mm,i,j);
            end
        end
    end
    let_plot1 = {'(a)','(b)','(c)','(d)','(e)','(f)'};
    let_plot= let_plot1
    f = figure;
    p = uipanel('Parent',f,'BorderType','none'); 
    p.TitlePosition = 'centertop'; 
    p.FontSize = 12;
    p.FontWeight = 'bold';
    p.BackgroundColor = [1 1 1];
    j = 1;
    env_labels={ 'N', 'P', 'Si', 'P*',...
        'SSS','SST', 'MLD', 'log(NPP)', 'PAR', 'log(chl)', 'pCO_{2}', 'Wind'}

    for i = 1:2:size(combined_mat,1)
        subplot(3,2,j,'Parent',p)
        tmp_mat = combined_mat;
        tmp_mat(combined_mat > 3 & combined_mat ~= 13) = 0; %for not significantly differente
        tmp_mat(combined_mat <= 3) = 1; %for significantly different
        tmp_mat(combined_mat == 13) = 3; %diagonal elements
        pcolor(squeeze(tmp_mat(i,2:end,2:end)))
        set(gca,'Ydir','reverse')
        for ii1 = 1:size(combined_mat,2)-2
            for ii2 = 1:size(combined_mat,2)-2
                if(ii1 ~= ii2)
                    text(ii2 + 0.2, ii1 + 0.5, num2str(combined_mat(i,ii1+1,ii2+1)),'FontSize',16)
                end
            end
        end


        xticks([1.5:1:10.5])
        yticks([1.5:1:10.5])
        xticklabels( {'(1) TRP','(2) HIL','(3) WIS','(4) SUS','(5) HIT ', '(6) MTR',...
        '(7) PEU', '(8) SMN'})
        xtickangle(90)
        yticklabels( {'(1) TRP','(2) HIL','(3) WIS','(4) SUS','(5) HIT ', '(6) MTR',...
        '(7) PEU', '(8) SMN'})
        colormap([[1 0 0];[1 1 1];[0 0 0]]);
        set(gca,'Ydir','reverse')
        descr1 = env_labels(i);
        if(i < 13)
        descr2 = env_labels(i+1);
        else
            descr2 = '';
        end

        set(findall(gcf,'-property','LineWidth'),'LineWidth',2)

        set(gca, 'XAxisLocation', 'top')
        text(9.5,4,descr1,'FontSize',25)
        h = text(4,10,descr2,'FontSize',25)
        text(-0.5,0,let_plot(j),'FontSize',25)
        j = j+1;
        axis square
        alpha(.5)
    end
end


