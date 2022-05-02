% =========================================================================
% Get niches of modeled species based on the presence absence output and
% our environmental data
% =========================================================================


% each species should have SST, N, P, Si, MLD, PAR, SSS, pCO2

load('/net/kryo/work/ursho/PhD/Projects/Biomes/Scripts/Biomes_PROOCE/Data/07Analysis/Environment/EnvData.mat')

nitrate_surf = permute(squeeze(nitrate(:,:,1,:)),[3,2,1]);
phosphate_surf = permute(squeeze(phosphate(:,:,1,:)),[3,2,1]);
silicate_surf = permute(squeeze(silicate(:,:,1,:)),[3,2,1]);
salinity_surf = permute(squeeze(salinity(:,:,1,:)),[3,2,1]);
temperature_surf = permute(squeeze(temperature(:,:,1,:)),[3,2,1]);
MLD_surf = permute(MLD,[3,2,1]);
PAR_surf = permute(PAR,[3,2,1]);
pCO2_surf = permute(pCO2,[3,2,1]);
wind_surf = permute(wind,[3,2,1]);

% =========================================================================
% Calculate niches for each species
% =========================================================================

%initialize vectors for each species
niches_modelled = ones(3,536,9).*NaN;

%get maps for each species and extract the values of environmental data
for i = 1:536
    spec_map = prepare2plotV2(No_nan_phyto_simple(:,[2 3 4 i+4]));
    %SST
    niches_modelled(1,i,1) = prctile(temperature_surf(spec_map == 1),25);
    niches_modelled(2,i,1) = prctile(temperature_surf(spec_map == 1),50);
    niches_modelled(3,i,1) = prctile(temperature_surf(spec_map == 1),75);
    %N
    niches_modelled(1,i,2) = prctile(nitrate_surf(spec_map == 1),25);
    niches_modelled(2,i,2) = prctile(nitrate_surf(spec_map == 1),50);
    niches_modelled(3,i,2) = prctile(nitrate_surf(spec_map == 1),75);
    %P
    niches_modelled(1,i,3) = prctile(phosphate_surf(spec_map == 1),25);
    niches_modelled(2,i,3) = prctile(phosphate_surf(spec_map == 1),50);
    niches_modelled(3,i,3) = prctile(phosphate_surf(spec_map == 1),75);
    %Si
    niches_modelled(1,i,4) = prctile(silicate_surf(spec_map == 1),25);
    niches_modelled(2,i,4) = prctile(silicate_surf(spec_map == 1),50);
    niches_modelled(3,i,4) = prctile(silicate_surf(spec_map == 1),75);
    %SSS
    niches_modelled(1,i,5) = prctile(salinity_surf(spec_map == 1),25);
    niches_modelled(2,i,5) = prctile(salinity_surf(spec_map == 1),50);
    niches_modelled(3,i,5) = prctile(salinity_surf(spec_map == 1),75);
    %MLD
    niches_modelled(1,i,6) = prctile(MLD_surf(spec_map == 1),25);
    niches_modelled(2,i,6) = prctile(MLD_surf(spec_map == 1),50);
    niches_modelled(3,i,6) = prctile(MLD_surf(spec_map == 1),75);
    %PAR
    niches_modelled(1,i,7) = prctile(PAR_surf(spec_map == 1),25);
    niches_modelled(2,i,7) = prctile(PAR_surf(spec_map == 1),50);
    niches_modelled(3,i,7) = prctile(PAR_surf(spec_map == 1),75);
    %pCO2
    niches_modelled(1,i,8) = prctile(pCO2_surf(spec_map == 1),25);
    niches_modelled(2,i,8) = prctile(pCO2_surf(spec_map == 1),50);
    niches_modelled(3,i,8) = prctile(pCO2_surf(spec_map == 1),75);
    %wind
    niches_modelled(1,i,9) = prctile(wind_surf(spec_map == 1),25);
    niches_modelled(2,i,9) = prctile(wind_surf(spec_map == 1),50);
    niches_modelled(3,i,9) = prctile(wind_surf(spec_map == 1),75);
    i
end

%%
load('Coded_values_SOM_11Oct2019_V3_species.mat', 'network_species')
network_species = unique([final_pairs(:,[1,3]);final_pairs(:,[2,3])],'rows')
[B,I] = sort(network_species(:,2),'ascend')

network_species = network_species(I,:)

my_species = [network_species(:,2),niches_modelled(1,network_species(:,1),1)',niches_modelled(2,network_species(:,1),1)',niches_modelled(3,network_species(:,1),1)',...
    niches_modelled(1,network_species(:,1),2)',niches_modelled(2,network_species(:,1),2)',niches_modelled(3,network_species(:,1),2)',...
    niches_modelled(1,network_species(:,1),3)',niches_modelled(2,network_species(:,1),3)',niches_modelled(3,network_species(:,1),3)',...
    niches_modelled(1,network_species(:,1),4)',niches_modelled(2,network_species(:,1),4)',niches_modelled(3,network_species(:,1),4)',...
    niches_modelled(1,network_species(:,1),5)',niches_modelled(2,network_species(:,1),5)',niches_modelled(3,network_species(:,1),5)',...
    niches_modelled(1,network_species(:,1),6)',niches_modelled(2,network_species(:,1),6)',niches_modelled(3,network_species(:,1),6)',...
    niches_modelled(1,network_species(:,1),7)',niches_modelled(2,network_species(:,1),7)',niches_modelled(3,network_species(:,1),7)',...
    niches_modelled(1,network_species(:,1),8)',niches_modelled(2,network_species(:,1),8)',niches_modelled(3,network_species(:,1),8)',...
    niches_modelled(1,network_species(:,1),9)',niches_modelled(2,network_species(:,1),9)',niches_modelled(3,network_species(:,1),9)'],


% network_species = unique([final_pairs(:,[1,3]);final_pairs(:,[2,3])],'rows')
% [B,I] = sort(network_species(:,2),'ascend')
% 
% network_species = network_species(I,:)
% network_species = unique(network_species(:,1))
% my_species2 = [network_species(:,1),niches_modelled(1,network_species(:,1),1)',niches_modelled(2,network_species(:,1),1)',niches_modelled(3,network_species(:,1),1)',...
%     niches_modelled(1,network_species(:,1),2)',niches_modelled(2,network_species(:,1),2)',niches_modelled(3,network_species(:,1),2)',...
%     niches_modelled(1,network_species(:,1),3)',niches_modelled(2,network_species(:,1),3)',niches_modelled(3,network_species(:,1),3)',...
%     niches_modelled(1,network_species(:,1),4)',niches_modelled(2,network_species(:,1),4)',niches_modelled(3,network_species(:,1),4)',...
%     niches_modelled(1,network_species(:,1),5)',niches_modelled(2,network_species(:,1),5)',niches_modelled(3,network_species(:,1),5)',...
%     niches_modelled(1,network_species(:,1),6)',niches_modelled(2,network_species(:,1),6)',niches_modelled(3,network_species(:,1),6)',...
%     niches_modelled(1,network_species(:,1),7)',niches_modelled(2,network_species(:,1),7)',niches_modelled(3,network_species(:,1),7)',...
%     niches_modelled(1,network_species(:,1),8)',niches_modelled(2,network_species(:,1),8)',niches_modelled(3,network_species(:,1),8)',...
%     niches_modelled(1,network_species(:,1),9)',niches_modelled(2,network_species(:,1),9)',niches_modelled(3,network_species(:,1),9)'],
%%

yvalues =  {'(1) TRP','(2) HIL','(3) WIS','(4) SUS','(5) HIT ', '(6) MTR',...
    '(7) PEU','(8) SMN'}
legend_names = {'(1) TRP','(2) HIL','(3) WIS','(4) SUS','(5) HIT ', '(6) MTR',...
    '(7) PEU','Median'}



f = figure;
p = uipanel('Parent',f,'BorderType','none'); 
%p.Title = 'Kruskal-Wallis Test'; 
p.TitlePosition = 'centertop'; 
p.FontSize = 12;
p.FontWeight = 'bold';
p.BackgroundColor = [1 1 1];


%plot sst for all biomes
data = my_species(:,2:4)
subplot(3,3,1,'Parent',p)
hold on
first_obs = my_species(1,1);
flag = 1
for n = 1:length(data)
    next_obs = my_species(n,1);
    if(flag) %first one
        pos = my_species(n,1);
        flag = 0;
    elseif(~flag & next_obs == first_obs)
        pos = pos + 0.0125;
    elseif(~flag & next_obs ~= first_obs)
        pos = my_species(n,1);
        first_obs = next_obs;
    end
    
    
    xneg = abs(data(n,1)-data(n,2));
    xpos = abs(data(n,3)-data(n,2));
    herr(my_species(n,1)) = errorbar(data(n,2),pos,0,0,xneg,xpos,'Color',comb_cmap(my_species(n,1),:))
    hh(n) = plot(data(n,2),[pos],'k*')
    
end
xlabel('SST')
ylim([0.5 8.5])

set(gca,'ytick',[])
set(gca,'yticklabel',[])


%plot salinity for all biomes
data = my_species(:,5:7)
subplot(3,3,2,'Parent',p)

hold on
first_obs = my_species(1,1);
flag = 1
for n = 1:length(data)
    next_obs = my_species(n,1);
    if(flag) %first one
        pos = my_species(n,1);
        flag = 0;
    elseif(~flag & next_obs == first_obs)
        pos = pos + 0.0125;
    elseif(~flag & next_obs ~= first_obs)
        pos = my_species(n,1);
        first_obs = next_obs;
    end
    
    
    xneg = abs(data(n,1)-data(n,2));
    xpos = abs(data(n,3)-data(n,2));
    errorbar(data(n,2),pos,0,0,xneg,xpos,'Color',comb_cmap(my_species(n,1),:))
    plot(data(n,2),[pos],'k*')
    
end
xlabel('N')
ylim([0.5 8.5])

set(gca,'ytick',[])
set(gca,'yticklabel',[])


%plot Temperature for all biomes
data = my_species(:,8:10)
subplot(3,3,3,'Parent',p)
hold on
first_obs = my_species(1,1);
flag = 1
for n = 1:length(data)
    next_obs = my_species(n,1);
    if(flag) %first one
        pos = my_species(n,1);
        flag = 0;
    elseif(~flag & next_obs == first_obs)
        pos = pos + 0.0125;
    elseif(~flag & next_obs ~= first_obs)
        pos = my_species(n,1);
        first_obs = next_obs;
    end
    
    
    xneg = abs(data(n,1)-data(n,2));
    xpos = abs(data(n,3)-data(n,2));
    errorbar(data(n,2),pos,0,0,xneg,xpos,'Color',comb_cmap(my_species(n,1),:))
    plot(data(n,2),[pos],'k*')
    
end
xlabel('P')
ylim([0.5 8.5])

set(gca,'ytick',[])
set(gca,'yticklabel',[])



%plot MLD for all biomes
data = my_species(:,11:13)
subplot(3,3,4,'Parent',p)

hold on
first_obs = my_species(1,1);
flag = 1
for n = 1:length(data)
    next_obs = my_species(n,1);
    if(flag) %first one
        pos = my_species(n,1);
        flag = 0;
    elseif(~flag & next_obs == first_obs)
        pos = pos + 0.0125;
    elseif(~flag & next_obs ~= first_obs)
        pos = my_species(n,1);
        first_obs = next_obs;
    end
    
    
    xneg = abs(data(n,1)-data(n,2));
    xpos = abs(data(n,3)-data(n,2));
    errorbar(data(n,2),pos,0,0,xneg,xpos,'Color',comb_cmap(my_species(n,1),:))
    plot(data(n,2),[pos],'k*')
    
end
xlabel('Si')
ylim([0.5 8.5])

set(gca,'ytick',[])
set(gca,'yticklabel',[])


%plot Pstar for all biomes
data = my_species(:,14:16)
subplot(3,3,5,'Parent',p)

hold on
first_obs = my_species(1,1);
flag = 1
for n = 1:length(data)
    next_obs = my_species(n,1);
    if(flag) %first one
        pos = my_species(n,1);
        flag = 0;
    elseif(~flag & next_obs == first_obs)
        pos = pos + 0.0125;
    elseif(~flag & next_obs ~= first_obs)
        pos = my_species(n,1);
        first_obs = next_obs;
    end
    
    
    xneg = abs(data(n,1)-data(n,2));
    xpos = abs(data(n,3)-data(n,2));
    errorbar(data(n,2),pos,0,0,xneg,xpos,'Color',comb_cmap(my_species(n,1),:))
    plot(data(n,2),[pos],'k*')
    
end
xlabel('SSS')
ylim([0.5 8.5])

set(gca,'ytick',[])
set(gca,'yticklabel',[])



%plot MLPAR for all biomes
data = my_species(:,17:19)
subplot(3,3,6,'Parent',p)

hold on
first_obs = my_species(1,1);
flag = 1
for n = 1:length(data)
    next_obs = my_species(n,1);
    if(flag) %first one
        pos = my_species(n,1);
        flag = 0;
    elseif(~flag & next_obs == first_obs)
        pos = pos + 0.0125;
    elseif(~flag & next_obs ~= first_obs)
        pos = my_species(n,1);
        first_obs = next_obs;
    end
    
    
    xneg = abs(data(n,1)-data(n,2));
    xpos = abs(data(n,3)-data(n,2));
    errorbar(data(n,2),pos,0,0,xneg,xpos,'Color',comb_cmap(my_species(n,1),:))
    plot(data(n,2),[pos],'k*')
    
end
xlabel('MLD')
ylim([0.5 8.5])

set(gca,'ytick',[])
set(gca,'yticklabel',[])

data = my_species(:,20:22)
subplot(3,3,7,'Parent',p)

hold on
first_obs = my_species(1,1);
flag = 1
for n = 1:length(data)
    next_obs = my_species(n,1);
    if(flag) %first one
        pos = my_species(n,1);
        flag = 0;
    elseif(~flag & next_obs == first_obs)
        pos = pos + 0.0125;
    elseif(~flag & next_obs ~= first_obs)
        pos = my_species(n,1);
        first_obs = next_obs;
    end
    
    
    xneg = abs(data(n,1)-data(n,2));
    xpos = abs(data(n,3)-data(n,2));
    errorbar(data(n,2),pos,0,0,xneg,xpos,'Color',comb_cmap(my_species(n,1),:))
    plot(data(n,2),[pos],'k*')
    
end
xlabel('PAR')
ylim([0.5 8.5])

set(gca,'ytick',[])
set(gca,'yticklabel',[])


data = my_species(:,23:25)
subplot(3,3,8,'Parent',p)

hold on
first_obs = my_species(1,1);
flag = 1
for n = 1:length(data)
    next_obs = my_species(n,1);
    if(flag) %first one
        pos = my_species(n,1);
        flag = 0;
    elseif(~flag & next_obs == first_obs)
        pos = pos + 0.0125;
    elseif(~flag & next_obs ~= first_obs)
        pos = my_species(n,1);
        first_obs = next_obs;
    end
    
    
    xneg = abs(data(n,1)-data(n,2));
    xpos = abs(data(n,3)-data(n,2));
    errorbar(data(n,2),pos,0,0,xneg,xpos,'Color',comb_cmap(my_species(n,1),:))
    plot(data(n,2),[pos],'k*')
    
end
xlabel('pCO_2')
ylim([0.5 8.5])

set(gca,'ytick',[])
set(gca,'yticklabel',[])

data = my_species(:,26:28)
subplot(3,3,9,'Parent',p)

hold on
first_obs = my_species(1,1);
flag = 1
for n = 1:length(data)
    next_obs = my_species(n,1);
    if(flag) %first one
        pos = my_species(n,1);
        flag = 0;
    elseif(~flag & next_obs == first_obs)
        pos = pos + 0.0125;
    elseif(~flag & next_obs ~= first_obs)
        pos = my_species(n,1);
        first_obs = next_obs;
    end
    
    
    xneg = abs(data(n,1)-data(n,2));
    xpos = abs(data(n,3)-data(n,2));
    herr(my_species(n,1)) = errorbar(data(n,2),pos,0,0,xneg,xpos,'Color',comb_cmap(my_species(n,1),:))
    hh(n) = plot(data(n,2),[pos],'k*')
    
end
xlabel('wind')
ylim([0.5 8.5])

set(gca,'ytick',[])
set(gca,'yticklabel',[])
legend([herr(1:7),hh(1)],legend_names),



%% 
% =========================================================================
% Map co-occurrences
% =========================================================================
final_pairs


for i = 1:8
    tmp = final_pairs(final_pairs(:,3) == i,1:2);
    co_map = ones(1,180,360).*0;
    for j= 1:size(tmp,1)
        dat_tmp = No_nan_phyto_simple(:,[tmp(j,:)+4]);
        dat_tmp = sum(dat_tmp,2);
        dat_tmp(dat_tmp < 2) = 0;
        dat_tmp(dat_tmp > 0) = 1;
        map_tmp = prepare2plot([No_nan_phyto_simple(:,[2 3 4]),dat_tmp]);
        map_tmp = sum(map_tmp,1,'omitnan');
        co_map = co_map+map_tmp;
    end
    co_map(isnan(corr_smooth_annual_map)) = NaN;
    plotSOM(co_map,1,NaN)
    title(i)
%     cmap = ametrine;
    colormap(comb_cmap);
    cc = colorbar;
%     caxis([0 150]);
    
end



