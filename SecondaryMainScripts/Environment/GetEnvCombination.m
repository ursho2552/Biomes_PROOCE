%% Find optimal environmental predictors

%assume three dimensions are sufficient (?)
%calculate all posible combinations of the 12 variables
%get the overlap for each combination, i.e. how many months are not
%significantly different (above 3 = not sign. diff.)

%calculate how many comparisons are not sign. diff., and how many
%overlapping non significantly different comparisons are between the three
%parameters. Overlap is defined here as two or more of the parameters not
%beeing significantly different for the same biome comparison

%the winning combination of environmental parameters should have the least
%amount of not significantly different comparisons, and the least amount of
%newly created non significantly different comparisons. If there are many
%results, then take the combination with the most overlap

%%
load('/net/kryo/work/ursho/Damiano_Presence_data/presence_absence_tables_ensemble_averages/Group_specific_background_approach/Data/07Analysis/Environment/')
load('MonthlyEnv.mat')
Standardized_env = Complete_env;
Surf_env = Standardized_env(:,[1 2 3 4 6 8 10 12 14 16 18 19 20 21 22 23 24]);
labels_surf = { 'N_{surf}', 'P_{surf}', 'Si_{surf}', 'P*_{surf}', 'N*_{surf}',...
    'SSS','SST', 'MLD', 'log(NPP)', 'PAR', 'log(Chl)', 'pCO_{2}', 'Wind'}
%use mat and take into account in which months two biomes cannot be
%differentiated between each other

tmp_mat = mat;
tmp_mat(:,:,1,:) = [];
tmp_mat(:,:,:,1) = [];

tmp_mat(:,:,:,end) = [];
tmp_mat(:,:,end,:) = [];


for i = 1:size(tmp_mat,3)
    for j = 1:i
        if(j<=i)
            tmp_mat(:,:,i,j) = 0;%NaN;
        end
    end
end
tmp_mat(isnan(tmp_mat)) = 0;
squeeze(tmp_mat(1,1,:,:))

%get matrix with possible combinations of environmental parameters 3D
v = 1:12; %12 env parameters
C = nchoosek(v,2);
C_boxes = ones(size(C,1),3)*NaN; %first  rows for combinations, 
%then number of overlap (number of boxes over 3 and the new boxes over 3)
C_boxes(:,1:2) = C;
pairings = ones(length(C_boxes),8)*NaN;
pairings_over = pairings;
predictors = [NaN, NaN];
for jj = 1:size(tmp_mat,3) %loop over each biome
    available = 1:size(tmp_mat,3);
    for i = 1:length(C)
        %take the first combination of parameters for the 12 months
        tmp_av = available;
        tmp_av(jj) = [];
        a1 = squeeze(tmp_mat(C(i,1),:,:,:));
        a2 = squeeze(tmp_mat(C(i,2),:,:,:));
%         a3 = squeeze(tmp_mat(C(i,3),:,:,:));
        %delete all other comparison, and only retain the comparisons for
        %biome jj
        a1(:,tmp_av,tmp_av) = 0;%NaN;
        a2(:,tmp_av,tmp_av) = 0;%NaN;
%         a3(:,tmp_av,tmp_av) = 0;%NaN;
        %add the three parameters
        a4 = a1 + a2;% + a3;
        %see if parameters overlap, if not signigicantly different in the same
        %month for more than one parameter, change to one, since it is
        %still not significantly different
        a4(a4 > 1) = 1;

        %get the sum
        sum_a4 = squeeze(sum(a4,1,'omitnan'));
         a4 = a1 + a2;% + a3;
        sum_a4_over = squeeze(sum(a4,1,'omitnan'));
        
        threshold_mat = sum_a4;
        %chenge values of thershold_mat to 1 if above 3 and 0 otherwise
        %threshold_mat(sum_a4 < 4) = 0;
        %threshold_mat(sum_a4 > 3) = 1;
        C_boxes(i,3) = sum(sum(threshold_mat,'omitnan'),'omitnan');
        pairings(i,jj) = sum(sum(threshold_mat,'omitnan'),'omitnan');
        pairings_over(i,jj) = sum(sum(sum_a4_over,'omitnan'),'omitnan') -  sum(sum(sum_a4,'omitnan'),'omitnan');
    end

end





%%
tmp_choice = [Surf_env,ones(length(Surf_env),1)*NaN];
%tmp_choice = [MLD_env,ones(length(MLD_env),1)*NaN];
tmp_choice(:,5+4) = [];%delete N*
boxtitle2 = { 'N', 'P', 'Si', 'P*',...
    'SSS','SST', 'MLD', 'log(NPP)', 'PAR', 'log(chl)', 'pCO_{2}', 'Wind'} 
%for each biome (col) in pairings find the minimum
predictors = [NaN,NaN, NaN, NaN];
for i = 1:size(pairings,2)
    %get the combination of parameters which is best able to differentiate
    %between biomes
    [r c] = find(pairings(:,i) == min(pairings(:,i)));
    dissim = NaN;
    if(length(r) ~= 1)
        for j = 1:size(C_boxes(r,1:2),1)
            [ R_spearman_annual_env, orig_annual_env,T_annual_env ] =...
                get_dendrogram('spearman',tmp_choice(:,4+C_boxes(r(j),1:2))',...
                boxtitle2(C_boxes(r(j),1:2)),2:3,'average');
            dissim = [dissim;R_spearman_annual_env(1,3)];
        end
        i
        
         
        close all
        dissim(1) = [];
        [r1 c1] = find(dissim == max(dissim));
        r = r(r1)
    end
    
    predictors = [predictors; [C_boxes(r,1).*0 + i,C_boxes(r,1:3)]];
end
predictors(1,:) = [];
env_labels = boxtitle2
env_labels(predictors(:,2:3))
predictors(:,end)


%%
predictors = [NaN,NaN, NaN NaN];
for i = 1:size(pairings,2)
    [r c] = find(pairings(:,i) == min(pairings(:,i)));
    
    
    predictors = [predictors; [C_boxes(r,1).*0 + i,C_boxes(r,1:2),pairings(r,i)]];
end
predictors(1,:) = [];
env_labels(predictors(:,2:end))

un_pred = unique([predictors(:,2);predictors(:,3)])
tmp = predictors(:,2:3);
for i = 1:length(un_pred)
    r = find(tmp == un_pred(i));
    env_labels(un_pred(i))
    length(r)
end
env_labels(predictors(:,2:3))




