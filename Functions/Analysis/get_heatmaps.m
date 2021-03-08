function [top_occ,I] = get_heatmaps(coverage, biome_names, n_clusters,labels_phyto)
% Function Produces heatmaps using the coverage data for the entire dataset, or divided by phylum. 

%{
Parameters:
    coverage (matrix): Matrix containing the mean area coverage of each species
        within each biome
    biome_names (vector(str)): Vector containing names of biomes
    n_clusters (int): Number of biomes/clusters in the monthly maps
    labels_phyto (matrix(str)): Matrix with taxonomic information of species

 
 Output:
    top_occ (matrix): Coverage of species sorted in descending order
    I (vector): IDs of species sorted in descending order of mean coverage

%}

    overall_occ = mean(coverage,1,'omitnan');%global_mean;%

    %sort the coverage
    [B, I] = sort(overall_occ,'descend');

    yvalues = ['global',biome_names];
    top_occ = coverage(:,I);

    all_data = [B;top_occ];
    figure
    hh = imagesc(all_data);
    cmap = ametrine;
    set(gca,'TickLength',[0 0])
    colormap(cmap)
    cc = colorbar;
    yticklabels(yvalues)
    xticklabels({})
    xlabel('Phytoplankton species')
    cc.YTick = 0:0.1:1;
    cc.TickLabels = {'0.0','0.1','0.2','0.3','0.4', '0.5','0.6','0.7','0.8','0.9','1.0'};
    ylabel(cc, 'Average normalized area coverage [-]')
    set(findall(gcf,'-property','FontSize'),'FontSize',30)

    R_species = corr([B;top_occ]','Type','Spearman'); %analyse the similarity


     [ R_mats, orig,T ] = get_dendrogram('spearman',...
         [B;top_occ], yvalues, n_clusters, 'weighted');

    for i = 1:size(R_species,1)
        R_species(i,i) = NaN;
    end
    mean(1-R_species,2,'omitnan')

    [mean(top_occ,2),std(top_occ,[],2)]


    jj = input('get core species across biomes? (Press Enter to continue)')

    %on genus level
    ann_occ = coverage;

    yvalues = ['global',biome_names];

    %%
    % =========================================================================
    % Core, subordinate and satellite species
    % =========================================================================

    CSH = (top_occ.*0)+0.5;

    CSH(top_occ > 0.9) = 1;
    CSH(top_occ <= 0.1) = 0;
    figure
    hh = heatmap(I,biome_names,CSH);
    hh.Colormap = winter(3);
    hh.FontSize = 20;

    %% Heat map per phylum
    %get phylum names

    phylum = unique(labels_phyto(:,3));
    sort_labels = labels_phyto(I,:);
    for i = 1:length(phylum)
        [r, ~] = find(sort_labels(:,3) == phylum{i});
        %get the rows of each phylum
        phylum_occ = top_occ(:,r);

        CSH = (phylum_occ.*0)+0.1;

        CSH(phylum_occ > 0.9) = 0.2;
        CSH(phylum_occ <= 0.1) = 0;

        %only keep those that are core in some biome
        [~, cc] = find(CSH == 0.2);
        if(~ isempty(cc))
             [ R_mats, orig,T] = get_dendrogram('cityblock',...
         phylum_occ, biome_names, n_clusters,'weighted')
            cc = unique(cc);
            length(cc)
            %get heatmap
            figure()
            xvalues = strtrim(string(char(sort_labels(r(cc),1))));

            hh = heatmap(xvalues,biome_names,phylum_occ(:,cc));
            hh.Colormap = winter(11);
        else
            figure()
            xvalues = strtrim(string(char(sort_labels(r,1))));

            hh = heatmap(xvalues,biome_names,phylum_occ);
            hh.Colormap = winter(11);
        end
        disp('The phylum is:')
        phylum{i}

        jj = input('next? (Press Enter to continue)')
    end

    %% Heat maps on genus level

    % get heat map on genus level
    clear genus
    lab_tmp = labels_phyto(:,1);
    for i = 1:length(lab_tmp)
        C = strsplit(lab_tmp{i});
        genus(i) = C(1);
    end

    unique_genus = unique(genus);
    genus_occ = ones(size(ann_occ,1),length(unique_genus))*NaN;
    genus_str = strtrim(string(char(genus)));
    for i = 1:length(unique_genus)
        [r, ~] = find(genus_str == unique_genus{i});
        genus_occ(:,i) = mean(ann_occ(:,r),2);
    end

    overal_genus = mean(genus_occ,1,'omitnan');
    std_overal_genus = std(genus_occ,1,'omitnan');
    [~, I2] = sort(overal_genus,'descend');

    max_overal_genus = overal_genus + std_overal_genus;
    min_overal_genus = overal_genus - std_overal_genus;

    minmax_overal_genus = [max_overal_genus;min_overal_genus];
    minmax_overal_genus = minmax_overal_genus(:,I2);
    %genus_occ([5 9],:) = [];
    genus_occ_sorted = genus_occ(:,I2);

    figure()
    xvalues = strtrim(string(char(unique_genus)));
    xvalues = xvalues(I2);

    hh = heatmap(xvalues,yvalues,[overal_genus(1,I2);genus_occ_sorted]);
    hh.Colormap = winter;

    R_genus = corrcoef(genus_occ_sorted'); %analyse the similarity
    R_genus_glob = corrcoef([genus_occ_sorted;overal_genus(1,I2)]');
    mean(genus_occ_sorted,2)
    std(genus_occ_sorted,[],2)
    
end