function [raw_map] = Clean_up_biomes( raw_map,weights,...
    area_map,fraction,conn,flag)
% Function creates spatially coherent biomes by reassigning small patches 
% to the surrounding most similar patch

%{
Parameters:
    raw_map (matrix): 2D matrix of aggregated maps with the dimension
        1 x 180 x 360 
    weights (matrix): Weights of trained neurons with dimension m x k, with
        m as the number of clusters and k the number of features
    area_map (matrix) : Map of the area per 1°-pixel
    fraction (float): Minimum area needed for a biome. This area is given as a
        percentage of area of global surface ocean where there are observations
    conn (int): Connectivity of 1°-pixels. Can be 4, i.e. only pixels with
        adjacent edges or 8, i.e. pixels with adjacent edges and corners
    flag (bool): Flag to visualize the uncertainty, i.e. all 1°-pixels that
        do not have a unique majority occurrence frequency of a label across
        the n months
 
 Output:
    raw_map (matrix): Aggregated maps after smoothing

%}

    %get a map denoting pixels that have at least one observation
    obs_map = mean(raw_map,1,'omitnan');
    %calculate total available area
    tot_area = sum(area_map(~isnan(obs_map)));
    %for each label in raw_map, get the fraction of area coverage
    sequence_labels = unique(raw_map(~isnan(raw_map)));
    area_per_label = NaN(length(sequence_labels),size(raw_map,1));
    for i = 1:length(sequence_labels)
        for j = 1:size(raw_map,1)
            area_per_label(i,j) = sum(area_map(raw_map(j,:,:) == sequence_labels(i)))./tot_area;
        end
    end
    %sort the mean area over the months, seasons,...
    [~, I] = sort(mean(area_per_label,2,'omitnan'),'descend');
    sequence_labels = sequence_labels(I);

    %get similarity matrix for each label
    matrix_sequence = ones(size(sequence_labels,1)).*NaN;
    for i = 1:size(sequence_labels,1)
        neuron_label = sequence_labels(i);
        D = pdist2(weights(neuron_label,:),weights,'cityblock');
        maxclust = size(weights,1);
        D(1,neuron_label) = NaN;
        D = [(1:length(D))',D'];
        %
        D(~ismember(D(:,1),sequence_labels),2) = NaN;
        %
        weights(~ismember(D(:,1),sequence_labels),:) = NaN;

        [B,I] = sort(D(:,2),'descend');

        Z = linkage(weights,'weighted','cityblock');
        Z(isnan(Z(:,3)),:) = [];
        Z(:,3) = [(maxclust+1):(maxclust+size(Z,1))]';
        [similarity_neurons] = find_most_similar(D, Z,neuron_label,maxclust, neuron_label);
        similarity_neurons(similarity_neurons < 1) = [];
        matrix_sequence(i,:) = similarity_neurons';
    end
    %matrix_sequence contains for each label in the first column, the most
    %similar labels in the next columns
    %% Create patch map for each label
    parfor k = 1:size(raw_map,1)
        for i = 1:size(matrix_sequence,1)

            [patch_map,r] = isolate_patch(raw_map(k,:,:),matrix_sequence(i,1),area_map,tot_area,fraction,conn);
            %r contains the IDs for the patches found for each neuron label,
            %reassing each of them until they are no longer found
            changed_map = raw_map(k,:,:);
            for j = 1:length(r)
                for m = 2:size(matrix_sequence,2)
                    %get number of pixels in the patch
                    sum_patch = sum(~isnan(patch_map(patch_map == r(j))));
                    %change label in the copy
                    changed_map(patch_map == r(j)) = matrix_sequence(i,m);
                    %calculate new patches
                    [tmp_patch_map,~] = isolate_patch(changed_map,matrix_sequence(i,m),area_map,tot_area,fraction,conn);
                    %get new patch ID

                    new_ID = unique(tmp_patch_map(patch_map == r(j)));
                    new_ID(isnan(new_ID)) = [];

                    if(~isempty(new_ID))
                        if(length(new_ID) > 1)
                            tbl = tabulate(tmp_patch_map(patch_map == r(j)));
                            [r_,~] = find(tbl(:,3) == max(tbl(:,3)));
                            new_ID = tbl(r_(1),1);
                        end

                        sum_patch_tmp = sum(~isnan(tmp_patch_map(tmp_patch_map == new_ID)));
                        if(sum_patch_tmp > sum_patch)
                            break;
                        else
                            changed_map(patch_map == r(j)) = matrix_sequence(i,1);
                        end
                    else
                        %patch disappeared
                        break;
                    end
                end
                %if the patch is the same as in the beginning delete it
                if(unique(changed_map(patch_map == r(j))) == matrix_sequence(i,1))
                    changed_map(patch_map == r(j)) = NaN;
                end
                %for visualizing uncertainty
                if(flag)
                    changed_map(patch_map == r(j)) = NaN;
                end
            end

                raw_map(k,:,:) = changed_map;


        end


    end
end