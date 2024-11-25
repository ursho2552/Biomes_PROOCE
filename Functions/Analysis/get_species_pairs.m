function [biome_pairs] = get_species_pairs(Scores_data,coverage_data)
% The function takes in the interaction strength between species pairs
% (adapted from Dunning 1993), and the area coverage of each species. The
% data is for one month.
% The function calculates and identifies the so called "indicator" species
% pairs in a given month, i.e. the species pairs that co-occur in 


%%
    % =========================================================================
    % Initialize matrices containing the scores of the pairs, the minimum area
    % and the ID of the pairs
    % =========================================================================
    PosScores = Scores_data;
    n_features = size(coverage_data,2);
    mat_pairs = NaN(size(PosScores,1),(n_features*(n_features-1)/2));
    mat_pairs_area = NaN(size(PosScores,1),((n_features*(n_features-1))/2));
    ID_pairs = NaN(2,((n_features*(n_features-1))/2));




    for n = 1:size(PosScores,1)
        cc = 1;
        for i = 1:size(PosScores,2)

            for j = i+1:size(PosScores,3)
                %get the value for each pairing
                mat_pairs(n,cc) = PosScores(n,i,j);
                if(~isnan(coverage_data(n,i)) && ~isnan(coverage_data(n,j)))  
                    % =========================================================
                    % Get the smallest of the two monthly coverages, since
                    % that's the minimum two species can co-occur; UHE 11 Oct
                    % 2019
                    % =========================================================
                    mat_pairs_area(n,cc) = min(coverage_data(n,i),coverage_data(n,j));
                end
                %get the ID of the species involved in the pairings
                ID_pairs(:,cc) = [i;j]; 
                cc = cc + 1;
            end
        end
    end
    %%
    % =========================================================================
    % Get the area threshold based on the monthly area coverage of individual
    % species; UHE 11 Oct 2019
    % =========================================================================
    upper_area = prctile(coverage_data,90,2);
    lower_area = prctile(coverage_data,10,2);

    % =========================================================================
    % Flag species pairs that are core (i.e. they have a positive interaction
    % and occur often > 90th percentile). Flag species pairs that are satellite
    % (i.e. they have a negative interaction or they do not occur often <= 10th
    % percentile)
    % =========================================================================

    core_pairs = mat_pairs.*0;
    satellite_pairs = mat_pairs.*0;
    for i = 1:size(upper_area,1)
        core_pairs(i,mat_pairs(i,:) > 0 & mat_pairs_area(i,:) > upper_area(i)) = 1;
        satellite_pairs(i,mat_pairs(i,:) < 0 | mat_pairs_area(i,:) <= lower_area(i)) = 1;
    end

    % =========================================================================
    % Get the sum of the matrices above to find "indicator" pairs
    % =========================================================================
    sum_cores = sum(core_pairs,1,'omitnan');   
    sum_sat = sum(satellite_pairs,1,'omitnan');

    % =========================================================================
    % An "indicator" pair must be a satellite pair in all other biomes. The
    % number of available biomes can be seen from "upper_area", as it will only
    % be NaN if and only if a biome is missing in one month
    % =========================================================================
    i = sum(~isnan(upper_area));

    [~, c_ind] = find(sum_cores == 1 & sum_sat >=i-1);

    biome_pairs = [NaN NaN NaN NaN];
    for ii = 1:length(c_ind)
        [r, ~] = find(core_pairs(:,c_ind(ii)) == 1);
        %check if species are present in most of the biome

        biome_pairs = [biome_pairs;[r, ID_pairs(1,c_ind(ii)), ID_pairs(2,c_ind(ii)), mat_pairs(r,c_ind(ii)) ]];

    end

    biome_pairs(1,:) = [];




end





