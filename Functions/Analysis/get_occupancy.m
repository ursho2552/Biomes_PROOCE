function [coverage_month, coverage] = get_occupancy( area_map, Obs_simple,n_clusters,monthly_maps)
% Function calculated the area coverage of each species on each month. 

%{
Parameters:
    area_map (matrix): Matrix with the area of each 1°-pixel with dimension
        1 x 180 x 360
    Obs_simple (matrix): Vectorized observations
    n_clusters (int): Number of biomes/clusters in the monthly maps
    monthly_maps (matrix): Matrix containing the monthly/seasonal/annual
        biomes. Thie dimensions should be m x 180 x 360. m = 12 for monthly
        biomes

 
 Output:
    coverage_month (matrix): Coverage of each species in each month/time
        step
    coverage (matrix): Coverage of each species averaged across all
        months/time step

%}

    %Calculate Occurrence in terms of covered area
    %use monthly boundaries of biomes
    %for each species create monthly maps of occurrence
    
    
    %for each species (column) in Obs_simple create a map of the presence,
    %and multiply it by the area_map (element wise)
    mean_species_map = NaN(size(Obs_simple,2)-5,size(monthly_maps,1),180,360);
    tic
    for m = 1:size(monthly_maps,1)
        for i = 5:size(Obs_simple,2)-1

            tmp_data = Obs_simple(Obs_simple(:,4) == m,[2 3 4 i]);
            tmp_map = prepare2plot(tmp_data);

            tmp_map = tmp_map.* area_map;

            mean_species_map(i-4,m,:,:) = tmp_map;
             i
        end
    end
    toc



    %from all maps, calculate the coverage for each species in each month
    tic
    coverage = NaN(size(monthly_maps,1),n_clusters,size(Obs_simple,2)-5);
    for m = 1:size(monthly_maps,1)
        for i =1:n_clusters
            area_biome = sum(area_map(monthly_maps(m,:,:) == i));
            for j = 1:size(mean_species_map,1)
                

                coverage(m,i,j) = sum(mean_species_map(j,m,monthly_maps(m,:,:) ==i),'omitnan')/area_biome;
            end
        end
    end
    toc
    coverage(coverage == 0) = NaN;
    coverage_month = coverage;
% =========================================================================
% Only consider the months where the biome is present
% =========================================================================
    coverage = squeeze(mean(coverage,1,'omitnan'));
    coverage(isnan(coverage)) = 0;
end