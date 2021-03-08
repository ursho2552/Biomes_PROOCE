function [smooth_agg_map, raw_agg_map] = aggregate_months(monthly_data,weights,area_map,Obs,frac,flag,conn)
% Function to aggregate the maps given into a single map based on the label
% occurrence frequency across time steps. Afterwards it smooths the
% aggregated maps to ensure that the minimum fraction of global ocean area
% is still met by all biomes on the aggregated time scale

%{
Parameters:
    monthly_data (matrix): 3D matrix of n monthly maps with the dimension
        n x 180 x 360 
    weights (matrix): Weights of trained neurons with dimension m x k, with
        m as the number of clusters and k the number of features
    area_map (matrix) : Map of the area per 1°-pixel
    Obs (matrix): Vectorized observations of the months in monthly_data,
        with the structure |ID|LON|LAT|MONTH|SPECIES1|SPECIES2|...|SPECIESK|BUFFER 
    frac (float): Minimum area needed for a biome. This area is given as a
        percentage of area of global surface ocean where there are observations
    flag (bool): Flag to visualize the uncertainty, i.e. all 1°-pixels that
        do not have a unique majority occurrence frequency of a label across
        the n months
    conn (int): Connectivity of 1°-pixels. Can be 4, i.e. only pixels with
        adjacent edges or 8, i.e. pixels with adjacent edges and corners
 
 Output:
    smooth_agg_map (matrix): Aggregated maps after smoothing biomes to have
        the minimum area requirement
    raw_agg_map (matrix): Aggregated maps without smoothing

%}


    raw_agg_map = NaN(1,180,360);
    for i = 1:180
        for j = 1:360
            tmp = monthly_data(:,i,j);

            tmp(isnan(tmp)) = [];
            if(~isempty(tmp))
                tbl = tabulate(tmp);
                [r, ~] = find(tbl(:,3) == max(tbl(:,3)));
                if(length(r) == 1 && tbl(r,3) > 50) 
                    raw_agg_map(1,i,j) = tbl(r(1),1);
                elseif(length(r) == 1 && tbl(r,3) <= 50 && flag == 0)
                    raw_agg_map(1,i,j) = tbl(r(1),1);
                elseif((length(r) > 1 && flag == 0))


                    obs_all = Obs(Obs(:,2) == j & Obs(:,3) == i, 5:end-1);
                    D_sum = NaN(length(r),1);
                    for k = 1:length(r)
                        D = pdist2(weights(r(k),:),obs_all(tmp ~= r(k),:));
                        D_sum(k) = sum(D,'omitnan');
                    end
                    [rr, ~] = find(D_sum == min(D_sum));

                        raw_agg_map(1,i,j) = tbl(r(rr(1)),1);

                end


            end
        end
    end


    [smooth_agg_map] = Clean_up_biomes( raw_agg_map,weights,...
        area_map,frac,conn,flag);





end


