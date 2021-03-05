function [smooth_agg_map, raw_agg_map] = aggregate_months(monthly_data,weights,area_map,Obs,frac,flag,conn)
% Function aggregates the maps given into a single map based on the label
% occurrence frequency across time steps. Afterwards it smooths the
% aggregated maps to ensure that the minimum fraction of global ocean area
% is still met by all biomes on the aggregated time scale


    raw_agg_map = ones(1,180,360).* NaN;
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
                    D_sum = ones(length(r),1).*NaN;
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


