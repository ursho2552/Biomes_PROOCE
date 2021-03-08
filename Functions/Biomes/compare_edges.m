function [new_map,tmp_map] = compare_edges(map,nan_map,h,conn)
% Function compares the edges of the 2D map, to ensure that biomes 
% expanding across the edges are recognized correctly

%{
Parameters:
    map (matrix): 2D matrix of with map of biomes/patches (1 x 180 x 360) 
    nan_map (matrix): 2D matrix of with map of missing values (1 x 180 x 360)
    h (int) : Deprecated. Use only h=1
    conn (int): Connectivity of 1°-pixels. Can be 4, i.e. only pixels with
        adjacent edges or 8, i.e. pixels with adjacent edges and corners
 
 Output:
    new_map (matrix): Map of biomes/patches after changing the labels, if a
        biome/patch expands across the edges of the matrix
    tmp_map (matrix): Deprecated. Contains copy of map without any changes

%}

    if(h == 4 || h == 5)

        %% New Version

        new_map = squeeze(map);
        tmp_map = squeeze(nan_map);
        %get left side
        left_edge = new_map(:,1);
        tmp_left = tmp_map(:,1);
        %append left side to the right
        new_map = [new_map,left_edge];
        %extend map by one row above and below (for loop later)
        upper = new_map(1,:)*NaN;
        new_map = [upper;new_map;upper];
        tmp_map = [upper;[tmp_map,tmp_left];upper];

        %start from the right side
        for i = size(new_map,2)-1:-1:1
            %take the left column
            left = new_map(:,i+1);
            %loop over each element in the left column
            for j = 2:length(left)-1
                tmp = left(j);
                %check if pixel is not zero
                if(tmp ~= 0)
                    %if the pixel is not a missing value
                    %if(~isnan(tmp))
                        for k = j-1:j+1
                            if(new_map(k,i) ~= 0 )
                                new_map(new_map == new_map(k,i)) = tmp;

                            elseif(new_map(k,i) == 0 && tmp_map(k,i) == 1)
                                new_map(k,i) = tmp;
                            end
                        end

                end
            end
        end

        new_map(tmp_map == 1) = NaN; 

        new_map(:,end) = [];
        new_map(1,:) = [];
        new_map(end,:) = [];



        tmp_map(:,end) = [];
        tmp_map(1,:) = [];
        tmp_map(end,:) = [];
    else
        %% Changed on 17/07/2019
        %if not already done
        new_map = squeeze(map);
        tmp_map = new_map;

        %get left side
        left_edge = new_map(:,1);
        %get right side
        right_edge = new_map(:,end);

        both_edges = [[NaN;right_edge;NaN],[NaN;left_edge;NaN]];
        pairings = [NaN, NaN];


        for i = 2:length(both_edges)-1
            if(both_edges(i,2) > 0 && ~isnan(both_edges(i,2)))
            %take the three numbers next to the ith pixel
                if(conn == 8)
                    right_tmp = both_edges(i-1:i+1,1);
                elseif(conn == 4)
                    right_tmp = both_edges(i,1);
                end
                %delete NaNs and zeros
                right_tmp(right_tmp == 0 | isnan(right_tmp)) = [];
                %check if empty
                if(~isempty(right_tmp))
                    right_tmp = unique(right_tmp);
                    for j = 1:length(right_tmp)
                        pairings = [pairings;both_edges(i,2),right_tmp(j)];
                    end
                end
            end
        end


        pairings(1,:) = [];

        if(~isempty(pairings))
            unique_pairs = unique(pairings, 'rows');
            %go through each element and change on map
            for i = 1:size(unique_pairs,1)
                %find all IDs with the same number at the first or second column
                [r1, ~] = find(unique_pairs(:,1) == unique_pairs(i,1));
                [r2, ~] = find(unique_pairs(:,2) == unique_pairs(i,2));

                to_change = [unique_pairs(r1,2);unique_pairs(r2,1)];
                for j = 1:length(to_change)
                    new_map(new_map == to_change(j)) = unique_pairs(i,1);
                end

            end
        end

    
    end

end