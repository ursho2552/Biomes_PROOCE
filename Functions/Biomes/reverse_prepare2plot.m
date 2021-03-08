function [classes] = reverse_prepare2plot(new_map)
%Function converts map data to vectorized data

%{
Parameters:
    new_map (matrix): 3D matrix of with map of biomes (m x 180 x 360)

 Output:
    classes (vetor): Vector of biome labels

%}


classes = [NaN NaN NaN NaN];
for m = 1:size(new_map,1)
    for lat = 180:-1:1
        for lon = 1:360
            classes = [classes;[lon, lat, m, new_map(m,lat,lon)]];
        end
    end
end
classes(isnan(classes(:,end)),:) = [];




end
