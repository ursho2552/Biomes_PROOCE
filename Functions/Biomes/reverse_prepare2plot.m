function [classes] = reverse_prepare2plot(new_map)
%Function converts map data to vector
%new_map contains the map (1x180x360)
%old_sequence contains the original/clustered classes


classes = [NaN NaN NaN NaN];
for m = 1:size(new_map,1)
    for lat = 180:-1:1
        for lon = 1:360
            classes = [classes;[lon, lat, m, new_map(m,lat,lon)]];
        end
    end
    m
end
classes(isnan(classes(:,end)),:) = [];




end
