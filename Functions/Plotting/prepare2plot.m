function [map_classes] = prepare2plot(classes)
%this function takes the information on location and time to create a
%matrix for plotting in a world map

%{
Parameters:
    classes (matrix): nx4 matrix containing Longitude_index|Latitude_index|Month|Class label

 
 Output:
    map_classes (matrix): mx180x360 matrix mapping the classes on a world
    map for each month m

%}

    %get the different months in the data
    index = unique(classes(:,3));

    map_classes = NaN(length(index),180,360);
    for m = 1:length(index)
        tmp = classes(classes(:,3) == index(m),:);
        for i = 1:length(tmp)
            map_classes(index(m),tmp(i,2),tmp(i,1)) = tmp(i,4);
        end
    end

    

end

