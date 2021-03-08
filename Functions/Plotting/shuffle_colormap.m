function [cmap_tmp] = shuffle_colormap(cmap)
% Function used to shuffle colormaps and try to get colors that differ most
% in case the steps are too similar to each other

%{
Parameters:
    cmap (matrix): Matrix encoding the colormap
 
 Output:
    cmap_tmp (matrix): Matrix encodeing the shuffled colormap

%}
    cmap_tmp = cmap;
    counter_l = 1;
    counter_h = size(cmap,1);
    for i = 1:size(cmap,1)

        if(mod(i,2) == 1) %if odd
            cmap_tmp(i,:) = cmap(counter_l,:);
            counter_l = counter_l + 1;


        else

            cmap_tmp(i,:) = cmap(counter_h,:);
            counter_h = counter_h -1;
        end
    end
end