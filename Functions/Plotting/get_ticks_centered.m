function [positions] = get_ticks_centered(n_classes)
% Function determines the middle point for each class in the colorbar. The 
% colorbar has limits 1 to n_classes.

%{
Parameters:
    n_classes (int): Label of the biome with the largest magnitude
 
 Output:
    positions (matrix): Matrix encoding the positions of the ticks for each
        label

%}

    %calculate the steps from, i.e. how large of an area does a label occupy
    step = (n_classes-1)/n_classes;
    %step = round(100*(n_classes/(n_classes+1))/2)/100;

    %calculate the boundaries for each label starting from 1
    tmp_positions = 1;
    for i = 1:n_classes
        tmp_positions = [tmp_positions,tmp_positions(end)+step];
    end

    %calculate positions as the mid point between the boundaries
    positions = NaN(1,n_classes);
    for i =1:length(tmp_positions)-1
        positions(i) = (tmp_positions(i) + tmp_positions(i+1))/2
    end



end