function [connected_edges,num_pairings] = find_most_connected(edges,val_edges)
    

    %see if there are duplicate edges and delete them
    
    [edges, ia, ic] = unique(edges,'rows','stable');

    %create a graph
    G = graph(edges(:,1)', edges(:,2)');

    %bin connected components 
    Connected_comp = conncomp(G);
    %Connected_comp contains the labels for each pair

    %search for the most connected chain
    num_components = [unique(Connected_comp);unique(Connected_comp)*NaN];
    for i = 1:length(num_components)
        num_components(2,i) = length(find(Connected_comp == num_components(1,i)));
    end
    %num_components contains the number of pairs for each label
    
        
    %search for maximum
    r = find(num_components(2,:) == max(num_components(2,:)));
    if(length(r) == 1)
        %in Connected_comp find the group r
        chain = find(Connected_comp == r);
        %get the pairs from edges
        col1 = ismember(edges(:,1),chain);
        %col2 = ismember(edges(:,2),chain)

        edges(col1 == 0,:,:) = [];
        connected_edges = edges;
        num_pairings = size(connected_edges,1);
    elseif(length(r) > 1)
        %get the value in the graph and pick the largest one
        val_sub = r.*NaN;
        for i = 1:length(r)
            
             chain = find(Connected_comp == r(i));
            %get the pairs from edges
            col1 = ismember(edges(:,1),chain);
            val_sub(i) = mean(val_edges(col1 == 1),'omitnan');
        end
        [rr, cc] = find(val_sub == max(val_sub,[],'omitnan'));
        
            
        chain = find(Connected_comp == r(rr));
        %get the pairs from edges
        col1 = ismember(edges(:,1),chain);
        %col2 = ismember(edges(:,2),chain)

        edges(col1 == 0,:,:) = [];
        connected_edges = edges;
        num_pairings = size(connected_edges,1);
    else
        num_pairings = 0;
        connected_edges = [NaN NaN NaN];
    end

end