function [sequence] = find_most_similar(D, Z,starter,maxclust, sequence)
%Function finds the most similar neuron for the starter

%get distance from starter to all other neurons
% D = get_distance_from_weights(weights,starter);
% %construct linkage matrix
% Z = get_linkage_from_weights(weights,maxclust,starter);

%read sequence from most similar to least similar according to the linkage
%matrix


%continue until sequence is large enough
%while(length(sequence) < size(D,1)+1)

%find starter

    [r c] = find(Z(:,1:2) == starter);




    temp_solution = Z(r,:);
    temp_solution(1,c) = NaN;
    temp_solution = temp_solution(~isnan(temp_solution));

    if(temp_solution(1) <= maxclust && ~ismember(temp_solution(1),sequence))
        sequence = [sequence;temp_solution(1)];
        %search for next
        starter = temp_solution(2);
        if(starter <= maxclust)
            sequence = [sequence;starter];
        else
            [sequence] = find_most_similar(D, Z,starter,maxclust, sequence);
        end
    elseif(temp_solution(1) > maxclust)
        % search for all rows with temp_solution(1),

        [r1 c1] = find(Z(:,3) == temp_solution(1));
        tmp_Z = Z(r1,1:2);
        clust_seq = NaN;
        clust_seq = [clust_seq;reshape(tmp_Z(tmp_Z <= maxclust),[],1)];
        tmp_Z(tmp_Z <= maxclust) = [];
        all_found = 0;
        while(~all_found)
            for i = 1:length(tmp_Z)
                [r2 c2] = find(Z == tmp_Z(i));

                clust_seq = [clust_seq;reshape(Z(r2,:),[],1)];
            end
            clust_seq = unique(clust_seq);
            clust_seq(end) = [];

            %check that no new clusters appear other than those listed in tmp_Z



            clust_seq(clust_seq > max(tmp_Z,[],'omitnan')) = [];
            clust_seq(clust_seq == temp_solution(1)) = [];
            clust_seq(ismember(clust_seq,tmp_Z)) = [];
            clust_seq = unique(clust_seq);
            if(~isempty(clust_seq(clust_seq > maxclust))) %(length(clust_seq(clust_seq > maxclust)) > 0) %
                tmp_Z = [tmp_Z, reshape(clust_seq(clust_seq > maxclust),1,[])];
            else
                all_found = 1;

                clust_seq(clust_seq > max(tmp_Z,[],'omitnan')) = [];
                clust_seq(clust_seq == temp_solution(1)) = [];
                clust_seq = unique(clust_seq);
            end
        end




        %check distances to determine sequence
        D_tmp = D(clust_seq,:);
        D_tmp = sort(D_tmp,2,'ascend');
        sequence = [sequence;D_tmp(:,1)];
        starter = temp_solution(2);
        if(starter <= maxclust)
            sequence = [sequence;starter];
        else
            [sequence] = find_most_similar(D, Z,starter,maxclust, sequence);
        end



    end

        
end
                   
                    
            
            
            
 