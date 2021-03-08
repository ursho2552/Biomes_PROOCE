function [mat] = get_observations_per_biome(env_map, biome_map, n_clusters)
%get matrix with all values for an environmental parameter found in each
%biome. 

%Prepare matrix that will contain the average values. rows are the
%different env parameters, row are the different clusters
mat = ones(180*360,n_clusters)*NaN;


    for j = 1:n_clusters
        tmp_env = env_map(1,biome_map(1,:,:) == j);
        
        for i =1:length(tmp_env(~isnan(tmp_env)))
            mat(i,j) = tmp_env(i);
        end
    end


end