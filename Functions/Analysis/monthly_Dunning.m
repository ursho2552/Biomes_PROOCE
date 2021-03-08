function [Scores] = monthly_Dunning(map,No_nan_phyto_simple,n_clusters)
% Calculate Dunning using monthly data

    data = map;

    %how many labels can be found in data?
    classes_lab = unique(data(~isnan(data)));
    %remove the ninth cluster from this analysis since it will not be
    %explored!!
    classes_lab(classes_lab==9) = [];


    %initialize matrix to store the scores
    Scores = NaN(n_clusters,12,size(No_nan_phyto_simple,2)-5,size(No_nan_phyto_simple,2)-5);
    % Scores = ones(n_clusters,size(No_nan_phyto_simple,2)-5,size(No_nan_phyto_simple,2)-5)*NaN;
    %Loop over each label and take all months together

    for i = 1:length(classes_lab)
        all_biome_phyto = ones(1,size(No_nan_phyto_simple,2))*NaN;
        for m = 1:size(data,1)
            %make a copy of the presence data for month m
            tmp_phyto = No_nan_phyto_simple(No_nan_phyto_simple(:,4) == m,:);
            %make a copy of the map
            tmp_map = data(m,:,:);
            %delete all other classes
            tmp_map(tmp_map ~= classes_lab(i)) = NaN;
            %get row and column of tmp_map where it is not missing
            [lat, lon] = find(squeeze(tmp_map) == classes_lab(i));
            %find the indeces in LatLon that are equal to lat and lon
            biome_phyto = zeros(length(lat),size(tmp_phyto,2));

            for j =1:length(lat)
                [r, c] = find(tmp_phyto(:,2) == lon(j) & tmp_phyto(:,3) == lat(j));
                biome_phyto(j,:) = tmp_phyto(r,:);

            end
            
            Scores(i,m,:,:) = Dunning(biome_phyto(:,5:end-1));
            m
        end
        i  
    end

end