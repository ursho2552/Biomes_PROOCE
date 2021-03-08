function kappa = cohensKappa(original_map, new_map)
% Function to calculate the area weighted Kappa index between two maps

%{
Parameters:
    original_map (matrix): Map data with dimension n x 180 x 360
    new_map (matrix): Map data with dimension n x 180 x 360
   
 Output:
    kappa (float): Kappa index

%}
    %load area map
    cd('/net/kryo/work/ursho/Damiano_Presence_data/presence_absence_tables_ensemble_averages/Group_specific_background_approach/Data/')
    load('Area_map.mat')
    %transform to vector
    yhat = reshape(new_map,[],1);
    y = reshape(original_map,[],1);
    
    %get weights (area)
    months = size(original_map,1);
    extended_area = ones(months,180,360).*NaN;
    for m=1:months
        extended_area(m,:,:) = area_map;
    end
    %Calculate confusion matrix with area values
    min_label = min(unique([unique(original_map(~isnan(original_map)));unique(new_map(~isnan(new_map)))]));

    max_label = max(unique([unique(original_map(~isnan(original_map)));unique(new_map(~isnan(new_map)))]));
    labels = min_label:max_label;
    
    n = length(labels);
    
%     other_labels = unique(new_map(~isnan(new_map)));
%     other_labels
%     n1 = min(length(other_labels),n);
    C_area = ones(n,n).*0;
    
    for i=1:n
        for j=1:n
            C_area(i,j) = sum(extended_area(original_map==labels(i) & new_map == labels(j)));
        end
    end
    Area_tot = sum(extended_area(~isnan(original_map)));
    C_area = C_area./Area_tot;
    
    C = ones(n,n).*0;
    for i=1:n
        for j=1:n
            C(i,j) = length(original_map(original_map == labels(i) & new_map == labels(j)));
        end
    end

    %Calculate confusion matrix without area
%     C = confusionmat(y, yhat); % compute confusion matrix
%     C
%     C_area
    %Multiply C with C_area to get an area weighted Kappa index
    C = C.*C_area;
    
    %Calculate overlap relative to randomness
    
    n = sum(C(:)); % get total N
    C = C./n; % Convert confusion matrix counts to proportion of n
    r = sum(C,2); % row sum
    s = sum(C); % column sum
    expected = r*s; % expected proportion for random agree
    po = sum(diag(C)); % Observed proportion correct
    pe = sum(diag(expected)); % Proportion correct expected
    kappa = (po-pe)/(1-pe); % Cohen's kappa
    
end


%% Test behaviour
% cohensKappa(orig_map,changed_data)
% 
% cohensKappa(orig_map,orig_map)
% 
% tmp_o = orig_map.*-1;
% [r,c] = find(squeeze(orig_map(1,:,:)) == 1);
% % tmp_o(1,r(1),c(1)) = 1;
% cohensKappa(orig_map,tmp_o)
