function [metric_1,metric_2,metric_3,overlap_pairs1,overlap_pairs2] = compare_overlap(ref_map,new_map,area_map,ref_weights,leak_weights)
%function [area_correspondence,mean_corr] = compare_overlap(ref_weights,ref_map, new_weights,new_map,area_map)

%% initialization part
%use spatial correspondence instead of weights
%get labels in ref_map
unique_ref = unique(ref_map(~isnan(ref_map)));
%get labels in new_map
unique_new = unique(new_map(~isnan(new_map)));
%get time scale of map, monthly (12), seasonal (3), annual (1)
mm = size(ref_map,1);
%initialize area coverage matrix
area_corr = zeros(mm,length(unique_ref),length(unique_new));

%Get difference between centroids

D = pdist2(ref_weights,leak_weights);


%Get spatial overlap
for m = 1:mm
    %for each label in the reference map
    for i = 1:length(unique_ref)
        %get the shared area coverage with labels in new_map
        for j = 1:length(unique_new)
            area_corr(m,i,j) = sum(area_map(new_map(m,:,:) == unique_new(j) & ref_map(m,:,:) == unique_ref(i)));%/sum(area_map( ref_map(m,:,:) == unique_ref(i)));
            
        end
    end
end
% get the sum (not mean) over all months, and convert to 2D matrix
area_corr = squeeze(sum(area_corr,1,'omitnan'));

tmp_area_corr = area_corr;
tmp_D = D;
overlap_pairs1= [unique_ref,unique_ref.*NaN];
overlap_pairs2 = [NaN,NaN];%overlap_pairs1;


%calculate both correspondences
map_corr_1 = ref_map.*0;
map_corr_2 = map_corr_1;
for i = 1:length(unique_ref)
%find the maximum value of the sum of shared area
    if(max(max(tmp_area_corr)) ~= 0)
        [r1, c1] = find(tmp_area_corr == max(max(tmp_area_corr)));
        overlap_pairs1(r1,2) = c1;
        tmp_area_corr(r1,:) = NaN;
        tmp_area_corr(:,c1) = NaN;
        map_corr_1(ref_map == overlap_pairs1(r1,1) & new_map == overlap_pairs1(r1,2)) = 1;
    end   
    
    [r2, c2] = find(tmp_D == min(min(tmp_D)));
    overlap_pairs2 = [overlap_pairs2;[r2 c2]];
    %delete the row and column found, since we only want 1:1
    %correspondences
    tmp_D(r2,:) = NaN;
    tmp_D(:,c2) = NaN;

    if(~isempty(r2))
        map_corr_2(ref_map == r2 & new_map == c2) = 1;
    end

end
overlap_pairs2(1,:) = [];
[B, I] = sort(overlap_pairs2(:,1),'ascend');
overlap_pairs2 = overlap_pairs2(I,:);
%% Analysis part
% get overlap using either only spatial features or centroids


%stability of spatial features based on spatial features
%distances in weights based on spatial overlap
metric_1 = NaN;
overlap_pairs1(isnan(overlap_pairs1(:,2)),:) = [];
for i = 1:size(overlap_pairs1,1)
    metric_1 = [metric_1;D(overlap_pairs1(i,1),overlap_pairs1(i,2))];
end
metric_1 = mean(metric_1,'omitnan');


%stability of spatial features basen on weights
%compare the centroids, i.e. get correspondence vector 
monthly_tmp = NaN;
for m = 1:size(ref_map,1)
    monthly_tmp = [monthly_tmp; sum(area_map( map_corr_1(m,:,:) == 1))/sum(area_map(~isnan(ref_map(m,:,:))))];
end

metric_2 = mean(monthly_tmp,'omitnan');

%combined --> first weight correspondence, then spatial
metric_3 = NaN;
monthly_tmp = NaN;
for m = 1:size(ref_map,1)
    monthly_tmp = [monthly_tmp; sum(area_map( map_corr_2(m,:,:) == 1))/sum(area_map(~isnan(ref_map(m,:,:))))];
end

metric_3 = mean(monthly_tmp,'omitnan');
    
    
end
