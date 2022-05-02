function [correspondence, maps] = get_correspondence_partitionings( partitioning1,partitioning2,area_map)
% Function Produces a correspondence matrix between two partitionings. 

%{
Parameters:
    partitioning1 (matrix): Matrix containing first partitioning
    partitioning2 (matrix): Matrix containing second partitioning
    area_map (matrix): Matrix containing area of each 1° pixel

 
 Output:
    correspondence (matrix): Matrix containing the correspondence between
    the first and the second partitioning
    maps (matrix): Matrix containing shared areas

%}
n = unique(partitioning1(~isnan(partitioning1)));
m = unique(partitioning2(~isnan(partitioning2)));

[t,y,x] = size(partitioning1);
if size(area_map,1) < t
    tmp = ones(t,y,x).*NaN;
    for i =1:t
        tmp(i,:,:) = area_map;
    end
end
area_map = tmp;
clear tmp

coverage = ones(length(n),length(m)).*NaN;

for i = 1:length(n)
    for j = 1:length(m)
        coverage(i,j) = sum(area_map(partitioning1 == n(i) & partitioning2 == m(j)));%/sum(area_map(partitioning1 == n(i)));
    end
end
        

%find minimum, until everything is 0 or NaN;
coverage(coverage == 0) = NaN;
tmp_coverage = coverage;

correspondence = [NaN NaN];
while(sum(sum(isnan(tmp_coverage))) < length(n)*length(m))

    [r c] = find(tmp_coverage == max(max(tmp_coverage,[],'omitnan'),[],'omitnan'));
    
    for i = 1:length(r)
        correspondence = [correspondence;[n(r(i)),m(c(i))]];
        tmp_coverage(r(i),:) = NaN;
        tmp_coverage(:,c(i)) = NaN;
    end

    
    
end
maps = partitioning1.*NaN;
correspondence(1,:) = [];

for i = 1:size(correspondence,1)
    maps(partitioning1 == correspondence(i,1) & partitioning2 == correspondence(i,2)) = correspondence(i,2);
end

% plotSOM(maps,1,NaN)

for i = 1:size(correspondence,1)
    maps(partitioning1 == correspondence(i,1) & partitioning2 == correspondence(i,2)) = correspondence(i,1);
end

% plotSOM(maps,1,9)

end