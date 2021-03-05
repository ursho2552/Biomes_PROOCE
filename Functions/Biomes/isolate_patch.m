function [patch_map,rr]  = isolate_patch(old_map,label, area_map, tot_area, fraction, conn)
%Function labels all patches and isolates them

tot = tot_area;
%define a threshold of the minimum area for a patch
threshold = tot*fraction/100;
 
    %initialize a variable containing the current highest label, since for
    %our definition a biome with the same label in another ocean sector
    %constitutes a "seperate biome"
    highest = 0;
    %initialize a matrix that will store the resulting patches
    patch_map = old_map*NaN;
    %loop over each ocean basin (Indo-Pacific and Atlantic)
    for h = 1%:2

        biome = old_map;

        %isolate a cluster/biome label and transform it to binary
        biome(biome ~= label) = 0;
        biome(biome == label) = 1;
        
        %biome now only has 1 and 0, add a NaN where land should be

        biome = squeeze(biome);
        biome(isnan(biome)) = 0;
%       label the different patches
%         
        [patches,~] = bwlabel(biome,conn);
         
%         size(patches)
%         figure
%         hold on
%         surf(squeeze(patches))
%         jj = input('dfd')
%         
%         % add land masses back
%         tmp = old_map*NaN;
%         tmp(1,:,:) = patches;
%         [patches, ~] = add_remove_land(tmp,index_map,Atl_boundary, 2);
%         patches = squeeze(patches);
        
        
        %compare edges since sphere and not flat surface. Specially for the
        %Pacific sector since in our matrix it is divided
        tmp_biome = old_map*NaN; %this is a help variable that was used
        %in the previous version
        tmp_biome(isnan(old_map)) = 1; 
        
        
        [patches,~] = compare_edges(patches,tmp_biome,1,conn); %use old version
        %for further analysis, and making sure that the results are what we
        %intended (plotting) change the 0's to missing values (looks nicer)
%         
%            patches_tmp_tmp = patches;
%            patches_tmp_tmp(patches > 5) = 0;
%         figure
%         hold on
%         surf(squeeze(patches_tmp_tmp))
%         jj = input('dfd')
        
        
        patches(patches == 0) = NaN;
        %get the unique values of the defined patch labels
%         figure
%         hold on
%         surf(squeeze(patches))
%         jj = input('dfd')
        labels_patch = unique(patches);
        labels_patch(isnan(labels_patch)) = [];
        %only proceed if there are patches in the current ocean sector
        if(~isempty(labels_patch))
            %get cases, since Pacific is first, define a new highest label    
            if(h == 1)
                highest = labels_patch(end);
            %for all other cases after Pacific sector, change the values of
            %the patch labels to be larger than those of the previous
            %analyzed sector
            else
                for j = 1:length(labels_patch)
                    patches(patches == labels_patch(j)) = labels_patch(j) + highest;
                end
                %define a new highest value
                highest = unique(patches(~isnan(patches)));
                highest = highest(end);
            end
        end
        %add the values of each ocean sector to the matrix defined above
        %that by the end of the loop over the ocean sectors should contain
        %all patches
        patch_map(1,~isnan(patches)) = patches(~isnan(patches));
    end
    %for plotting and checking if the interim results are as
    %expected/"correct"
%     patch_map(isnan(patch_map)) = 0;
%     patches = squeeze(patch_map);
%     patch_map(1,:,:) = patches;
%     patch_map(patch_map == 0) = NaN; 
    %whole_map now only contains the patches
    %plotSOM(patch_map,13,latchl,lonchl,NaN)
    
    labels_whole = unique(patch_map(~isnan(patch_map)));    
    %initialize vector that stores the area and the label
    size_whole = ones(length(labels_whole),2)*NaN;
    size_whole(:,1) = labels_whole;
    %loop over each label in whole_map and store the size (area)

    for j = 1:length(labels_whole)
        %get area of each patch
        size_whole(j,2) = sum(area_map(patch_map == labels_whole(j)));
    end
    %find the patches that are below the threshold
    r = find(size_whole(:,2) < threshold);
    r_not = find(size_whole(:,2) >= threshold);
    if(isempty(r))
        patch_map = patch_map*NaN;
        rr = NaN;
    else
        rr = labels_whole(r);
        %delete all patches that are not found in rr
        if(~isempty(r_not))
            for j = 1:length(r_not)
                patch_map(patch_map == labels_whole(r_not(j))) = NaN;
            end
        end
        

        
    end
    
end    
    