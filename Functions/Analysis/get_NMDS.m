function [Y_dis,stress,disparities,spread,convex_hull_points,points_biome,location_biome ] = get_NMDS(net, monthly_maps,classes, n_clusters,months, ID_maps)
% Function that performs the Non-metric multi-dimensional scaling on a
% given dataset

%{
Parameters:
    net (Network): Trained Self-organizing map
    monthly_maps (matrix): 3D matrix of biomes with dimension m x 180 X 360
        where m is the number of time steps (e.g. months, seasons, annual) 
    classes (vector): Labels of the observations
    n_clusters (int): Number of clusters
    months (vector): Indeces in monthly_maps that should be considered
    ID_maps (matrix,optional): 3D matrix of IDs of observations with dimension 
        m x 180 X 360.

 Output:
    Y_dis (matrix): Position of observatins in the 2-dimensional space
    stress (float): Minimized stress of mapping
    disparities (vector): Monotonic transformation of the dissimilarities
    spread (matrix): Maximum spread in the convex hull in the 2-dimensional
        space in x- and y-direction
    convex_hull_points (matrix): Position of corners of convex hulls
    points_biome (matrix): Position of points of biomes in 2D space
    location_biome (matrix): Position of IDs in 2D space

%}

% =========================================================================
% Initialize the matrices to store the coordinates of the convex hul, the
% points found for each biome, and the colormap
% =========================================================================

    convex_hull_points = [NaN NaN NaN];
    points_biome = [NaN NaN NaN];
    
    cmap = morgenstemning(12);
    [cmap_tmp] = shuffle_colormap(cmap);
    [cmap_tmp1] = shuffle_colormap(cmap_tmp);


    cmap = ametrine(12);
    [cmap_tmp] = shuffle_colormap(cmap);
    [cmap_tmp2] = shuffle_colormap(cmap_tmp);

    cmap = isolum(12);
    [cmap_tmp] = shuffle_colormap(cmap);
    [cmap_tmp3] = shuffle_colormap(cmap_tmp);

    comb_cmap = [cmap_tmp2([2,1,3],:);cmap_tmp1(7,:);cmap_tmp3(5,:);cmap_tmp1(8,:);cmap_tmp1(9,:);cmap_tmp1(4,:)];

    %% test with month

    % =========================================================================
    % Create monthly maps of the associated neuron to each pixel, and shift the
    % Southern Hemisphere by six months
    % =========================================================================

    [neuron_maps] = prepare2plot( classes);
    tmp_corrected_neuron = neuron_maps*NaN;
    for i =1:12
        if(i<7)
            j = mod(i+6,13);
        else
            j = mod(i+6,13) +1;
        end
        %i and j are the indices of the months that need to be combined
        tmp_corrected_neuron(i,91:end,:) = neuron_maps(i,91:end,:);%raw_monthly_maps(i,91:end,:);%smooth_map(i,91:end,:);
        tmp_corrected_neuron(i,1:90,:) = neuron_maps(j,1:90,:); % raw_monthly_maps(j,1:90,:); %smooth_map(j,1:90,:); % 
    end
    corrected_neuron = tmp_corrected_neuron(months,:,:);


    % =========================================================================
    % Load in SOM data, calculate dissimilarities and assess goodness of NMDS
    % analysis with standard tools
    % =========================================================================

    data = net.IW{1};
    dissimilarities = pdist(data,'cityblock');

    opts = statset('Display','final');

    [Y_dis,stress,disparities]  = mdscale(dissimilarities,2);
    spread = NaN(n_clusters,2);
    figure
    hold on;

    distances = pdist(Y_dis,'cityblock');
    [~,ord] = sortrows([disparities(:) dissimilarities(:)]);
    plot(dissimilarities,distances,'bo', ...
         dissimilarities(ord),disparities(ord),'r.-', ...
         [0 450],[0 450],'k-')
    xlabel('Dissimilarities')
    ylabel('Distances/Disparities')
    legend({'Distances' 'Disparities' '1:1 Line'},...
           'Location','NorthWest');

    %% NMDS 

    % =========================================================================
    % NMDS part. Use the neuron maps from above and the seasonally corrected
    % monthly maps
    % =========================================================================
    axis tight manual % this ensures that getframe() returns a consistent size

    filename = 'NMDS.gif';
    months = {'Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec'};


    for m = 1:12

        figure
        fig = gcf;
        hold on;
        location_biome = [NaN NaN NaN];
        available_labels = unique(monthly_maps(m,~isnan(monthly_maps(m,:,:))));
        %Don't do cluster 9 as we do not analyze it!
        available_labels(available_labels == 9) = [];


            for ii = 1:length(available_labels)
                %get the label
                i = available_labels(ii);
                %get all the neurons found in one biome and store as tmp, do it
                %individually for every month, is slower, but needed to match
                %location_tmp (I think)
            %     tmp = NaN;
            %     for mm = 1:1%2
            %         tmp = [tmp;corrected_neuron(mm,monthly_maps(mm,:,:) == i)'];
            %     end
            %     tmp(1) = [];
                tmp = corrected_neuron(m,monthly_maps(m,:,:) == i)';
                if(nargin == 6)
                    %get a vector with the ID for the pixels
                    location_tmp = [NaN NaN];
                    for mm = 1:1%12
                        tmp2 = ID_maps(mm,monthly_maps(mm,:,:) == i); 
                        location_tmp = [location_tmp;[tmp2'.*0 + mm, tmp2']];
                    end
                    location_tmp(1,:) = [];
                end


                %get all the 2D values of all neurons transformed, and consider the
                %frequency of occurrence of each neuron (i.e. not only unique neurons)
                Y = Y_dis(tmp,:);
                % =====================================================================
                % Check that these neurons are really in the biome; UHE 19 AUG 2019
                % 
                % =====================================================================
            %     in_biome_vec =  [ID_vec(location_tmp(:,2),2:3) location_tmp(:,1), location_tmp(:,1).*0 + 1];
            %     [in_biome_mat] = prepare2plotV2( in_biome_vec);
            %     in_biome_mat = sum(in_biome_mat,1,'omitnan');
            %     plotSOM(in_biome_mat,1,NaN);
            %     i
            %     jj = input('df')
                % =====================================================================
                % Code above works fine only neurons found in a biome are taken
                % =====================================================================


                %separate into components x and y  (first dimension and second
                %dimension)
                x = Y(:,1);
                y = Y(:,2);
                XY = [x,y];

                %get some thresholds to define envelope
                p95 = prctile(XY,90,1); 
                p5 = prctile(XY,10,1); 

%                 p95 = prctile(XY,100,1); 
%                 p5 = prctile(XY,1,1); 

                if(nargin == 6)
                    %get a vector with the ID for the pixels
                    location_tmp(XY(:,1) > p95(1) | XY(:,2) > p95(2) | XY(:,1) < p5(1) | XY(:,2) < p5(2),:) = [];
            %         location_tmp(XY(:,1) < p5(1) | XY(:,2) < p5(2),:) = [];
                    %store the found IDs into location biome, which contains the biome
                    %label, the month, and the ID in the columns
                    location_biome = [location_biome;[location_tmp(:,1).*0 + i,location_tmp]];

                    % =====================================================================
                    % Check again that these neurons are really in the biome; UHE 19
                    % Aug 2019
                    % 
                    % =====================================================================
            %         tmp_tmp = location_biome(location_biome(:,1) == i,2:3);
            %         in_biome_vec =  [ID_vec(tmp_tmp(:,2),2:3) tmp_tmp(:,1), tmp_tmp(:,1).*0 + 1];
            %         [in_biome_mat] = prepare2plotV2( in_biome_vec);
            %         in_biome_mat = sum(in_biome_mat,1,'omitnan');
            %         plotSOM(in_biome_mat,1,NaN);
            %         i
            %         jj = input('df')
                    % =====================================================================
                    % Code above still works fine only neurons found in a biome are taken
                    % =====================================================================


                end



                XY(XY(:,1) > p95(1) | XY(:,2) > p95(2) | XY(:,1) < p5(1) | XY(:,2) < p5(2),:) = [];
            %     XY(XY(:,1) < p5(1) | XY(:,2) < p5(2),:) = [];



                if(size(unique(XY,'rows'),1) > 2)
                    k = convhull(XY(:,1),XY(:,2));
                    h0(i) = plot(XY(k,1),XY(k,2),'LineWidth',4,'Color',comb_cmap(i,:));
                    convex_hull_points = [convex_hull_points ; [XY(k,1).*0 + i, XY(k,1), XY(k,2)]];        
                end 

                spread(i,1) = abs(max(XY(:,1))- min(XY(:,1)));
                spread(i,2) =  abs(max(XY(:,2))- min(XY(:,2)));
                %h1(i) = plot(x,y,'+','LineWidth',1,'Color',comb_cmap(i,:));
                h1(i) = plot(x,y,'o','LineWidth',1,'Color',comb_cmap(i,:),'LineWidth',4);
                h2(i) = plot(median(Y(:,1)),median(Y(:,2)),...
                            'k*','LineWidth',8); 

                points_biome = [points_biome;[XY(:,1).*0 + i, XY]];


            end
            xlabel('Dimension 1')
            ylabel('Dimension 2')
            title(months(m))
            %
            legend_names = {'(1) TRP','(2) HIL','(3) WIS','(4) SUS','(5) HIT ', '(6) MTR',...
                '(7) PEU','(8) SMN','Cluster Mean'}
            tmp = legend_names([available_labels,9]);
            legend([h1(available_labels), h2(1)],tmp);

            xlim([-200 200])
            ylim([-150 250])
            convex_hull_points(1,:) = [];
            points_biome(1,:) = [];
            if(nargin == 6)
                location_biome(1,:) = [];
            end
            set(fig, 'Position', get(0, 'Screensize'));
            drawnow 

            % Capture the plot as an image 
            frame = getframe(fig); 
            im = frame2im(frame); 
            [imind,cm] = rgb2ind(im,256); 

            % Write to the GIF File 
            if m == 1 
              imwrite(imind,cm,filename,'gif', 'Loopcount',inf); 
            else 
              imwrite(imind,cm,filename,'gif','WriteMode','append'); 
            end 


    end



end
