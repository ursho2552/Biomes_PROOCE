%% Validation 

%==========================================================================
% In this script we first perform some dimensionality reduction using PCA,
% on the trained neurons, and then cluster the neurons into 2 to 100
% clusters, i.e. we associate each class label to a "cluster" label. We 
% then applay the trained SOM to the validation set, i.e. we determine the 
% class label for each observation in the validation set. We use the 
% previous clustering to determine the "cluster" label in the validation
% set and compute the error.
%==========================================================================

%load all variables needed
cd('/net/kryo/work/ursho/PhD/Projects/Biomes/Scripts/Biomes_PROOCE/Data/')
load('00Probabilities/Simple_sort_Data.mat')
load('HelpVariables.mat')

%define directory to store validation experiments
folder_results = '/net/kryo/work/ursho/PhD/Projects/Biomes/Scripts/Biomes_PROOCE/Data/03CrossValidation/Validation/';
%define directory where trained SOMs are stored
folder_SOM_data = '/net/kryo/work/ursho/PhD/Projects/Biomes/Scripts/Biomes_PROOCE/Data/03CrossValidation/SOMs/';

%define the "names" of leave-out data
fraqs = [10,20,30,50];
%define number of folds per leave-out experiment
Num_boundaries = [1,10;1,5;1,3;1,2];

%uppper bound on the number of clusters
n = 100;


for j = 1:length(fraqs)
    fr = fraqs(j);
    Num = Num_boundaries(j,2);
    %We tested different metrics but took the most used one in the end
    difference = NaN(1,Num,n);
    FOM = difference;
    difference_alt = difference;
    
    
    for ts = 1:Num
        cd(folder_SOM_data)
        stri = horzcat('CV_Fr_',int2str(fr),'_Num_',int2str(ts),'.mat');
        load(stri)


        %======================================================================
        % Dimensionality reduction using PCA 
        %======================================================================
        tic
        original_weights = noise_net.IW{1};
        %original_weights = bsxfun(@minus,original_weights,mean(original_weights));
        [~,~,latent,~,~] = pca(original_weights);

        %Kaiser's rule
        [r,c] = find(latent > 1); 
        clearvars latent
        %[coeff,score,latent,tsquared,explained]
        [coeff,score,~,~,~] = pca(original_weights,'NumComponents',r(end));
        toc

        %agglomerate trained neurons into 2 to 100 clusters
        dat_net = original_weights*coeff;%score;
        %cluster trained neurons into n_clusters clusters
        [R_training,T_matrix_training,~] = DaviesBouldinDendrogram('cityblock', dat_net,classes_noise(:,1),...
             2,0,'weighted',0); 

        %T_matrix_training contains for each column the associated cluster of each trained neuron    
        max_clust = length(R_training)-1;

        % =====================================================================
        % Get original distribution of neurons
        % =====================================================================


        ind = No_nan_phyto_simple(training_flag == 0,1); 
        ind_tr = No_nan_phyto_simple(training_flag == 1,1);

        tmp_simple = No_nan_phyto_simple(:,2:4);
        tmp_simple(ind,:) = [];

        [neuron_maps] = prepare2plot( [tmp_simple,classes_noise(:,1)]);

        [ validation_classes ] = Apply_SOM( No_nan_phyto_simple(ind,:), noise_net );
        
%             %Alternative        
%             Analog_data = No_nan_phyto_simple(ind,:);
%             validation_classes = Analog_data(:,1).*NaN;
% 
%             for iii = 1:size(Analog_data,1)
%                 tmp = nansum(abs(Analog_data(iii,5:end-1) - noise_net.IW{1}),2);
%                 [r,c] = find(tmp == min(tmp));
%                 validation_classes(iii) = r;
%                
%             end


        %choose number of cluster

        disp('Clustering...')

        for k = 1:99%from 2 to 100 (as defined by function DaviesBouldinDendrogram)
            %number of clusters
            n_clusters = k+1;
            if(n_clusters > max_clust)
                break;
            end

            % =================================================================
            % Calculate average cluster neuron by clustering the training
            % classes, and averaging using the frequency of occurrence of
            % each class
            % =================================================================
            [ new_training_classes] = reduce_classes( T_matrix_training(:,n_clusters-1),...
                   [tmp_simple(:,1:2),classes_noise] ); 

            % get map of clusters
            training_map_monthly = prepare2plot([tmp_simple,new_training_classes(:,end)]);

            %calculate weights

            weights = noise_net.IW{1};          
            new_weights = NaN(n_clusters,size(weights,2));

            for ii = 1:n_clusters
                tmp = neuron_maps(training_map_monthly == ii);
                new_weights(ii,:) = mean(weights(tmp,:),1,'omitnan');
            end

            % =================================================================
            % Cluster validation classes based on training
            % =================================================================

             [ new_validation_classes] = reduce_classes( T_matrix_training(:,n_clusters-1),...
                     [No_nan_phyto_simple(ind,2:3),validation_classes] );

            %==================================================================
            % Expand new_weights by one dummy neurons that is zero in all
            % places, use this one if new_validation classes is NaN
            % 18. Apr 2019
            %==================================================================
            new_weights = [new_weights;new_weights(1,:).*0];

            new_validation_classes(isnan(new_validation_classes(:,3)),3) = size(new_weights,1);

            obs_avg = No_nan_phyto_simple(ind,5:end-1) - new_weights(new_validation_classes(:,3),:);
            difference(1,ts,k) = mean(mean(abs(obs_avg),1,'omitnan'),'omitnan'); %fendereski
            difference_alt(1,ts,k) = sum(sum(abs(obs_avg),2,'omitnan'),'omitnan');%as quantization error
            FOM(1,ts,k) = sqrt(mean(sum(obs_avg.^2,2,'omitnan'),'omitnan'));%figure of merit
            %NOTE: magnitude of quantization error form and fendereski is
            %different, but slope change is equal. For FOM, the slope change is
            %less large, but similar --> problem of large dimensions aggarwal
            %et al. 

        end

        cd(folder_results)
        save(horzcat('CV_error_Fr_',int2str(fr)),'difference','FOM',...
        'difference_alt')
    end

end

    

%% Figure 1
% We use the Oliver approach: optimal number of classes is 
% the point where adding one more class decreases the error by less than 1%
% for three consecutive clusterings

Num_boundaries = [1,10;1,5;1,3;1,2];


 styles_line = ['s-','-o','+-','d-'];
 cv_tr = [10,8,9,9];
for test = 1:4
    Num = Num_boundaries(test,2);
    
    if Num == 10
        load('CV_error_Fr_10.mat')
    elseif Num == 5
        load('CV_error_Fr_20.mat')
    elseif Num == 3
        load('CV_error_Fr_30.mat')
    elseif Num == 2
        load('CV_error_Fr_50.mat')
    end
    
    difference(difference==0) = NaN;
    
    mean_difference = mean(squeeze(difference),1,'omitnan');

    dydx = mean_difference.*NaN;
    for i = 2:length(mean_difference)-1
        dydx(i) = (mean_difference(i+1) - mean_difference(i-1))./(2*nanmax(mean_difference));
    end

    thresh = abs(dydx(1:13));
        
    h = 2:length(mean_difference)+1; 
    figure()
    hold on;

    xlabel('Number of clusters')

    yyaxis left
    ylabel('Mean validation error') 
    plot(h,mean_difference);

    yyaxis right
    ylabel('Relative change of validation error') 
    plot(h,abs(dydx),'+-');

    set(gca,'XTick', h(1:5:99));
    hline = refline(0,0.01)
    hline.Color = 'k';

    hline.LineStyle = '-.';

    vline = line([cv_tr(test) cv_tr(test)], [0 0.16]);
    vline.Color = 'k';
    vline.LineStyle = '-.';
    grid on

    xlim([0 100]);

    % add line at tipping point
    vline = line([cv_tr(test) cv_tr(test)], [0 0.16]);
    vline.Color = 'k';
    vline.LineStyle = '-.';
    set(findall(gcf,'-property','FontSize'),'FontSize',25)
    set(findall(gcf,'-property','LineWidth'),'LineWidth',3)
    abs(dydx(1:15))    
       
end