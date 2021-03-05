%% Self-Organizing Map

% =========================================================================
% This main script trains the SOMs with different number of neurons, and
% epochs and finds the "optimal" setup out of all the tested setups. First
% we test for the optimal number of neurons, and then for the optimal
% number of epochs. For the former we use the standard number of epochs of
% 200.
% =========================================================================

%load dataset
cd('/net/kryo/work/ursho/Damiano_Presence_data/presence_absence_tables_ensemble_averages/Group_specific_background_approach/Data/00Probabilities/')
load('Simple_sort_Data.mat')

% =========================================================================
% For each month count the number of observations, i.e. the number of 
% 1°x1°-pixels, and not the number of features, i.e. number of species 
% modeled.
% =========================================================================

for i =1:12
    [r, ~] = find(No_nan_phyto_simple(:,4) == i);
    range_phyto(i) = length(r);
%     mm_bio(i) = r(end)+1;
end

% =========================================================================
% Define the total number of neurons based on the average value in 
% range_phyto using the rule of thumb as presented in Vesanto and Alhoniemi
% (2000). 
% =========================================================================

bounds = mean(5 * sqrt(range_phyto));

% =========================================================================
% Calculate the square root to get one side (i.e. assume the neurons are
% positioned in a squared matrix. Be aware that there are other setups out
% there.
% =========================================================================

approx = sqrt(bounds);
approx = round(approx);

% =========================================================================
% I defined a range of number of neurons for testing the
% stability/robustness/accuracy of our choice
% =========================================================================

min_approx = approx(1)-floor(approx(1)/4);
max_approx = approx(1)+floor(approx(1)/4);
num_neurons = [(5:3:min_approx-2),(min_approx+1:2:approx),(approx+2:2:max_approx),max_approx+14:20:70];
num_neurons = [num_neurons; num_neurons]';



%load data needed
cd('/net/kryo/work/ursho/Damiano_Presence_data/presence_absence_tables_ensemble_averages/Group_specific_background_approach/Data/00Probabilities/')
load('Transformed_CompleteSuitePhyto.mat')


%allocate matrices to store the computational time (might be interesting to
%see if a setup is faster or not
computing_time = error_neurons;

for i = 1:size(num_neurons,1)
    
        %get dimensions of neuron lattice    
        d1 = num_neurons(i,1);
        d2 = num_neurons(i,2);
        tic
        %train SOM
        [classes, net] = My_SOM( Transformed_phyto, d1,d2, 200,'mandist' );
        
        computing_time(i) = toc
        cd('/net/kryo/work/ursho/Damiano_Presence_data/presence_absence_tables_ensemble_averages/Group_specific_background_approach/Data/01NeuronsError')
        save(horzcat('Single_run',num2str(i)),'classes','net','i')
        
        %calculate error
        tic
        [ qe, te, total_error ] = get_total_error( Transformed_phyto,classes, net,'mandist' ); 
        toc
        
        %store qe and te to later compare their fraction to the total error
        %save the error
        save(horzcat('Single_run_error',num2str(i)),'total_error','qe','te','i','computing_time')

end
%save(horzcat('Single_run',num2str(i)),'total_error','qe','te','classes','net','i','computing_time')

%% Plot error-metric as a function of number of neurons

%merge single runs into a single one
error_neurons = ones(1,size(num_neurons,1)).*NaN;
for i=1:size(num_neurons,1)
    load(horzcat(file,'Single_run_',int2str(i),'_error.mat'));
    error_neurons(i) = total_error;
end


%find the best lattice dimension based on the figure
dim1= strcat(int2str(num_neurons(:,1)),' x ',int2str(num_neurons(:,2)));
label = dim1;
n4 = size(num_neurons,1);
h = num_neurons(:,1).*num_neurons(:,2)
x = 1:n4;
y = error_neurons(1:end);
%calculate first derivative
dydx = y*NaN;
for i =2:length(y)-1
dydx(i) = (y(i+1) - y(i-1))/(h(i+1)-h(i-1))
end


figure()
hold on;
yyaxis left
xlabel('Number of Neurons')
plot(h,y)
ylabel('Total error')
yyaxis right
plot(h,dydx)
ylabel('Total error change')
set(gca,'XTick', h)
set(gca,'XTickLabel',[label])
hold off;

% Optimal dimension appears to be 31 x 31 neurons
optimal_dim = [31,31];

%% Find optimal number of epochs

epoch = [1, 5, 10, 20, 50, 100, 200, 300, 400, 500, 700, 1000];
error_epoch =ones(size(epoch,2),1)*NaN;

computing_time_ep = zeros(1,size(epoch,2));

for i = 1:size(epoch,2)
        
        d1 = optimal_dim(1);
        d2 = optimal_dim(2);
        tic
        [classes_ep, net_ep] = My_SOM( Transformed_phyto, d1,d2, epoch(i),'mandist' );   
        computing_time_ep(i) = toc
        %save each SOM run
        cd('/net/kryo/work/ursho/Damiano_Presence_data/presence_absence_tables_ensemble_averages/Group_specific_background_approach/Data/02EpochError')
        save(horzcat('Single_run_ep',num2str(i)),'classes_ep','net_ep','i')
        tic
        [ qe, te, total_error_ep ] = get_total_error( Transformed_phyto,classes_ep, net_ep,'mandist' ); 
        toc
        error_epoch(i) = total_error_ep;
        
        %store qe and te to later compare their fraction to the total error
        %save the error
        save(horzcat('Single_run_error_ep',num2str(i)),'total_error_ep','qe','te','i','computing_time_ep')
    
end

x = epoch(1:end);
y = error_epoch(1:end);

figure()
hold on;
plot(x,y)
xlabel('Epochs')
ylabel('Total error')
hold off;

%The optimal is 200
optimal_epoch = 200;

%% Optimal setup run
%==========================================================================
% Train final SOM with optimal dimension and epoch if not already done
% In our case this is not needed since optimal epoch is 200 and we already 
% calculated this above!!
%==========================================================================


% tic
% d1 = optimal_dim(1);
% d2 = optimal_dim(2);
% [optimal_classes, optimal_net] = My_SOM( Transformed_phyto, d1,d2, optimal_epoch,'mandist' );
% toc
% %copy classes from SOM into the lat-lon combination
% [ qe, te, total_error ] = get_total_error( Transformed_phyto,optimal_classes, optimal_net,'mandist' )
% classes_month = [Transformed_phyto(:,2:3),optimal_classes]; %lon|lat|class
% new_data_ann = prepare2plot(classes_month,mm_bio,LatLon,1,0);




