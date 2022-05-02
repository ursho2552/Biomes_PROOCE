%% Create 10 partitionings for all leave out ratios

% =========================================================================
% This script divides the data into folds for the cross validation test. We
% use a seed for reproducibility
% =========================================================================


%load data
cd('/net/kryo/work/ursho/PhD/Projects/Biomes/Scripts/Biomes_PROOCE/Data/00Probabilities')
load('Simple_sort_Data.mat')



N = size(No_nan_phyto_simple,1);
%define the fractions. Note that these are only "names" used and not all
%accurately define the fraction of data used, e.g. 30 does not stand for
%30% but for 1/3 of the data or 33.333....%
fraqs = [10,20,30,50];

%choose a seed so that experiment can be repeated
seed = 7
   
for f = 1:4
    fr = fraqs(f);
    if(fr == 1)
        disp('1/100...')
        rng(seed)
        part = cvpartition(N,'KFold',100);
        idxtrain = [training(part,1),training(part,2),training(part,3),...
            training(part,4),training(part,5),training(part,6),...
            training(part,7),training(part,8),training(part,9),training(part,10)];
    elseif(fr == 5)
        disp('1/20...')
        rng(seed)
        part = cvpartition(N,'KFold',20);
        idxtrain = [training(part,1),training(part,2),training(part,3),...
            training(part,4),training(part,5),training(part,6),...
            training(part,7),training(part,8),training(part,9),training(part,10)];
    elseif(fr == 10)
        disp('1/10...')
        rng(seed)
        part = cvpartition(N,'KFold',10);
        idxtrain = [training(part,1),training(part,2),training(part,3),...
            training(part,4),training(part,5),training(part,6),...
            training(part,7),training(part,8),training(part,9),training(part,10)];
    elseif(fr == 20)
        disp('1/5...')
        rng(seed)
        part1 = cvpartition(N,'KFold',5);


        idxtrain = [training(part1,1),training(part1,2),training(part1,3),...
            training(part1,4),training(part1,5)];

    elseif(fr == 30)
        disp('1/3...')
        rng(seed)
        part1 = cvpartition(N,'KFold',3);


        idxtrain = [training(part1,1),training(part1,2),training(part1,3)];

    elseif(fr == 50)
        disp('1/2...')
        rng(seed)
        part1 = cvpartition(N,'KFold',2);


        idxtrain = [training(part1,1),training(part1,2)];


    end
    cd('/net/kryo/work/ursho/PhD/Projects/Biomes/Scripts/Biomes_PROOCE/Data/03CrossValidation/Folds')
    save(horzcat('Data_partitioning_cross_validation_fr_',int2str(fr),'_Seed_',int2str(seed)),'idxtrain');
end

