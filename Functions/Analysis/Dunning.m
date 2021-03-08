%% Who with whom where?

%look at all observations in a biome

%for each observation compare it to all others --> Sorensen-Dice
%coefficient

%take average value over all observations

%based on: Accurate Methods for the Statistics of Surprise and Coincidence
%by Ted Dunning

% for each two species compute the score based on the checkerboard approach
% put the score in a N by N matrix, where N is the total number of species

%table looks as follows:
%{ 
     ______________________________
    | A and B    | not A and B     |
    |____________|_________________|
    | A and not B| not A and not B |
    |____________|_________________|

where you count how many times it occurs per biome
%}


function [Score] = Dunning(presence)
    %how many species
    N = size(presence,2);
    N_tokens = size(presence,1);
    %prepare matrix
    Score = NaN(N);
    %go over each species
    for i = 1:N-1
        %since symmetrical only from i to N
        for j = i+1:N
            %take ith column and jth column
            dat = [presence(:,i),presence(:,j)];
            %create checkerboard
            
            check = NaN(2,2);            
            check(1,1) = length(find(dat(:,1) == 1 & dat(:,2) == 1));
            check(1,2) = length(find(dat(:,1) == 0 & dat(:,2) == 1));
            check(2,1) = length(find(dat(:,1) == 1 & dat(:,2) == 0));
            check(2,2) = length(find(dat(:,1) == 0 & dat(:,2) == 0));
                         
            %calculate log-likelihood ratio score
            k = reshape(check,[1 4]);
           
            Observed = check(1,1);
            Expected = length(find(dat(:,1) == 1))*length(find(dat(:,2) == 1))/N_tokens;
            tmp = 2 * sum(k,'omitnan') * (Shannon(k) -...
                Shannon(sum(check,2,'omitnan')) - Shannon(sum(check,1,'omitnan')));
            Score(i,j) = tmp;
            
            %Add Evert et al 2008 to distinguish between positive and negative associations
            %Absolute value expresses significance
            if(Observed < Expected)

                Score(i,j) = -tmp;
            end

            Score(j,i) = Score(i,j);
            
        end
    end
    
    function [H] = Shannon(vec)
        %input has to be a vector       
        n = sum(vec,'omitnan');
        H = sum(vec./n .* log(vec./n + (vec == 0)),'omitnan');
        
            
    end


end

 
 
 
 
 
 
 
 
 