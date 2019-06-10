% Local Feature Stencil Code
% Written by James Hays for CS 143 @ Brown / CS 4476/6476 @ Georgia Tech

% 'features1' and 'features2' are the n x feature dimensionality features
%   from the two images.
% If you want to include geometric verification in this stage, you can add
% the x and y locations of the features as additional inputs.
%
% 'matches' is a k x 2 matrix, where k is the number of matches. The first
%   column is an index in features1, the second column is an index
%   in features2. 
% 'Confidences' is a k x 1 matrix with a real valued confidence for every
%   match.
% 'matches' and 'confidences' can empty, e.g. 0x2 and 0x1.
function [matches, confidences] = match_features(features1, features2)

% This function does not need to be symmetric (e.g. it can produce
% different numbers of matches depending on the order of the arguments).

% To start with, simply implement the "ratio test", equation 4.18 in
% section 4.1.3 of Szeliski. For extra credit you can implement various
% forms of spatial verification of matches.

% Placeholder that you can delete. Random matches and confidences
% num_features = min(size(features1, 1), size(features2,1));
% matches = zeros(num_features, 2);
% matches(:,1) = randperm(num_features); 
% matches(:,2) = randperm(num_features);
% confidences = rand(num_features,1);
%% Manually Computing Euclidean Distance
features1_size = size(features1);
features2_size = size(features2);

matches = [];
confidences = [];
features11 = features1(:,:,1);
features12 = features1(:,:,2);
features21 = features2(:,:,1);
features22 = features2(:,:,2);

for i = 1:features1_size(1)
    distances11 =  bsxfun(@minus,features21, features11(i,:));
    distances11 = distances11.^2;
    distances11 = sqrt(sum(distances11'));
    
    distances21 =  bsxfun(@minus,features22, features11(i,:));
    distances21 = distances21.^2;
    distances21 = sqrt(sum(distances21'));
    
  
    distances12 =  bsxfun(@minus,features21, features12(i,:));
    distances12 = distances12.^2;
    distances12 = sqrt(sum(distances12'));
    
    distances22 =  bsxfun(@minus,features22, features12(i,:));
    distances22 = distances22.^2;
    distances22 = sqrt(sum(distances22'));

    
    distances = [distances11 distances12 distances21 distances22];
    [distances, index] = sort(distances);
    
    distances(distances==0) = 1e-7;
    
    NNDR = distances(1)/distances(2);
    if (NNDR>0.9)
        continue;
    end
    
    in_f2 = mod(index(1),features2_size(1));
    if in_f2==0
        in_f2 = features2_size(1);
    end
    matches = [matches; [i, in_f2]];
    confidences = [confidences; (1-NNDR)];
    
end
%% KNNSearch Algorithm
% 
% features1_size = size(features1);
% 
% matches = [];
% confidences = [];
% features11 = features1(:,:,1);
% features21 = features2(:,:,1);
% 
% [two_indices1, two_distances1] = knnsearch(features21, features11, 'k',2);
% 
% for i = 1:features1_size(1)
% 
%     NNDR = two_distances1(i,1)/two_distances1(i,2);
%     if (NNDR>0.9)
%         continue;
%     end
%     
%     in_f2 = two_indices1(i,1);
%     matches = [matches; [i, in_f2]];
%     confidences = [confidences; (1-NNDR)];
% end
%%
% Sort the matches so that the most confident onces are at the top of the
% list. You should not delete this, so that the evaluation
% functions can be run on the top matches easily.
[confidences, ind] = sort(confidences, 'descend');
matches = matches(ind,:);