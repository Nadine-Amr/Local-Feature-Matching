% Local Feature Stencil Code
% Written by James Hays for CS 143 @ Brown / CS 4476/6476 @ Georgia Tech

% Returns a set of feature descriptors for a given set of interest points. 

% 'image' can be grayscale or color, your choice.
% 'x' and 'y' are nx1 vectors of x and y coordinates of interest points.
%   The local features should be centered at x and y.
% 'descriptor_window_image_width', in pixels, is the local feature descriptor width. 
%   You can assume that descriptor_window_image_width will be a multiple of 4 
%   (i.e., every cell of your local SIFT-like feature will have an integer width and height).
% If you want to detect and describe features at multiple scales or
% particular orientations, then you can add input arguments.

% 'features' is the array of computed features. It should have the
%   following size: [length(x) x feature dimensionality] (e.g. 128 for
%   standard SIFT)

function [features] = get_features(image, x, y, descriptor_window_image_width)

% To start with, you might want to simply use normalized patches as your
% local feature. This is very simple to code and works OK. However, to get
% full credit you will need to implement the more effective SIFT descriptor
% (See Szeliski 4.1.2 or the original publications at
% http://www.cs.ubc.ca/~lowe/keypoints/)

% Your implementation does not need to exactly match the SIFT reference.
% Here are the key properties your (baseline) descriptor should have:
%  (1) a 4x4 grid of cells, each descriptor_window_image_width/4. 'cell' in this context
%    nothing to do with the Matlab data structue of cell(). It is simply
%    the terminology used in the feature literature to describe the spatial
%    bins where gradient distributions will be described.
%  (2) each cell should have a histogram of the local distribution of
%    gradients in 8 orientations. Appending these histograms together will
%    give you 4x4 x 8 = 128 dimensions.
%  (3) Each feature should be normalized to unit length
%
% You do not need to perform the interpolation in which each gradient
% measurement contributes to multiple orientation bins in multiple cells
% As described in Szeliski, a single gradient measurement creates a
% weighted contribution to the 4 nearest cells and the 2 nearest
% orientation bins within each cell, for 8 total contributions. This type
% of interpolation probably will help, though.

% You do not have to explicitly compute the gradient orientation at each
% pixel (although you are free to do so). You can instead filter with
% oriented filters (e.g. a filter that responds to edges with a specific
% orientation). All of your SIFT-like feature can be constructed entirely
% from filtering fairly quickly in this way.

% You do not need to do the normalize -> threshold -> normalize again
% operation as detailed in Szeliski and the SIFT paper. It can help, though.

% Another simple trick which can help is to raise each element of the final
% feature vector to some power that is less than one.

%Placeholder that you can delete. Empty features.
% features = zeros(size(x,1), 128, 'single');
%% Relatively Complicated Descriptor

% temp = x;
% x = y;
% y = temp;
% 
% image = im2double(image);
% if length(size(image)) == 3
%     image = rgb2gray(image);
% end
% sobelx = [-1, 0, 1; -2, 0, 2; -1, 0, 1];
% sobely = [1, 2, 1; 0, 0, 0; -1, -2, -1];
% 
% gaussian_sobelx = imgaussfilt(sobelx);
% gaussian_sobely = imgaussfilt(sobely);
% 
% I_x = imfilter(image, gaussian_sobelx, 'conv');
% I_y = imfilter(image, gaussian_sobely, 'conv');
% 
% % Alternative Method That Produces The Same Results
% % filterx = [0, 0, 0; -1, 0, 1; 0, 0, 0];
% % filtery = [0, 1, 0; 0, 0, 0; 0, -1, 0];
% % 
% % I_x = imfilter(image, filterx, 'conv');
% % I_y = imfilter(image, filtery, 'conv');
%
% grad_mag = (I_x.^2 + I_y.^2).^(0.5);
% grad_dir = atan2d(I_y,I_x);
% grad_dir(grad_dir<0) = grad_dir(grad_dir<0) + 360;
% 
% gauss_window = fspecial('gaussian', [descriptor_window_image_width,descriptor_window_image_width], 1.5);
% half_window = descriptor_window_image_width/2;
% quart_window = half_window/2;
% features = NaN(length(x),128,2);
%
% %%
%
% for i = 1:length(x)
%     subset = grad_dir(x(i)-half_window+1:x(i)+half_window,y(i)-half_window+1:y(i)+half_window);
%     mag_window = (grad_mag(x(i)-half_window+1:x(i)+half_window,y(i)-half_window+1:y(i)+half_window))/(sum(sum(grad_mag)));
%     subset_weights = gauss_window.*mag_window;
%     subset1 = ((floor(subset/10))*10)+5;
%     subset1(subset1==360) = 0;
%     
% % Alternative Method
% %     window = reshape(subset1,[1,descriptor_window_image_width*descriptor_window_image_width]);
% %     [counts,binLocations] = hist(window,5:10:355);
% 
%     counts = zeros(1,36);
%     binLocations = 5:10:355;
%     b = 1;
%     for bin = 1:36
%         counts(b) = sum(sum((subset1==binLocations(bin)).*subset_weights));
%         b = b + 1;
%     end
%     
%     sorted_counts = sort(counts);
%     dom_orien1 = binLocations(counts==sorted_counts(end));
%     subset11 = subset - dom_orien1(1);
% 
%     if (length(dom_orien1)>1)
%         dom_orien2 = dom_orien2(2);
%         subset12 = subset - dom_orien2(1);
%     elseif ((sorted_counts(end-1)/sorted_counts(end))>=0.8)
%         dom_orien2 = binLocations(counts==sorted_counts(end-1));
%         subset12 = subset - dom_orien2(1);
%     else
%         dom_orien2 = -1;
%         subset12 = zeros(descriptor_window_image_width,descriptor_window_image_width);
%     end
%     
%     subset21 = (floor(subset11/45))*45;
%     subset21(subset21==360) = 0;
%     subset22 = (floor(subset12/45))*45;
%     subset22(subset22==360) = 0;
%     
%     f1 = [];
%     for j = 1:quart_window:descriptor_window_image_width
%         for k = 1:quart_window:descriptor_window_image_width
%             new_window1 = reshape(subset21(j:j+quart_window-1,k:k+quart_window-1),[1,quart_window*quart_window]);
%             [counts2,binLocations2] = hist(new_window1,0:45:315);
%             f1 = [f1 counts2];
%         end
%     end
%     
%     f1 = f1/sum(f1);
%     features(i,:,1) = f1;
%     
%     if (dom_orien2~=-1)
%         f2 = [];
%         for j = 1:quart_window:descriptor_window_image_width
%             for k = 1:quart_window:descriptor_window_image_width
%                 new_window2 = reshape(subset22(j:j+quart_window-1,k:k+quart_window-1),[1,quart_window*quart_window]);
%                 [counts3,binLocations3] = hist(new_window2,0:45:315);
%                 f2 = [f2 counts3];
%             end
%         end
%         
%         f2 = f2/sum(f2);
%         features(i,:,2) = f2;
%     end
% end
%% Simplified Descriptor

temp = x;
x = y;
y = temp;

image = im2double(image);
if length(size(image)) == 3
    image = rgb2gray(image);
end
sobelx = [-1, 0, 1; -2, 0, 2; -1, 0, 1];
sobely = [1, 2, 1; 0, 0, 0; -1, -2, -1];

gaussian_sobelx = imgaussfilt(sobelx);
gaussian_sobely = imgaussfilt(sobely);

I_x = imfilter(image, gaussian_sobelx, 'conv');
I_y = imfilter(image, gaussian_sobely, 'conv');

grad_dir = atan2d(I_y,I_x);
grad_dir(grad_dir<0) = grad_dir(grad_dir<0) + 360;

half_window = descriptor_window_image_width/2;
quart_window = half_window/2;
features = NaN(length(x),128,2);

%%

for i = 1:length(x)
    subset = grad_dir(x(i)-half_window+1:x(i)+half_window,y(i)-half_window+1:y(i)+half_window);
 
    subset = (floor(subset/45))*45;
    subset(subset==360) = 0;
    
    f1 = [];
    for j = 1:quart_window:descriptor_window_image_width
        for k = 1:quart_window:descriptor_window_image_width
            new_window1 = reshape(subset(j:j+quart_window-1,k:k+quart_window-1),[1,quart_window*quart_window]);
            [counts2,binLocations2] = hist(new_window1,0:45:315);
            f1 = [f1 counts2];
        end
    end
    
    f1 = f1/sum(f1);
    features(i,:,1) = f1;
end
end









