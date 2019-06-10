% Local Feature Stencil Code
% Written by James Hays for CS 143 @ Brown / CS 4476/6476 @ Georgia Tech

% Returns a set of interest points for the input image

% 'image' can be grayscale or color, your choice.
% 'descriptor_window_image_width', in pixels.
%   This is the local feature descriptor width. It might be useful in this function to 
%   (a) suppress boundary interest points (where a feature wouldn't fit entirely in the image, anyway), or
%   (b) scale the image filters being used. 
% Or you can ignore it.

% 'x' and 'y' are nx1 vectors of x and y coordinates of interest points.
% 'confidence' is an nx1 vector indicating the strength of the interest
%   point. You might use this later or not.
% 'scale' and 'orientation' are nx1 vectors indicating the scale and
%   orientation of each interest point. These are OPTIONAL. By default you
%   do not need to make scale and orientation invariant local features.
function [x, y, confidence, scale, orientation] = get_interest_points(image, descriptor_window_image_width)

% Implement the Harris corner detector (See Szeliski 4.1.1) to start with.
% You can create additional interest point detector functions (e.g. MSER)
% for extra credit.

% If you're finding spurious interest point detections near the boundaries,
% it is safe to simply suppress the gradients / corners near the edges of
% the image.

% The lecture slides and textbook are a bit vague on how to do the
% non-maximum suppression once you've thresholded the cornerness score.
% You are free to experiment. Here are some helpful functions:
%  BWLABEL and the newer BWCONNCOMP will find connected components in 
% thresholded binary image. You could, for instance, take the maximum value
% within each component.
%  COLFILT can be used to run a max() operator on each sliding window. You
% could use this to ensure that every interest point is at a local maximum
% of cornerness.

% Placeholder that you can delete -- random points
% x = ceil(rand(500,1) * size(image,2));
% y = ceil(rand(500,1) * size(image,1));

%% Harris Corner Detector

image = im2double(image);
if length(size(image)) == 3
    image = rgb2gray(image);
end

[im_rows, im_columns] = size(image);
half_window = descriptor_window_image_width/2;

sobelx = [-1, 0, 1; -2, 0, 2; -1, 0, 1];
sobely = [1, 2, 1; 0, 0, 0; -1, -2, -1];

gaussian_sobelx = imgaussfilt(sobelx);
gaussian_sobely = imgaussfilt(sobely);

I_x = imfilter(image, gaussian_sobelx, 'conv');
I_y = imfilter(image, gaussian_sobely, 'conv');
I_xx = I_x.*I_x;
I_yy = I_y.*I_y;
I_xy = I_x.*I_y;

S_xx = imgaussfilt(I_xx);
S_yy = imgaussfilt(I_yy);
S_xy = imgaussfilt(I_xy);
%%
% hyper-parameters
alpha = 0.04;
threshold = 0.01;

R = ((S_xx.*S_yy) - (S_xy.*S_xy)) - (alpha*((S_xx+S_yy).^2));
% Alternative Method
% R = ((S_xx.*S_yy) - (S_xy.*S_xy))./(S_xx+S_yy);

R(R<threshold) = 0;
R_padded = padarray(R,[1 1],'both');
R_new_padded = zeros(size(R_padded));

for i = 2+half_window:im_rows+1-half_window
    for j = 2+half_window:im_columns+1-half_window
        subset = R_padded(i-1:i+1,j-1:j+1);
        max_subset = max(max(subset));
        if (R_padded(i,j)==max_subset)
            R_new_padded(i,j) = R_padded(i,j);
        end
    end
end

%%
% Uncoment to view
% imshow(image)
% hold on

x = [];
y = [];
for i = 2+half_window:im_rows+1-half_window
    for j = 2+half_window:im_columns+1-half_window
        if(R_new_padded(i,j)~=0)
            x = [x; j-1];
            y = [y; i-1];
        end       
    end
end

% Uncoment to view
% scatter(x,y)

%% SIFT Detector

% image = im2double(image);
% if length(size(image)) == 3
%     image = rgb2gray(image);
% end
% 
% image_blur1 = image;
% im_size = size(image_blur1);
% [im_rows, im_columns] = size(image);
% 
% image_blur2 = imgaussfilt(image_blur1,1);
% image_blur3 = imgaussfilt(image_blur1,2);
% image_blur4 = imgaussfilt(image_blur1,3);
% image_blur5 = imgaussfilt(image_blur1,4);
% 
% LoG1 = abs(image_blur1 - image_blur2);
% LoG2 = abs(image_blur2 - image_blur3);
% LoG3 = abs(image_blur3 - image_blur4);
% LoG4 = abs(image_blur4 - image_blur5);
% 
% LoG1_padded = padarray(LoG1,[1 1],'both');
% LoG2_padded = padarray(LoG2,[1 1],'both');
% LoG3_padded = padarray(LoG3,[1 1],'both');
% LoG4_padded = padarray(LoG4,[1 1],'both');
% 
% step31 = zeros(im_size + 2);
% step32 = zeros(im_size + 2);
% 
% for i = 2:im_size(1)-1
%     for j = 2:im_size(2)-1
%         subset1 = LoG1_padded(i-1:i+1,j-1:j+1);
%         max1 = max(max(subset1));
%         min1 = min(min(subset1));
%         subset2 = LoG2_padded(i-1:i+1,j-1:j+1);
%         max2 = max(max(subset2));
%         min2 = min(min(subset2));
%         subset3 = LoG3_padded(i-1:i+1,j-1:j+1);
%         max3 = max(max(subset3));
%         min3 = min(min(subset3));
%         max_LoG2 = max(([max1,max2,max3]));
%         min_LoG2 = min(([min1,min2,min3]));
%         if ((LoG2_padded(i,j)==max_LoG2)||(LoG2_padded(i,j)==min_LoG2))
%             step31(i,j) = LoG2_padded(i,j);
%         end
%         
%         subset4 = LoG4_padded(i-1:i+1,j-1:j+1);
%         max4 = max(max(subset4));
%         min4 = min(min(subset4));
%         max_LoG3 = max(([max2,max3,max4]));
%         min_LoG3 = min(([min2,min3,min4]));
%         if ((LoG3_padded(i,j)==max_LoG3)||(LoG3_padded(i,j)==min_LoG3))
%             step32(i,j) = LoG3_padded(i,j);
%         end
%     end
% end
% 
% threshold = 0.02;
% 
% step31(LoG2_padded<threshold) = 0;
% step32(LoG3_padded<threshold) = 0;

% %%

% % Uncoment to view
% % imshow(image_blur1)
% % hold on
% 
% x = [];
% y = [];
% half_window = descriptor_window_image_width/2;
% 
% for i = 2+half_window:im_rows+1-half_window
%     for j = 2+half_window:im_columns+1-half_window
%         if(step31(i,j)~=0)
%             x = [x; j-1];
%             y = [y; i-1];
%         end
%         
%         if(step32(i,j)~=0)
%             x = [x; j-1];
%             y = [y; i-1];
%         end
%        
%     end
% end
% 
% % Uncoment to view
% % scatter(x,y)
end

