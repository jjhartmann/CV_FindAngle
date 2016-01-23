%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calc the angles, orientation, and position of pencils. 

% Get the image. 
img = imread('images/SixCrossed.JPG');
imshow(img)

%convert the image to a gray scale image
bw = double(double(img(:,:,1)) + double(img(:, :, 2)) + double(img(:, :, 3)))/(3*255);
imshow(bw);

% TODO: Smoothing of the image. Needs to promote the pencils more. 
% laplacian of gaussian
h = fspecial('log',26,7);
bg=filter2(h, bw);
surf(h);  %Display the mask
%Can't display the result of the convolution directly because it has
%negative values.
%Add an offset and then scale to 0-1 range for display purposes
bgs=bg+.5;
bgs=bgs-min(bgs(:));
bgs=bgs/max(bgs(:));
imshow(bgs);

% Threshold
th = .6;
th_image = zeros(size(bgs)); % create a black image with size of bw
th_index = find(bgs<th); % Finds all the index that satisfy the threshold
th_image(th_index) = 1; % Takes all the indices and assigns one to them, creating a black and white file from the gray scale. 
imshow(th_image)

imgedg = edge(bgs, 'canny');
imshow(imgedg)
