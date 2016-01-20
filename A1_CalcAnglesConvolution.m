%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calc the angles, orientation, and position of pencils. 

% Get the image. 
img = imread('images/CrossedPencilsA.JPG');
imshow(img)

%convert the image to a gray scale image
bw = double(double(img(:,:,1)) + double(img(:, :, 2)) + double(img(:, :, 3)))/(3*255);
imshow(bw);

% TODO: Smoothing of the image. Needs to promote the pencils more. 
% convert image to a threshold black and white image
th = 110/255;
th_image = -ones(size(bw)); % create a black image with size of bw
th_index = find(bw<th); % Finds all the index that satisfy the threshold
th_image(th_index) = 1; % Takes all the indices and assigns one to them, creating a black and white file from the gray scale. 

imshow(th_image)
