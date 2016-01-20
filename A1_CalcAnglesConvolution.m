%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calc the angles, orientation, and position of pencils. 

% Get the image. 
img = imread('images/CrossedPencilsA.JPG');
imshow(img)

%convert the image to a gray scale image
bw = double(double(img(:,:,1)) + double(img(:, :, 2)) + double(img(:, :, 3)))/(3*255);
imshow(bw);

