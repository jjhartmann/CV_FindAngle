%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calc the angles, orientation, and position of pencils. 

% Get the image. 
img = imread('images/Red_Green_Pencils.JPG');
imshow(img)

% Convert RGB image to chosen color space
RGB = im2double(img);
cform = makecform('srgb2lab', 'AdaptedWhitePoint', whitepoint('D65'));
I = applycform(RGB,cform);
imshow(RGB)

% Define thresholds for channel 1 based on histogram settings
channel1Min = 13.262;
channel1Max = 67.444;

% Define thresholds for channel 2 based on histogram settings
channel2Min = -1.205;
channel2Max = 14.341;

% Define thresholds for channel 3 based on histogram settings
channel3Min = -13.473;
channel3Max = 13.090;

% Create mask based on chosen histogram thresholds
BW = (I(:,:,1) >= channel1Min ) & (I(:,:,1) <= channel1Max) & ...
    (I(:,:,2) >= channel2Min ) & (I(:,:,2) <= channel2Max) & ...
    (I(:,:,3) >= channel3Min ) & (I(:,:,3) <= channel3Max);

% Invert mask
BW = ~BW;

% Initialize output masked image based on input image.
maskedRGBImage = RGB;

% Set background pixels where BW is false to zero.
maskedRGBImage(repmat(~BW,[1 1 3])) = 0;
maskedGray = rgb2gray(maskedRGBImage);
imshow(maskedGray)

% %convert the image to a gray scale image
% bw = double(1.2 * double(img(:,:,1)) + 1.2 * double(img(:, :, 2)) + .2 * double(img(:, :, 3)))/(3*255);
% imshow(bw);

% Convert to binary image
% Threshold
th = 100;
th_image = im2bw(maskedGray, th/255);
imshow(th_image)

% Clean up image based on segments
th_image = bwareaopen(th_image, 100);
imshow(th_image)

% Clean the border of the image from artifacts. 
th_image = imclearborder(th_image);
imshow(th_image)

% TODO: Smoothing of the image. Needs to promote the pencils more. 
% laplacian of gaussian
h = fspecial('log',15,3);
bg=filter2(h, th_image);
surf(h);  %Display the mask
%Can't display the result of the convolution directly because it has
%negative values.
%Add an offset and then scale to 0-1 range for display purposes
bgs=bg+.5;
bgs=bgs-min(bgs(:));
bgs=bgs/max(bgs(:));
imshow(bgs);

%Find the edges of the image using laplacian of gaussian. 
imgedg = edge(bgs, 'log');
imshow(imgedg)

%Crop the outer ridges of the image. 
[width, height] = size(imgedg);
imgedg = imcrop(imgedg,[30,30,height-60, width-60]);
img = imcrop(img,[30,30,height-60, width-60]);
imshow(imgedg)


% Compute hough transform. 
[H, theta, rho] = hough(imgedg);

% Display transform data
figure, imshow(imadjust(mat2gray(H)), [], 'XData', theta, 'YData', rho,'InitialMagnification', 'fit');
xlabel('\theta (degrees)'), ylabel('\rho');
axis on, axis normal, hold on;
colormap(hot)

% Find peaks in teh Hough transform
peaks = houghpeaks(H, 3, 'threshold', ceil(0.7*max(H(:))));

% Plot on colormap
x = theta(peaks(:, 2));
y = rho(peaks(:,1));
plot(x,y,'s', 'color', 'black')
 
% Find lines using houghlines
lines = houghlines(imgedg, theta, rho, peaks, 'FillGap', 5, 'MinLength', 7);

% Plot lines on original image. 
figure, imshow(img), hold on
max_len = 0 % Todo: expand for the top n lines
xy_long = 0;

for k = 1:length(lines)
    xy = [lines(k).point1; lines(k).point2];
    plot(xy(:, 1), xy(:,2), 'Linewidth', 2, 'Color', 'green');
    
    % Plot end points
    plot(xy(1,1), xy(1,2), 'x', 'LineWidth', 2, 'Color', 'yellow');
    plot(xy(2,1), xy(2,2), 'x', 'LineWidth', 2, 'Color', 'red');
    
    % Find the longest lines
    % TODO: find n longest lines
    len = norm(lines(k).point1 - lines(k).point2);
    if (len > max_len)
        max_len = len;
        xy_long = xy;
    end
    
end

% hightlight the n longest lines
plot(xy_long(:, 1), xy_long(:,2), 'LineWidth', 2, 'Color', 'red');
    
    
    