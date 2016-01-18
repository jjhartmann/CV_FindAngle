%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Find all angles in a image. 

% Load Image. 
image = imread('images/CrossedPencilsA.JPG');
imshow(image)

% Create a gray scale image
img_gray = rgb2gray(image);
imshow(img_gray)

% Get Threshold and convert to Black and white
threshold = 65/100;
BW = im2bw(img_gray, threshold);
BW = BW; % Complement
imshow(BW)

dim = size(BW);

% Determine penicl
row = 200;
col = find(BW(row,:), 1);

boundary = bwtraceboundary(BW, [row, col], 'S', 8, 100, 'counter');

imshow(image); hold on;
plot(boundary(:, 2), boundary(:, 1), 'g', 'LineWidth', 2);
hold off;

% Determine the angle. 
line = polyfit(boundary(:, 2), boundary(:, 1), 1);
cx = 0;
cy = 0;
ex = 50;
ey = line(1) * 50;

dx = ex - cx;
dy = ey - cy;
theta = atan(dy/dx);
angle = theta * 180/pi;

disp(['Angle: ', num2str(abs(angle))])

