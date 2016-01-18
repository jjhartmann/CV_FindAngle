%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Find angles of pencils using Hough Functions

% Load an image.
img = imread('images/SixCrossed.JPG');
imshow(img)

% Gray Scale
img_g = rgb2gray(img);



th = 120;
th_image = zeros(size(img_g)); % create a black image with size of bw
th_index = find(img_g<th); % Finds all the index that satisfy the threshold
th_image(th_index) = 1; % Takes all the indices and assigns one to them, creating a black and white file from the gray scale. 
imshow(th_image)

BW = edge(th_image, 'canny');
imshow(BW)

% Compute hough transform. 
[H, theta, rho] = hough(BW);

% Display transform data
figure, imshow(imadjust(mat2gray(H)), [], 'XData', theta, 'YData', rho,'InitialMagnification', 'fit');
xlabel('\theta (degrees)'), ylabel('\rho');
axis on, axis normal, hold on;
colormap(hot)

% Find peaks in teh Hough transform
peaks = houghpeaks(H, 3, 'threshold', ceil(0.3*max(H(:))));

% Plot on colormap
x = theta(peaks(:, 2));
y = rho(peaks(:,1));
plot(x,y,'s', 'color', 'black')
 
% Find lines using houghlines
lines = houghlines(BW, theta, rho, peaks, 'FillGap', 5, 'MinLength', 7);

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
    
    
    
    
    



