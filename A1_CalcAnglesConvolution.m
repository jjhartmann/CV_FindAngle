%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calc the angles, orientation, and position of pencils. 

%% Get the image.
sel = input('Choose test method: \n1: OnePencilA.JPG\n2: Red_Green_Pencils.JPG\n3: CrossedPencilsA.JPG\n4: Three-on-capret.JPG\n5: Touching.JPG\n6: SixCorssed.JPG\n\n>> ');

switch sel
    case 1 
        image_text = 'images/OnePencilA.JPG';
    case 2  
        image_text = 'images/Red_Green_Pencils.JPG';
    case 3  
        image_text = 'images/CrossedPencilsA.JPG';
    case 4  
        image_text = 'images/Three-on-Carpet.JPG';
    case 5 
        image_text = 'images/Touching.JPG';
    case 6
        image_text = 'images/SixCrossed.JPG';
    otherwise print('Not Valid input'), quit();
end

img = imread(image_text);
imshow(img)

%% Convert RGB image to HVS
RGB = im2double(img);
cform = makecform('srgb2lab', 'AdaptedWhitePoint', whitepoint('D65'));
I = applycform(RGB,cform);

%% Determine if the image is light or dark
GrayImageM = rgb2gray(RGB);
imgMedian = median(GrayImageM(:));

if (imgMedian < 0.55)
%% Define thresholds for channel 1 based on histogram settings
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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

else
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Apply gaussian smoothing
    img = imgaussfilt3(img, 2);
    imshow(img)

    % Convert RGB image to chosen color space
    I = img;

    % Define thresholds for channel 1 based on histogram settings
    channel1Min = 170.000;
    channel1Max = 255.000;

    % Define thresholds for channel 2 based on histogram settings
    channel2Min = 167.000;
    channel2Max = 255.000;

    % Define thresholds for channel 3 based on histogram settings
    channel3Min = 162.000;
    channel3Max = 255.000;

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
    imshow(maskedRGBImage)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end


maskedGray = rgb2gray(maskedRGBImage);
imshow(maskedGray)

%% %convert the image to a gray scale image
% bw = double(1.2 * double(img(:,:,1)) + 1.2 * double(img(:, :, 2)) + .2 * double(img(:, :, 3)))/(3*255);
% imshow(bw);

% Convert to binary image
% Threshold
% th = 100;
% th_image = im2bw(maskedGray, th/255);
% imshow(th_image)
% 
% % Clean up image based on segments
% th_image = bwareaopen(th_image, 100);
% imshow(th_image)
% 
% % Clean the border of the image from artifacts. 
% th_image = imclearborder(th_image);
% imshow(th_image)

% % TODO: Smoothing of the image. Needs to promote the pencils more. 
% % laplacian of gaussian
% h = fspecial('log',15,3);
% bg=filter2(h, th_image);
% surf(h);  %Display the mask
% %Can't display the result of the convolution directly because it has
% %negative values.
% %Add an offset and then scale to 0-1 range for display purposes
% bgs=bg+.5;
% bgs=bgs-min(bgs(:));
% bgs=bgs/max(bgs(:));
% imshow(bgs);

%% Find the edges of the image using laplacian of gaussian. 
imgedg = edge(maskedGray, 'log');
imshow(imgedg)

maskedGray(find(maskedGray)) = 1;

smooth = imdilate(imgedg,strel('disk',1));
imshow(smooth)

%% Find Regions and display
[B, L] = bwboundaries(smooth, 'holes');

%% TEST
% L = maskedGray;
%%

imshow(label2rgb(L))

%% Clean the border of the image from artifacts. 
L = imclearborder(L);
imshow(L)

numRegions = max(L(:));

%% Clean up image based on segments
L = bwareaopen(L, 400);
imshow(L)

%% Thin the image
L_thin = bwmorph(L, 'thin', inf);
imshow(L_thin)
imgedg = L_thin;

%% % Find shapes based on eccentricity
% stats = regionprops(L, 'all');
% shapes = [stats.Eccentricity];
% pencils_index = find(shapes > 0.98);
% 
% 
% % Location of pencils
% pencils = stats(pencils_index);
% locations = [pencils.Centroid];
% 
% 
% figure;
% imshow(th_image); hold on;
% for k = 1:(length(locations)/2)
%     plot(locations((k*2)-1), locations((k*2)), 'x', 'LineWidth', 2, 'Color', 'red'); 
% end


%% Crop the outer ridges of the image. 
[width, height] = size(imgedg);
imgedg = imcrop(imgedg,[30,30,height-60, width-60]);
img = imcrop(img,[30,30,height-60, width-60]);
imshow(imgedg)

%smooth image
imgedg = imdilate(imgedg,strel('disk',1));
imshow(imgedg)

%% Compute hough transform. 
[H, theta, rho] = hough(imgedg);

% Display transform data
figure, imshow(imadjust(mat2gray(H)), [], 'XData', theta, 'YData', rho,'InitialMagnification', 'fit');
xlabel('\theta (degrees)'), ylabel('\rho');
axis on, axis normal, hold on;
colormap(hot)

%% Find peaks in teh Hough transform
% change peaks based on blobs. 
peaks = houghpeaks(H, 6, 'threshold', 10);

%% Plot on colormap
x = theta(peaks(:, 2));
y = rho(peaks(:,1));
plot(x,y,'s', 'color', 'black')
 
%% Find lines using houghlines
lines = houghlines(imgedg, theta, rho, peaks, 'FillGap', 30, 'MinLength', 100);

%% Eliminate lines that are simlar. ie some epsilon difference between them.
for k = 1:length(lines)
    deltaPoint = lines(k).point2 - lines(k).point1;
    lines(k).Length = sqrt(deltaPoint(1)^2 + deltaPoint(2)^2);
end

% Sort by length
linesField = fieldnames(lines);
linesCell = struct2cell(lines);
linesSZ = size(linesCell);

linesCell = reshape(linesCell, linesSZ(1), []);
linesCell = linesCell';
[trash, idx] = sort(cell2mat(linesCell(:,5)) , 'descend');

% Put back in struct
lines = lines(idx);

%% Remove duplicates. 

lineLen = length(lines);
lenEpsilon = 40;
xyEpsilon = 30;

k = 1;
while k <= lineLen
    
    xy1a = lines(k).point1;
    xy2a = lines(k).point2;
    % Iterate over the next K+1 lines;
    i = k + 1;
    while i <= lineLen
        
        xy1b = lines(i).point1;
        xy2b = lines(i).point2;
       
        % compare end points. 
        deltaLength = abs(lines(k).Length - lines(i).Length); 
        deltaXY1 = (abs(xy1a(1) - xy1b(1)) + abs(xy1a(2) - xy1b(2)));
        deltaXY2 = (abs(xy2a(1) - xy2b(1)) + abs(xy2a(2) - xy2b(2)));

        % remove the line
        if (deltaXY1 <= xyEpsilon || deltaXY2 <= xyEpsilon)
            % check to make sure it fits on line.
            f1 = polyfit([xy1a(1), xy2a(1)], [xy1a(2), xy2a(2)], 1);
            f2 = polyfit([xy1b(1), xy2b(1)], [xy1b(2), xy2b(2)], 1);
            
            testX = 0;
            if (deltaXY1 <= xyEpsilon)
                testX = lines(i).point2(1);
            else
                testX = lines(i).point1(1);
            end
            
            pY1 = polyval(f1, testX);
            pY2 = polyval(f2, testX);
            
            deltaP = abs(pY1 - pY2);
            
            % remove line if a fit
            if (deltaP <= xyEpsilon)
                lines(:,i)=[];
                i = i - 1;
                lineLen = lineLen - 1;
            end
        end  

        i = i + 1;
    end
    
    k = k + 1;
end

%% Find orientation



%% Plot lines on original image. 
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
    
    % Print output
    disp(['Pencil ', num2str(k) ':'])
    disp(['Length: ', num2str(lines(k).Length)])
    disp(['Angle: ', num2str(lines(k).theta)])
    
    % find position
    Lmidpoint = [(xy(1,1) + xy(2,1))/2, (xy(1,2) + xy(2,2))/2];
    plot(Lmidpoint(1), Lmidpoint(2), 'o', 'LineWidth', 2, 'Color', 'blue');
   
    disp(['Position: [', num2str(Lmidpoint(1)), ', ', num2str(Lmidpoint(2)), ']'])
    disp(' ')
    
end

%% hightlight the n longest lines
plot(xy_long(:, 1), xy_long(:,2), 'LineWidth', 2, 'Color', 'red');
hold off;



    
    
    