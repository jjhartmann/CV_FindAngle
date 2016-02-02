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

%% Determine if the image is light or dark
GrayImageM = rgb2gray(img);
imgMedian = median(GrayImageM(:));

if (imgMedian < 150)
%% Define thresholds for channel 1 based on histogram settings
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Convert RGB image to LAB
    RGB = im2double(img);
    cform = makecform('srgb2lab', 'AdaptedWhitePoint', whitepoint('D65'));
    I = applycform(RGB,cform);

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
    RGB = im2double(img);
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

%% Convert to gray scale mask
maskedGray = rgb2gray(maskedRGBImage);
imshow(maskedGray)

%% Find the edges of the image using laplacian of gaussian. 
imgedg = edge(maskedGray, 'log');
imshow(imgedg)

smooth = imdilate(imgedg,strel('disk',1));
imshow(smooth)

%% Find Regions and display
[B, L] = bwboundaries(smooth, 'holes');

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

%% Crop the outer ridges of the image. 
[width, height] = size(L_thin);
L_thin = imcrop(L_thin,[30,30,height-60, width-60]);
img = imcrop(img,[30,30,height-60, width-60]);
imshow(L_thin)

%smooth image
L_smooth = imdilate(L_thin,strel('disk',1));
imshow(L_smooth)

%% Compute hough transform. 
[H, theta, rho] = hough(L_smooth);

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
lines = houghlines(L_smooth, theta, rho, peaks, 'FillGap', 30, 'MinLength', 100);

%% Eliminate lines that are simlar. ie some epsilon difference between them.
for k = 1:length(lines)
    deltaPoint = lines(k).point2 - lines(k).point1;
    lines(k).Length = sqrt(deltaPoint(1)^2 + deltaPoint(2)^2);
    
    % Create angles in degrees of lines
    lines(k).LineDegree = (atan(double(deltaPoint(2)/deltaPoint(1))) * double(180/pi)) + 90;
    
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
    
    % find position
    Lmidpoint = [(xy(1,1) + xy(2,1))/2, (xy(1,2) + xy(2,2))/2];
    plot(Lmidpoint(1), Lmidpoint(2), 'o', 'LineWidth', 2, 'Color', 'blue');
    
    % Print output
    disp(['Pencil ', num2str(k) ':'])
    disp(['  Length:          ', num2str(lines(k).Length)])
    disp(['  Angle:           ', num2str(atan((lines(k).point2(1) - Lmidpoint(1))/(lines(k).point2(2) - Lmidpoint(2))) * (180/pi) + 90)])
    disp(['  Center Position: [', num2str(Lmidpoint(1)), ', ', num2str(Lmidpoint(2)), ']'])
    disp(['  End Position:    [', num2str(lines(k).point1), ']'])
    disp(['  Begin Position:  [',  num2str(lines(k).point2), ']'])

    disp(' ')
    
end

%% hightlight the n longest lines
plot(xy_long(:, 1), xy_long(:,2), 'LineWidth', 2, 'Color', 'red');
hold off;



    
    
    