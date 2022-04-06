% Master Thesis: Real-Time Stereo Vision     Wim Abbeloos    May 2010
% Karel de Grote-Hogeschool University College, Belgium
%
% FAST MATLAB STEREO MATCHING ALGORITHM (SAD, Sum of absolute intensity differences)
% Description: This function performs the computationally expensive step of
% matching two rectified and undistorted stereo images.  The output is a 
% dense disparity map.  If camera parameters are known, this allows for 
% three dimensional reconstruction.
%
% Please note this function requires the Image Processing Toolbox!
%
% [depthmap, dcost, pcost, wcost] = stereomatch(imgleft, imgright, windowsize, disparity, smoothing)
%
% The standard images included are from
% [1] 	D. Scharstein and R. Szeliski. A taxonomy and evaluation of dense two-frame stereo correspondence algorithms.
% International Journal of Computer Vision, 47(1/2/3):7-42, April-June 2002.
% [2] 	D. Scharstein and R. Szeliski. High-accuracy stereo depth maps using structured light.
% In IEEE Computer Society Conference on Computer Vision and Pattern Recognition (CVPR 2003), volume 1, pages 195-202, Madison, WI, June 2003. 

function [depthmap, dcost, pcost, wcost] = stereomatch_s(imgleft, imgright, windowsize, disparity, smoothing)

% Set Parameters
WS  = uint16(windowsize);               % Set window size, must be uneven
WS2 = uint16( ( WS - 1 ) / 2 );         % Half window
D   = uint16(disparity)+1;              % number of disparities

% Read image sizes
heightL = uint16( size( imgleft, 1 ) );    heightR = uint16( size( imgright, 1 ) );
widthL  = uint16( size( imgleft, 2 ) );    widthR  = uint16( size( imgright, 2 ) );
if ( heightL ~= heightR  ||  widthL ~= widthR )
    error('Height and width of left and right image must be equal');
end

% Initialization
hd = waitbar(0, 'please wait ....'); %new
pcost = zeros( heightL, widthL, D, 'uint8' );
wcost = zeros( heightL, widthL, D, 'single' );
dmap  = zeros( heightL, widthL, 'uint8' );
dcost = zeros( heightL, widthL, 'single' );
h = zeros(WS,WS,'double'); h(1,1) = 1; h(1,WS) = -1; h(WS,1) = -1; h(WS,WS) = 1;

% Calculate pixel cost
imgleft = single(imgleft); %new
imgright = single(imgright); %new
for Dc = 1 : D
   maxL = widthL + 1 - Dc; 
   %pcost(:, Dc : widthL, Dc ) = imabsdiff( imgright( :, 1 : maxL), imgleft( :, Dc : widthL) );
    absDiff = abs( imgright( :, 1 : maxL) - imgleft( :, Dc : widthL) ); %new
    pcost(:, Dc : widthL, Dc ) = mean(absDiff, 3); %new
end

% Calculate integral cost
icost = single(pcost);
icost = cumsum( cumsum( icost ), 2 );

% Calculate window cost
%wcost = imfilter(icost,h,'same','symmetric');
icost_expend = cat(2, repmat(icost(:,1,:), 1, WS2), icost, repmat(icost(:,end,:), 1, WS2)); %new
icost_expend = cat(1, repmat(icost_expend(1,:,:), WS2,1), icost_expend, repmat(icost_expend(end,:,:), WS2,1)); %new
for Dc = 1 : D %new
    waitbar(double(Dc)/double(D), hd) %new
    temp = filter2(single(h), icost_expend(:, :, Dc)); %new
    wcost(:,:,Dc) = temp(WS2+1:end-WS2, WS2+1:end-WS2); %new
end %new

% Search disparity value
[ dcost(:,D+WS2:widthL), dmap(:,D+WS2:widthL)] = min( wcost(:,D+WS2:widthL,:),[], 3 );
for j=WS2+1:D+WS2
    [ dcost(:,j), dmap(:,j)] = min( wcost(:, j, 1 : (j - WS2) ),[], 3 );
end

% Adjust disparity map
warning off;
depthmap = double(dmap-1);

wx = double(widthL)/10;
for y=1:wx
    depthmap(:,y) = depthmap(:,y) + mean(depthmap(:))*(wx - y)/(wx);
end

x= sort(depthmap(:)); % sort the depth map, new
x_1_percentile= x(round(0.01*numel(x)));  % calculate 1st percentile, new
depthmap(find(depthmap < x_1_percentile)) = x_1_percentile; %new
x_99_percentile= x(round(0.99*numel(x)));  % calculate 99th percentile, new
depthmap(find(depthmap > x_99_percentile)) = x_99_percentile; %new

if smoothing ==1 %depth map blurring  %new
    depthmap = filter2(ones(7)/49, depthmap);  %new
end  %new

delete(hd) %new
warning on;

