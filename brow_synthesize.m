function [ZleftBrow bigLeftBrowRect ZrightBrow bigRightBrowRect] = brow_synthesize(X, leftBrow, rightBrow)

addpath(genpath(pwd));

bin = 8;
Angle = 360;
L = 3;

neighborsCount = 30;
param.lambda = 0.001;

% Synthesize left brow
load './leftBrow.mat';

% Middle point
leftMiddlePoint = [round((min(leftBrow(1,:))+max(leftBrow(1,:)))/2) round((max(leftBrow(2,:))+min(leftBrow(2,:)))/2)];

% Left -> right
bigLeftBrowRect(1,1) = leftMiddlePoint(1,1) - round(globalLeftBrowRectWidth/2);
bigLeftBrowRect(1,2) = leftMiddlePoint(1,1) + round(globalLeftBrowRectWidth/2);
% Top -> bottom
bigLeftBrowRect(1,3) = leftMiddlePoint(1,2) - round(globalLeftBrowRectHeight/2);
bigLeftBrowRect(1,4) = leftMiddlePoint(1,2) + round(globalLeftBrowRectHeight/2);

XleftBrow = X(bigLeftBrowRect(1,3):bigLeftBrowRect(1,4),bigLeftBrowRect(1,1):bigLeftBrowRect(1,2),1);
leftBrowROI(1,1) = bigLeftBrowRect(1,3);
leftBrowROI(2,1) = bigLeftBrowRect(1,4);
leftBrowROI(3,1) = bigLeftBrowRect(1,1);
leftBrowROI(4,1) = bigLeftBrowRect(1,2);
leftBrowPhog = anna_phog(X(:,:,1),bin,Angle,L,leftBrowROI);

leftBrowIndex = flann_load_index('./leftBrow_index.index', leftBrowPhogs);
[result, dists] = flann_search(leftBrowIndex, leftBrowPhog, neighborsCount, leftBrowParameters);
flann_free_index(leftBrowIndex);
for i=1:neighborsCount
    d = leftBrowImages{result(i,1)};
    imwrite(d, ['./patch/leftBrow' num2str(uint8(i)) '.jpg']);
    d = rgb2ycbcr(uint8(d));
    d = d(:,:,1);
    leftBrowDic(:,i) = double(d(:));
end

alpha = mexLasso(double(XleftBrow(:)),double(leftBrowDic),param);
ZleftBrow = reshape(leftBrowDic*alpha,size(XleftBrow));

%figure,imshow(XleftBrow);
%figure,imshow(uint8(ZleftBrow));

XleftBrow(:,:,2:3) = X(bigLeftBrowRect(1,3):bigLeftBrowRect(1,4),bigLeftBrowRect(1,1):bigLeftBrowRect(1,2),2:3);
%figure,imshow(ycbcr2rgb(XleftBrow));
ZleftBrow(:,:,2:3) = XleftBrow(:,:,2:3);
%figure,imshow(ycbcr2rgb(uint8(ZleftBrow)));

%{
im = zeros(size(XleftBrow,1), size(XleftBrow,2)*2+5, 3);
im(:,1:size(XleftBrow,2),:) = ycbcr2rgb(XleftBrow);
im(:,size(XleftBrow,2)+6:size(XleftBrow,2)*2+5,:) = ycbcr2rgb(uint8(ZleftBrow));
figure,imshow(uint8(im));
%}

%ZleftBrow(:,:,1) = global_color_correction(XleftBrow(:,:,1), ZleftBrow(:,:,1));

imwrite(ycbcr2rgb(XleftBrow), './patch/XleftBrow.jpg');
imwrite(ycbcr2rgb(uint8(ZleftBrow)), './patch/ZleftBrow.jpg');

im = zeros(size(XleftBrow,1), size(XleftBrow,2)*2+5, 3);
im(:,1:size(XleftBrow,2),:) = ycbcr2rgb(XleftBrow);
im(:,size(XleftBrow,2)+6:size(XleftBrow,2)*2+5,:) = ycbcr2rgb(uint8(ZleftBrow));
figure,imshow(uint8(im));

% Synthesize right brow
load './rightBrow.mat';

% Middle point
rightMiddlePoint = [round((min(rightBrow(1,:))+max(rightBrow(1,:)))/2) round((max(rightBrow(2,:))+min(rightBrow(2,:)))/2)];

% Left -> right
bigRightBrowRect(1,1) = rightMiddlePoint(1,1) - round(globalRightBrowRectWidth/2);
bigRightBrowRect(1,2) = rightMiddlePoint(1,1) + round(globalRightBrowRectWidth/2);
% Top -> bottom
bigRightBrowRect(1,3) = rightMiddlePoint(1,2) - round(globalRightBrowRectHeight/2);
bigRightBrowRect(1,4) = rightMiddlePoint(1,2) + round(globalRightBrowRectHeight/2);

XrightBrow = X(bigRightBrowRect(1,3):bigRightBrowRect(1,4),bigRightBrowRect(1,1):bigRightBrowRect(1,2),1);
rightBrowROI(1,1) = bigRightBrowRect(1,3);
rightBrowROI(2,1) = bigRightBrowRect(1,4);
rightBrowROI(3,1) = bigRightBrowRect(1,1);
rightBrowROI(4,1) = bigRightBrowRect(1,2);
rightBrowPhog = anna_phog(X(:,:,1),bin,Angle,L,rightBrowROI);

rightBrowIndex = flann_load_index('./rightBrow_index.index', rightBrowPhogs);
[result, dists] = flann_search(rightBrowIndex, rightBrowPhog, neighborsCount, rightBrowParameters);
flann_free_index(rightBrowIndex);

for i=1:neighborsCount
    d = rightBrowImages{result(i,1)};
    imwrite(d, ['./patch/rightBrow' num2str(uint8(i)) '.jpg']);
    d = rgb2ycbcr(uint8(d));
    d = d(:,:,1);
    rightBrowDic(:,i) = double(d(:));
end

alpha = mexLasso(double(XrightBrow(:)),double(rightBrowDic),param);
ZrightBrow = reshape(rightBrowDic*alpha,size(XrightBrow));

%figure,imshow(XrightBrow);
%figure,imshow(uint8(ZrightBrow));

XrightBrow(:,:,2:3) = X(bigRightBrowRect(1,3):bigRightBrowRect(1,4),bigRightBrowRect(1,1):bigRightBrowRect(1,2),2:3);
%figure,imshow(ycbcr2rgb(XrightBrow));
ZrightBrow(:,:,2:3) = XrightBrow(:,:,2:3);
%figure,imshow(ycbcr2rgb(uint8(ZrightBrow)));

%{
im = zeros(size(XrightBrow,1), size(XrightBrow,2)*2+5, 3);
im(:,1:size(XrightBrow,2),:) = ycbcr2rgb(XrightBrow);
im(:,size(XrightBrow,2)+6:size(XrightBrow,2)*2+5,:) = ycbcr2rgb(uint8(ZrightBrow));
figure,imshow(uint8(im));
%}

%ZrightBrow(:,:,1) = global_color_correction(XrightBrow(:,:,1), ZrightBrow(:,:,1));

imwrite(ycbcr2rgb(XrightBrow), './patch/XrightBrow.jpg');
imwrite(ycbcr2rgb(uint8(ZrightBrow)), './patch/ZrightBrow.jpg');

im = zeros(size(XrightBrow,1), size(XrightBrow,2)*2+5, 3);
im(:,1:size(XrightBrow,2),:) = ycbcr2rgb(XrightBrow);
im(:,size(XrightBrow,2)+6:size(XrightBrow,2)*2+5,:) = ycbcr2rgb(uint8(ZrightBrow));
figure,imshow(uint8(im));

%{
% Shrinkage left brow
% Left -> right
smallLeftBrowRect(1,1) = floor(min(leftBrow(1,:))) - 1;
smallLeftBrowRect(1,2) = ceil(max(leftBrow(1,:))) + 1;
% Top -> bottom
smallLeftBrowRect(1,3) = floor(min(leftBrow(2,:))) - 1;
smallLeftBrowRect(1,4) = ceil(max(leftBrow(2,:))) + 1;

Ztmp = X;
Ztmp(bigLeftBrowRect(1,3):bigLeftBrowRect(1,4),bigLeftBrowRect(1,1):bigLeftBrowRect(1,2),:) = ZleftBrow(:,:,:);
ZleftBrow = Ztmp(smallLeftBrowRect(1,3):smallLeftBrowRect(1,4),smallLeftBrowRect(1,1):smallLeftBrowRect(1,2),:)

% Shrinkage right brow
% Left -> right
smallRightBrowRect(1,1) = floor(min(rightBrow(1,:))) - 1;
smallRightBrowRect(1,2) = ceil(max(rightBrow(1,:))) + 1;
% Top -> bottom
smallRightBrowRect(1,3) = floor(min(rightBrow(2,:))) - 1;
smallRightBrowRect(1,4) = ceil(max(rightBrow(2,:))) + 1;

Ztmp = X;
Ztmp(bigRightBrowRect(1,3):bigRightBrowRect(1,4),bigRightBrowRect(1,1):bigRightBrowRect(1,2),:) = ZrightBrow(:,:,:);
ZrightBrow = Ztmp(smallRightBrowRect(1,3):smallRightBrowRect(1,4),smallRightBrowRect(1,1):smallRightBrowRect(1,2),:)
%}
