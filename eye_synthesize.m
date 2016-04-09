function [ZleftEye leftEyeRect ZrightEye rightEyeRect] = eye_synthesize(X, leftEye, rightEye)

addpath(genpath(pwd));

bin = 8;
Angle = 360;
L = 3;

neighborsCount = 30;
param.lambda = 0.001;

% Synthesize left eye
load './leftEye.mat';

% Middle point
leftMiddlePoint = [round((min(leftEye(1,:))+max(leftEye(1,:)))/2) round((max(leftEye(2,:))+min(leftEye(2,:)))/2)];

% Left -> right
leftEyeRect(1,1) = leftMiddlePoint(1,1) - round(globalLeftEyeRectWidth/2);
leftEyeRect(1,2) = leftMiddlePoint(1,1) + round(globalLeftEyeRectWidth/2);
% Top -> bottom
leftEyeRect(1,3) = leftMiddlePoint(1,2) - round(globalLeftEyeRectHeight/2);
leftEyeRect(1,4) = leftMiddlePoint(1,2) + round(globalLeftEyeRectHeight/2);

XleftEye = X(leftEyeRect(1,3):leftEyeRect(1,4),leftEyeRect(1,1):leftEyeRect(1,2),1);
leftEyeROI(1,1) = leftEyeRect(1,3);
leftEyeROI(2,1) = leftEyeRect(1,4);
leftEyeROI(3,1) = leftEyeRect(1,1);
leftEyeROI(4,1) = leftEyeRect(1,2);
leftEyePhog = anna_phog(X(:,:,1),bin,Angle,L,leftEyeROI);

leftEyeIndex = flann_load_index('./leftEye_index.index', leftEyePhogs);
[result, dists] = flann_search(leftEyeIndex, leftEyePhog, neighborsCount, leftEyeParameters);
flann_free_index(leftEyeIndex);

for i=1:neighborsCount
    d = leftEyeImages{result(i,1)};
    imwrite(d, ['./patch/leftEye' num2str(uint8(i)) '.jpg']);
    d = rgb2ycbcr(uint8(d));
    d = d(:,:,1);
    leftEyeDic(:,i) = double(d(:));
end

alpha = mexLasso(double(XleftEye(:)),double(leftEyeDic),param);
ZleftEye = reshape(leftEyeDic*alpha,size(XleftEye));

%figure,imshow(XleftEye);
%figure,imshow(uint8(ZleftEye));

XleftEye(:,:,2:3) = X(leftEyeRect(1,3):leftEyeRect(1,4),leftEyeRect(1,1):leftEyeRect(1,2),2:3);
%figure,imshow(ycbcr2rgb(XleftEye));
ZleftEye(:,:,2:3) = XleftEye(:,:,2:3);
%figure,imshow(ycbcr2rgb(uint8(ZleftEye)));

imwrite(ycbcr2rgb(XleftEye), './patch/XleftEye.jpg');
imwrite(ycbcr2rgb(uint8(ZleftEye)), './patch/ZleftEye.jpg');

im = zeros(size(XleftEye,1), size(XleftEye,2)*2+5, 3);
im(:,1:size(XleftEye,2),:) = ycbcr2rgb(XleftEye);
im(:,size(XleftEye,2)+6:size(XleftEye,2)*2+5,:) = ycbcr2rgb(uint8(ZleftEye));
figure,imshow(uint8(im));

% Synthesize right eye
load './rightEye.mat';

% Middle point
rightMiddlePoint = [round((min(rightEye(1,:))+max(rightEye(1,:)))/2) round((max(rightEye(2,:))+min(rightEye(2,:)))/2)];

% Left -> right
rightEyeRect(1,1) = rightMiddlePoint(1,1) - round(globalRightEyeRectWidth/2);
rightEyeRect(1,2) = rightMiddlePoint(1,1) + round(globalRightEyeRectWidth/2);
% Top -> bottom
rightEyeRect(1,3) = rightMiddlePoint(1,2) - round(globalRightEyeRectHeight/2);
rightEyeRect(1,4) = rightMiddlePoint(1,2) + round(globalRightEyeRectHeight/2);

XrightEye = X(rightEyeRect(1,3):rightEyeRect(1,4),rightEyeRect(1,1):rightEyeRect(1,2),1);
rightEyeROI(1,1) = rightEyeRect(1,3);
rightEyeROI(2,1) = rightEyeRect(1,4);
rightEyeROI(3,1) = rightEyeRect(1,1);
rightEyeROI(4,1) = rightEyeRect(1,2);
rightEyePhog = anna_phog(X(:,:,1),bin,Angle,L,rightEyeROI);

rightEyeIndex = flann_load_index('./rightEye_index.index', rightEyePhogs);
neighborsCount = 30;
[result, dists] = flann_search(rightEyeIndex, rightEyePhog, neighborsCount, rightEyeParameters);
flann_free_index(rightEyeIndex);

for i=1:neighborsCount
    d = rightEyeImages{result(i,1)};
    imwrite(d, ['./patch/rightEye' num2str(uint8(i)) '.jpg']);
    d = rgb2ycbcr(uint8(d));
    d = d(:,:,1);
    rightEyeDic(:,i) = double(d(:));
end

alpha = mexLasso(double(XrightEye(:)),double(rightEyeDic),param);
ZrightEye = reshape(rightEyeDic*alpha,size(XrightEye));

%figure,imshow(XrightEye);
%figure,imshow(uint8(ZrightEye));

XrightEye(:,:,2:3) = X(rightEyeRect(1,3):rightEyeRect(1,4),rightEyeRect(1,1):rightEyeRect(1,2),2:3);
%figure,imshow(ycbcr2rgb(XrightEye));
ZrightEye(:,:,2:3) = XrightEye(:,:,2:3);
%figure,imshow(ycbcr2rgb(uint8(ZrightEye)));

imwrite(ycbcr2rgb(XrightEye), './patch/XrightEye.jpg');
imwrite(ycbcr2rgb(uint8(ZrightEye)), './patch/ZrightEye.jpg');

im = zeros(size(XrightEye,1), size(XrightEye,2)*2+5, 3);
im(:,1:size(XrightEye,2),:) = ycbcr2rgb(XrightEye);
im(:,size(XrightEye,2)+6:size(XrightEye,2)*2+5,:) = ycbcr2rgb(uint8(ZrightEye));
figure,imshow(uint8(im));