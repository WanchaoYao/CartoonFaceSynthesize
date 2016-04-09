function [Zmouth mouthRect] = mouth_synthesize(X, mouth)

addpath(genpath(pwd));

bin = 8;
Angle = 360;
L = 3;

% Synthesize mouth
load './mouth.mat';

neighborsCount = 30;
param.lambda = 0.001;

% Middle point
middlePoint = [round((mouth(1,4)+mouth(1,10))/2) round((mouth(2,4)+mouth(2,10))/2)];
% Left -> right
mouthRect(1,1) = middlePoint(1,1) - round(globalMouthRectWidth/2);
mouthRect(1,2) = middlePoint(1,1) + round(globalMouthRectWidth/2);
% Top -> bottom
mouthRect(1,3) = middlePoint(1,2) - round(globalMouthRectHeight/2);
mouthRect(1,4) = middlePoint(1,2) + round(globalMouthRectHeight/2);

Xmouth = X(mouthRect(1,3):mouthRect(1,4),mouthRect(1,1):mouthRect(1,2),1);
mouthROI(1,1) = mouthRect(1,3);
mouthROI(2,1) = mouthRect(1,4);
mouthROI(3,1) = mouthRect(1,1);
mouthROI(4,1) = mouthRect(1,2);
mouthPhog = anna_phog(X(:,:,1),bin,Angle,L,mouthROI);

%[~,threshold] = edge(Xmouth,'canny',[],1.0);
%k = 0.50;
%E = edge(Xmouth,'canny',k*threshold,1.0);
%figure,imshow(E)

mouthIndex = flann_load_index('./mouth_index.index', mouthPhogs);
[result, dists] = flann_search(mouthIndex, mouthPhog, neighborsCount, mouthParameters);
flann_free_index(mouthIndex);

for i=1:neighborsCount
    d = mouthImages{result(i,1)};
    imwrite(d, ['./patch/mouth' num2str(uint8(i)) '.jpg']);
    d = rgb2ycbcr(uint8(d));
    d = d(:,:,1);
    mouthDic(:,i) = double(d(:));
    %d = edge(mouthImages{i},'canny');
    %mouthEdgeDic(:,i) = double(d(:));
end

alpha = mexLasso(double(Xmouth(:)),double(mouthDic),param);
Zmouth = reshape(mouthDic*alpha,size(Xmouth));

%figure,imshow(Xmouth);
%figure,imshow(uint8(Zmouth));

Xmouth(:,:,2:3) = X(mouthRect(1,3):mouthRect(1,4),mouthRect(1,1):mouthRect(1,2),2:3);
%figure,imshow(ycbcr2rgb(Xmouth));
Zmouth(:,:,2:3) = Xmouth(:,:,2:3);
%figure,imshow(ycbcr2rgb(uint8(Zmouth)));

imwrite(ycbcr2rgb(Xmouth), './patch/Xmouth.jpg');
imwrite(ycbcr2rgb(uint8(Zmouth)), './patch/Zmouth.jpg');

im = zeros(size(Xmouth,1), size(Xmouth,2)*2+5, 3);
%size(im)
%size(Xmouth)
im(:,1:size(Xmouth,2),:) = ycbcr2rgb(Xmouth);
im(:,size(Xmouth,2)+6:size(Xmouth,2)*2+5,:) = ycbcr2rgb(uint8(Zmouth));
figure,imshow(uint8(im));
