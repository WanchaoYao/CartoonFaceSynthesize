function [Znose noseRect] = nose_synthesize(X, nose)

addpath(genpath(pwd));

bin = 8;
Angle = 360;
L = 3;

neighborsCount = 30;
param.lambda = 0.001;

% Synthesize nose
load './nose.mat';

% Middle point
middlePoint = [round((min(nose(1,:))+max(nose(1,:)))/2) round((max(nose(2,:))+((nose(2,7)+nose(2,8))/2))/2)];

% Left -> right
noseRect(1,1) = middlePoint(1,1) - round(globalNoseRectWidth/2);
noseRect(1,2) = middlePoint(1,1) + round(globalNoseRectWidth/2);
% Top -> bottom
noseRect(1,3) = middlePoint(1,2) - round(globalNoseRectHeight/2);
noseRect(1,4) = middlePoint(1,2) + round(globalNoseRectHeight/2);

Xnose = X(noseRect(1,3):noseRect(1,4),noseRect(1,1):noseRect(1,2),1);
noseROI(1,1) = noseRect(1,3);
noseROI(2,1) = noseRect(1,4);
noseROI(3,1) = noseRect(1,1);
noseROI(4,1) = noseRect(1,2);
nosePhog = anna_phog(X(:,:,1),bin,Angle,L,noseROI);

%[~,threshold] = edge(Xnose,'canny',[],1.0);
%k = 0.50;
%E = edge(Xnose,'canny',k*threshold,1.0);
%figure,imshow(E)

noseIndex = flann_load_index('./nose_index.index', nosePhogs);
[result, dists] = flann_search(noseIndex, nosePhog, neighborsCount, noseParameters);
flann_free_index(noseIndex);

for i=1:neighborsCount
    d = noseImages{result(i,1)};
    imwrite(d, ['./patch/nose' num2str(uint8(i)) '.jpg']);
    d = rgb2ycbcr(uint8(d));
    d = d(:,:,1);
    noseDic(:,i) = double(d(:));
    %d = edge(noseImages{i},'canny');
    %noseEdgeDic(:,i) = double(d(:));
end

alpha = mexLasso(double(Xnose(:)),double(noseDic),param);
Znose = reshape(noseDic*alpha,size(Xnose));

%figure,imshow(Xnose);
%figure,imshow(uint8(Znose));

Xnose(:,:,2:3) = X(noseRect(1,3):noseRect(1,4),noseRect(1,1):noseRect(1,2),2:3);
%figure,imshow(ycbcr2rgb(Xnose));
Znose(:,:,2:3) = Xnose(:,:,2:3);
%figure,imshow(ycbcr2rgb(uint8(Znose)));

imwrite(ycbcr2rgb(Xnose), './patch/Xnose.jpg');
imwrite(ycbcr2rgb(uint8(Znose)), './patch/Znose.jpg');

im = zeros(size(Xnose,1), size(Xnose,2)*2+5, 3);
im(:,1:size(Xnose,2),:) = ycbcr2rgb(Xnose);
im(:,size(Xnose,2)+6:size(Xnose,2)*2+5,:) = ycbcr2rgb(uint8(Znose));
figure,imshow(uint8(im));