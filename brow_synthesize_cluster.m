function [ZleftBrow leftBrowRect ZrightBrow rightBrowRect] = ...
    brow_synthesize_cluster(X, leftBrow, rightBrow,...
    maxY,minY,maxCb,minCb,maxCr,minCr,meanFaceValue)

addpath(genpath(pwd));

% Synthesize left brow
% Left -> right
leftBrowRect(1,1) = floor(min(leftBrow(1,:))) - 2;
leftBrowRect(1,2) = ceil(max(leftBrow(1,:))) + 2;
% Top -> bottom
leftBrowRect(1,3) = floor(min(leftBrow(2,:))) - 2;
leftBrowRect(1,4) = ceil(max(leftBrow(2,:))) + 2;

ZleftBrow = zeros(leftBrowRect(1,4)-leftBrowRect(1,3)+1,...
    leftBrowRect(1,2)-leftBrowRect(1,1)+1,3);
meanLeftBrowValue = zeros(1,1,3);
leftBrowCount = 0;
for i = 1:leftBrowRect(1,4)-leftBrowRect(1,3)+1
    for j = 1:leftBrowRect(1,2)-leftBrowRect(1,1)+1
        ii = i+leftBrowRect(1,3)-1;
        jj = j+leftBrowRect(1,1)-1;
        if (X(ii,jj,1) >= minY && X(ii,jj,1) <= maxY && ...
            X(ii,jj,2) >= minCb && X(ii,jj,2) <= maxCb && ... 
            X(ii,jj,3) >= minCr && X(ii,jj,3) <= maxCr)
            ZleftBrow(i,j,:) = meanFaceValue(1,1,:);
        else
            leftBrowCount = leftBrowCount + 1;
            meanLeftBrowValue(1,1,:) = meanLeftBrowValue(1,1,:) + ...
                double(X(ii,jj,:));
        end
    end
end
meanFaceValue
meanLeftBrowValue(1,1,:) = round(meanLeftBrowValue(1,1,:)/leftBrowCount);

for i = 1:leftBrowRect(1,4)-leftBrowRect(1,3)+1
    for j = 1:leftBrowRect(1,2)-leftBrowRect(1,1)+1
        ii = i+leftBrowRect(1,3)-1;
        jj = j+leftBrowRect(1,1)-1;
        if (X(ii,jj,1) >= minY && X(ii,jj,1) <= maxY && ...
            X(ii,jj,2) >= minCb && X(ii,jj,2) <= maxCb && ... 
            X(ii,jj,3) >= minCr && X(ii,jj,3) <= maxCr)
            ZleftBrow(i,j,:) = meanFaceValue(1,1,:);
        else
            ZleftBrow(i,j,:) = meanLeftBrowValue(1,1,:);
        end
    end
end

XleftBrow = X(leftBrowRect(1,3):leftBrowRect(1,4),leftBrowRect(1,1):leftBrowRect(1,2),1:3);

im = zeros(size(XleftBrow,1), size(XleftBrow,2)*2+5, 3);
im(:,1:size(XleftBrow,2),:) = ycbcr2rgb(XleftBrow);
im(:,size(XleftBrow,2)+6:size(XleftBrow,2)*2+5,:) = ycbcr2rgb(uint8(ZleftBrow));
figure,imshow(uint8(im));

% Synthesize right brow
% Left -> right
rightBrowRect(1,1) = floor(min(rightBrow(1,:))) - 2;
rightBrowRect(1,2) = ceil(max(rightBrow(1,:))) + 2;
% Top -> bottom
rightBrowRect(1,3) = floor(min(rightBrow(2,:))) - 2;
rightBrowRect(1,4) = ceil(max(rightBrow(2,:))) + 2;

ZrightBrow = zeros(rightBrowRect(1,4)-rightBrowRect(1,3)+1,...
    rightBrowRect(1,2)-rightBrowRect(1,1)+1,3);
meanRightBrowValue = zeros(1,1,3);
rightBrowCount = 0;
for i = 1:rightBrowRect(1,4)-rightBrowRect(1,3)+1
    for j = 1:rightBrowRect(1,2)-rightBrowRect(1,1)+1
        ii = i+rightBrowRect(1,3)-1;
        jj = j+rightBrowRect(1,1)-1;
        if (X(ii,jj,1) >= minY && X(ii,jj,1) <= maxY && ...
            X(ii,jj,2) >= minCb && X(ii,jj,2) <= maxCb && ... 
            X(ii,jj,3) >= minCr && X(ii,jj,3) <= maxCr)
            ZrightBrow(i,j,:) = meanFaceValue(1,1,:);
        else
            rightBrowCount = rightBrowCount + 1;
            meanRightBrowValue(1,1,:) = meanRightBrowValue(1,1,:) + ...
                double(X(ii,jj,:));
        end
    end
end

meanRightBrowValue(1,1,:) = round(meanRightBrowValue(1,1,:)/rightBrowCount);

for i = 1:rightBrowRect(1,4)-rightBrowRect(1,3)+1
    for j = 1:rightBrowRect(1,2)-rightBrowRect(1,1)+1
        ii = i+rightBrowRect(1,3)-1;
        jj = j+rightBrowRect(1,1)-1;
        if (X(ii,jj,1) >= minY && X(ii,jj,1) <= maxY && ...
            X(ii,jj,2) >= minCb && X(ii,jj,2) <= maxCb && ... 
            X(ii,jj,3) >= minCr && X(ii,jj,3) <= maxCr)
            ZrightBrow(i,j,:) = meanFaceValue(1,1,:);
        else
            ZrightBrow(i,j,:) = meanRightBrowValue(1,1,:);
        end
    end
end

XrightBrow = X(rightBrowRect(1,3):rightBrowRect(1,4),rightBrowRect(1,1):rightBrowRect(1,2),1:3);

im = zeros(size(XrightBrow,1), size(XrightBrow,2)*2+5, 3);
im(:,1:size(XrightBrow,2),:) = ycbcr2rgb(XrightBrow);
im(:,size(XrightBrow,2)+6:size(XrightBrow,2)*2+5,:) = ycbcr2rgb(uint8(ZrightBrow));
figure,imshow(uint8(im));