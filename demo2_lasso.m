function demo2_lasso(photoPath)

tic

addpath(genpath(pwd));

neighborsCount = 20;

% Load data and index
load './demo2.mat'
index = flann_load_index('./demo2_index.index', trainingSet);

[X faceProfile leftEye rightEye leftBrow rightBrow nose mouth angle] = pre_process(photoPath, true, false);

if isempty(X)
    disp('Invalid input photo!');
    return;
end

toc
disp('Face dectection done!');

X = rgb2ycbcr(X);

[clustered,C] = pixel_cluster(X,3);

noOverlapPatchSize = 16;
noOverlapPatchesRowCount = ylength/noOverlapPatchSize;
noOverlapPatchesColCount = xlength/noOverlapPatchSize;
noOverlapPatches = cell(noOverlapPatchesRowCount,noOverlapPatchesColCount);

for i = 1:noOverlapPatchesRowCount
	for j = 1:noOverlapPatchesColCount
		noOverlapPatches{i,j} = [noOverlapPatchSize*(i-1)+1;noOverlapPatchSize*(i-1)+noOverlapPatchSize;...
            noOverlapPatchSize*(j-1)+1;noOverlapPatchSize*(j-1)+noOverlapPatchSize];
	end
end

% Get name and extension of the input photo
[~, name, ext] = fileparts(photoPath);

% Compute PHOG over all patches the input photo
D = my_phog(X(:,:,1),bin,Angle,L,noOverlapPatches);

% Construct test set
testSet = zeros(size(phogs{1}{1,1},1),noOverlapPatchesRowCount*noOverlapPatchesColCount);
for i = 1:noOverlapPatchesRowCount
    for j = 1:noOverlapPatchesColCount
        testSet(:,(i-1)*noOverlapPatchesColCount+j) = D{i,j}; 
    end
end

disp('Start searching nearest neighbors.');
toc

[result, dists] = flann_search(index, testSet, neighborsCount, parameters);
flann_free_index(index);

disp('Start synthesizing cartoon.');
toc

Dic = zeros(noOverlapPatchSize*noOverlapPatchSize, neighborsCount);
Zp = cell(noOverlapPatchesRowCount, noOverlapPatchesColCount);
Xp = cell(noOverlapPatchesRowCount, noOverlapPatchesColCount);
for m = 1:noOverlapPatchesRowCount
	for n = 1:noOverlapPatchesColCount
        % Construct dictionary
        patchIndexes = result(:,(m-1)*noOverlapPatchesColCount+n);
        for count = 1:neighborsCount
            patchIndex = patchIndexes(count,1);
            % For test
            % patchIndex
            patchIndex = patchIndex + patchRowCount * patchColCount + patchColCount;
            k = floor(patchIndex/(patchRowCount*patchColCount));
            patchIndex = patchIndex - k * patchRowCount * patchColCount;
			i = floor(patchIndex/patchColCount);
			j = patchIndex - i * patchColCount;
            if 0 == j
                i = i - 1;
                j = patchColCount;
            end
            if 0 >= i
                k = k -1;
                i = patchRowCount + i;
            end
            % For test    
            % k,i,j
            if (k-1)*patchRowCount*patchColCount+(i-1)*patchColCount+j ~= patchIndexes(count,1)
                disp('Error!');
                pause
            end
            patch = cartoons{k}(patches{i,j}(1,1):patches{i,j}(2,1),patches{i,j}(3,1):patches{i,j}(4,1),1);
            Dic(:,count) = patch(:);
            
            %{
            if m == noOverlapPatchesRowCount - 1 && n == 2
                patchIndex, k, i, j, patch
                figure, imshow(uint8(patch));
                drawnow;
                pause;
            end
            %}
        end
        
        Xp{m,n} = X(noOverlapPatches{m,n}(1,1):noOverlapPatches{m,n}(2,1),noOverlapPatches{m,n}(3,1):noOverlapPatches{m,n}(4,1),1);
        
        % Lasso
        param.lambda = 0.001;
        alpha = mexLasso(double(Xp{m,n}(:)),double(Dic),param);

        Zp{m,n} = reshape(Dic*alpha,noOverlapPatchSize,noOverlapPatchSize);
        
        %if m == 1 && n == 2
        %Xp{m,n}, Zp{m,n}, alpha
        %end
        %pause
	end
end

% Show/Save and compare photo patches and synthesized cartoon patches
for m = 1:noOverlapPatchesRowCount
	for n = 1:noOverlapPatchesColCount
		A(noOverlapPatchSize*(m-1)+1:noOverlapPatchSize*m,noOverlapPatchSize*(n-1)+1:noOverlapPatchSize*n) = Xp{m,n};
		B(noOverlapPatchSize*(m-1)+1:noOverlapPatchSize*m,noOverlapPatchSize*(n-1)+1:noOverlapPatchSize*n) = Zp{m,n};
	end
end
% Show
h = figure;
set(h,'name','Comparison between photo patches and synthesized cartoon patches');
subplot(1,2,1),imshow(uint8(A)),title('Patches of photo');
subplot(1,2,2),imshow(uint8(B)),title('Patches of synthesized cartoon');
%truesize;
% Save
im = frame2im(getframe(h));
imwrite(im,['./synthesis_lasso/',name,'_patchesComparison',ext],ext(2:length(ext)));

%pause

disp('Start putting together all patches.');
toc

% Put together all patches to synthesize cartoon Z
Z = zeros(ylength,xlength,3);
for m = 1:noOverlapPatchesRowCount
	for n = 1:noOverlapPatchesColCount
		Z(noOverlapPatches{m,n}(1,1):noOverlapPatches{m,n}(2,1),noOverlapPatches{m,n}(3,1):noOverlapPatches{m,n}(4,1),1) = Zp{m,n};
	end
end

% Color-Preserving Colorization
Z(:,:,2:3) = X(:,:,2:3);
figure, imshow(ycbcr2rgb(uint8(Z)));

Z_ori = Z;

% Synthesize mouth
[Zmouth mouthRect] = mouth_synthesize(X, mouth);

% Synthesize nose
[Znose noseRect] = nose_synthesize(X, nose);

% Synthesize eyes
[ZleftEye leftEyeRect ZrightEye rightEyeRect] = eye_synthesize(X, leftEye, rightEye);

% Synthesize brows
[ZleftBrow leftBrowRect ZrightBrow rightBrowRect] = brow_synthesize(X, leftBrow, rightBrow);

% Draw face
xv = [faceProfile(1,:) faceProfile(1,1)];
yv = [faceProfile(2,:) faceProfile(2,1)];
faceMiddleH = floor(max(min(yv(1:8)),min(yv(9:16))));
faceMaxH = ceil(max(yv));
eyeMinH = min([rightEyeRect(1,3) leftEyeRect(1,3)]);
faceMiddleH = max([faceMiddleH eyeMinH]);
faceMiddleH
faceMaxH
faceProfile(1,1)
faceProfile(1,17)
partClustered = clustered(faceMiddleH:faceMaxH,faceProfile(1,1):faceProfile(1,17),4);
faceCluster = mode(double(partClustered(:)));
correctFaceCount = 0;
meanValue = zeros(1,1,3);
maxY = 0;
minY = 255;
maxCb = 0;
minCb = 255;
maxCr = 0;
minCr = 255;
for i = faceMiddleH:faceMaxH
    for j = floor(faceProfile(1,1)):ceil(faceProfile(1,17))
        %{
        if i > mouthRect(1,3) && i < mouthRect(1,4) && ...
                j > mouthRect(1,1) && j < mouthRect(1,2)
            continue;
        end
        if i > noseRect(1,3) && i < noseRect(1,4) && ...
                j > noseRect(1,1) && j < noseRect(1,2)
            continue;
        end
        if i > leftEyeRect(1,3) && i < leftEyeRect(1,4) && ...
                j > leftEyeRect(1,1) && j < leftEyeRect(1,2)
            continue;
        end
        if i > rightEyeRect(1,3) && i < rightEyeRect(1,4) && ...
                j > rightEyeRect(1,1) && j < rightEyeRect(1,2)
            continue;
        end
        if i > leftBrowRect(1,3) && i < leftBrowRect(1,4) && ...
                j > leftBrowRect(1,1) && j < leftBrowRect(1,2)
            continue;
        end
        if i > rightBrowRect(1,3) && i < rightBrowRect(1,4) && ...
                j > rightBrowRect(1,1) && j < rightBrowRect(1,2)
            continue;
        end
        %}
        [in on] = inpolygon(j,i,xv,yv);
        if (in == 1 || on == 1 ) && faceCluster == clustered(i,j,4)
            meanValue(1,1,:) = meanValue(1,1,:) + double(X(i,j,:));
            correctFaceCount = correctFaceCount+1;
            if maxCb < X(i,j,2)
                maxCb = X(i,j,2);
            end
            if maxCr < X(i,j,3)
                maxCr = X(i,j,3);
            end
            if maxY < X(i,j,1)
                maxY = X(i,j,1);
            end
            if minCb > X(i,j,2)
                minCb = X(i,j,2);
            end
            if minCr > X(i,j,3)
                minCr = X(i,j,3);
            end
            if minY > X(i,j,1)
                minY = X(i,j,1);
            end
        end
    end
end
meanValue(1,1,:) = meanValue(1,1,:) / correctFaceCount;
meanValue = round(meanValue);

faceMinH = round(faceMiddleH-1.2*(faceMaxH-faceMiddleH));
if faceMinH <= 0
    faceMinH = 1;
end

%{
maxY,minY
maxCb,minCb
maxCr,minCr
faceMinH
faceMiddleH
faceMaxH
pause
%}

faceBoard = zeros(ylength,xlength);
faceBoard(:,:,1) = meanValue(1,1,1);
faceBoard(:,:,2) = meanValue(1,1,2);
faceBoard(:,:,3) = meanValue(1,1,3);

%imMask = zeros(ylength,xlength);

for i = faceMinH:faceMaxH
    for j = floor(faceProfile(1,1)):ceil(faceProfile(1,17))
        if i > faceMiddleH
            [in on] = inpolygon(j,i,xv,yv);
            %if (in == 1 || on == 1 ) && clustered(i,j,4) == faceCluster
            if (in == 1 || on == 1 ) && ...
                    (X(i,j,1) >= minY && X(i,j,1) <= maxY && ...
                    X(i,j,2) >= minCb && X(i,j,2) <= maxCb && ...
                    X(i,j,3) >= minCr && X(i,j,3) <= maxCr || ...
                    clustered(i,j,4) == faceCluster)
                Z(i,j,:) = meanValue(1,1,:);
                %imMask(i,j) = 1;
            end
        else
            %if clustered(i,j,4) == faceCluster
            if X(i,j,1) >= minY && X(i,j,1) <= maxY && ...
                    X(i,j,2) >= minCb && X(i,j,2) <= maxCb && ... 
                    X(i,j,3) >= minCr && X(i,j,3) <= maxCr || ...
                    clustered(i,j,4) == faceCluster
                Z(i,j,:) = meanValue(1,1,:);
                %imMask(i,j) = 1;
            end
        end
    end
end

% Draw chin profile
faceProfile = round(faceProfile);
for j = 6:11
    bottom = faceProfile(1,j):1:faceProfile(1,j+1);
    funcBottom = polyfit(faceProfile(1,j:j+1), faceProfile(2,j:j+1), 1);  %1´ÎÄâºÏ
    resultBottom = round(polyval(funcBottom, bottom, 1));
    for i = 1:length(bottom)
        Z(resultBottom(i),bottom(i),1) = 63;
    end
end

imwrite(ycbcr2rgb(uint8(Z)), './patch/face.jpg');

% Synthesize brows
%[ZleftBrow leftBrowRect ZrightBrow rightBrowRect] = brow_synthesize_cluster(X, leftBrow, rightBrow,...
%    maxY,minY,maxCb,minCb,maxCr,minCr,meanValue);

XleftBrow = X(leftBrowRect(1,3):leftBrowRect(1,4),leftBrowRect(1,1):leftBrowRect(1,2),1);
%{
[clusteredLeftBrow,clusteredLeftBrowC] = pixel_cluster(XleftBrow(:,:,1),3);
minC = round(min(clusteredLeftBrowC(:)));
for i = 1:size(ZleftBrow,1)
    for j = 1:size(ZleftBrow,2)
        if clusteredLeftBrow(i,j,1) == minC
            Z(i+leftBrowRect(1,3)-1,j+leftBrowRect(1,1)-1,1) = ZleftBrow(i,j,1);
        end
    end
end
%}

XrightBrow = X(rightBrowRect(1,3):rightBrowRect(1,4),rightBrowRect(1,1):rightBrowRect(1,2),1);
%{
[clusteredRightBrow,clusteredRightBrowC] = pixel_cluster(XrightBrow(:,:,1),3);
minC = round(min(clusteredRightBrowC(:)));
for i = 1:size(ZrightBrow,1)
    for j = 1:size(ZrightBrow,2)
        if clusteredRightBrow(i,j,1) == minC
            Z(i+rightBrowRect(1,3)-1,j+rightBrowRect(1,1)-1,1) = ZrightBrow(i,j,1);
        end
    end
end
%}

Z = double(ycbcr2rgb(uint8(Z)));
faceBoard = double(ycbcr2rgb(uint8(faceBoard)));
Zmouth = double(ycbcr2rgb(uint8(Zmouth)));
Znose = double(ycbcr2rgb(uint8(Znose)));
ZleftEye = double(ycbcr2rgb(uint8(ZleftEye)));
ZrightEye = double(ycbcr2rgb(uint8(ZrightEye)));
ZleftBrow = double(ycbcr2rgb(uint8(ZleftBrow)));
ZrightBrow = double(ycbcr2rgb(uint8(ZrightBrow)));

figure, imshow(uint8(Z));
figure, imshow(uint8(faceBoard));
%pause;

%Z(:,:,1) = poissonSolverSeamlessCloning1(faceBoard(:,:,1), Z(:,:,1), imMask, [0 0]);
%Z(:,:,2) = poissonSolverSeamlessCloning1(faceBoard(:,:,2), Z(:,:,2), imMask, [0 0]);
%Z(:,:,3) = poissonSolverSeamlessCloning1(faceBoard(:,:,3), Z(:,:,3), imMask, [0 0]);

%imwrite(uint8(faceBoard), '/Users/user/Downloads/face.jpg');
%imwrite(uint8(Znose), '/Users/user/Downloads/nose.jpg');

% Change mouth
%size(Zmouth)
[imRegion imNew(:,:,1)] = ...
	my_poissonSolverSeamlessCloning(Zmouth(:,:,1),Z(:,:,1),...
    [2 size(Zmouth,2)-2 2 size(Zmouth,1)-2],[mouthRect(1,1)+1 mouthRect(1,3)+1]);
Z(mouthRect(1,3)+2:mouthRect(1,4)-1,mouthRect(1,1)+2:mouthRect(1,2)-1,1) = imRegion;
[imRegion imNew(:,:,2)] = ...
    my_poissonSolverSeamlessCloning(Zmouth(:,:,2),Z(:,:,2),...
    [2 size(Zmouth,2)-2 2 size(Zmouth,1)-2],[mouthRect(1,1)+1 mouthRect(1,3)+1]);
Z(mouthRect(1,3)+2:mouthRect(1,4)-1,mouthRect(1,1)+2:mouthRect(1,2)-1,2) = imRegion;
[imRegion imNew(:,:,3)] = ...
    my_poissonSolverSeamlessCloning(Zmouth(:,:,3),Z(:,:,3),...
    [2 size(Zmouth,2)-2 2 size(Zmouth,1)-2],[mouthRect(1,1)+1 mouthRect(1,3)+1]);
Z(mouthRect(1,3)+2:mouthRect(1,4)-1,mouthRect(1,1)+2:mouthRect(1,2)-1,3) = imRegion;
ZmouthSave = Z(mouthRect(1,3)+2:mouthRect(1,4)-1,mouthRect(1,1)+2:mouthRect(1,2)-1,:);

% Change nose
Z(noseRect(1,3)+2:noseRect(1,4)-1,noseRect(1,1)+2:noseRect(1,2)-1,1) = ...
    my_poissonSolverSeamlessCloning(Znose(:,:,1),Z(:,:,1),...
    [2 size(Znose,2)-2 2 size(Znose,1)-2],[noseRect(1,1)+1 noseRect(1,3)+1]);
Z(noseRect(1,3)+2:noseRect(1,4)-1,noseRect(1,1)+2:noseRect(1,2)-1,2) = ...
    my_poissonSolverSeamlessCloning(Znose(:,:,2),Z(:,:,2),...
    [2 size(Znose,2)-2 2 size(Znose,1)-2],[noseRect(1,1)+1 noseRect(1,3)+1]);
Z(noseRect(1,3)+2:noseRect(1,4)-1,noseRect(1,1)+2:noseRect(1,2)-1,3) = ...
    my_poissonSolverSeamlessCloning(Znose(:,:,3),Z(:,:,3),...
    [2 size(Znose,2)-2 2 size(Znose,1)-2],[noseRect(1,1)+1 noseRect(1,3)+1]);
ZnoseSave = Z(noseRect(1,3)+2:noseRect(1,4)-1,noseRect(1,1)+2:noseRect(1,2)-1,:);

% Change eyes
Z(leftEyeRect(1,3)+2:leftEyeRect(1,4)-1,leftEyeRect(1,1)+2:leftEyeRect(1,2)-1,1) = ...
    my_poissonSolverSeamlessCloning(ZleftEye(:,:,1),Z(:,:,1),...
    [2 size(ZleftEye,2)-2 2 size(ZleftEye,1)-2],[leftEyeRect(1,1)+1 leftEyeRect(1,3)+1]);
Z(leftEyeRect(1,3)+2:leftEyeRect(1,4)-1,leftEyeRect(1,1)+2:leftEyeRect(1,2)-1,2) = ...
    my_poissonSolverSeamlessCloning(ZleftEye(:,:,2),Z(:,:,2),...
    [2 size(ZleftEye,2)-2 2 size(ZleftEye,1)-2],[leftEyeRect(1,1)+1 leftEyeRect(1,3)+1]);
Z(leftEyeRect(1,3)+2:leftEyeRect(1,4)-1,leftEyeRect(1,1)+2:leftEyeRect(1,2)-1,3) = ...
    my_poissonSolverSeamlessCloning(ZleftEye(:,:,3),Z(:,:,3),...
    [2 size(ZleftEye,2)-2 2 size(ZleftEye,1)-2],[leftEyeRect(1,1)+1 leftEyeRect(1,3)+1]);
ZleftEyeSave = Z(leftEyeRect(1,3)+2:leftEyeRect(1,4)-1,leftEyeRect(1,1)+2:leftEyeRect(1,2)-1,:);

Z(rightEyeRect(1,3)+2:rightEyeRect(1,4)-1,rightEyeRect(1,1)+2:rightEyeRect(1,2)-1,1) = ...
    my_poissonSolverSeamlessCloning(ZrightEye(:,:,1),Z(:,:,1),...
    [2 size(ZrightEye,2)-2 2 size(ZrightEye,1)-2],[rightEyeRect(1,1)+1 rightEyeRect(1,3)+1]);
Z(rightEyeRect(1,3)+2:rightEyeRect(1,4)-1,rightEyeRect(1,1)+2:rightEyeRect(1,2)-1,2) = ...
    my_poissonSolverSeamlessCloning(ZrightEye(:,:,2),Z(:,:,2),...
    [2 size(ZrightEye,2)-2 2 size(ZrightEye,1)-2],[rightEyeRect(1,1)+1 rightEyeRect(1,3)+1]);
Z(rightEyeRect(1,3)+2:rightEyeRect(1,4)-1,rightEyeRect(1,1)+2:rightEyeRect(1,2)-1,3) = ...
    my_poissonSolverSeamlessCloning(ZrightEye(:,:,3),Z(:,:,3),...
    [2 size(ZrightEye,2)-2 2 size(ZrightEye,1)-2],[rightEyeRect(1,1)+1 rightEyeRect(1,3)+1]);
ZrightEyeSave = Z(rightEyeRect(1,3)+2:rightEyeRect(1,4)-1,rightEyeRect(1,1)+2:rightEyeRect(1,2)-1,:);

% Change brows
boxSrc = [2 size(ZleftBrow,2)-2 2 size(ZleftBrow,1)-2];
posDest = [leftBrowRect(1,1)+1 leftBrowRect(1,3)+1];
Z(leftBrowRect(1,3)+2:leftBrowRect(1,4)-1,leftBrowRect(1,1)+2:leftBrowRect(1,2)-1,1) = ...
    my_poissonSolverSeamlessCloning(ZleftBrow(:,:,1),Z(:,:,1),boxSrc,posDest);
Z(leftBrowRect(1,3)+2:leftBrowRect(1,4)-1,leftBrowRect(1,1)+2:leftBrowRect(1,2)-1,2) = ...
    my_poissonSolverSeamlessCloning(ZleftBrow(:,:,2),Z(:,:,2),boxSrc,posDest);
Z(leftBrowRect(1,3)+2:leftBrowRect(1,4)-1,leftBrowRect(1,1)+2:leftBrowRect(1,2)-1,3) = ...
    my_poissonSolverSeamlessCloning(ZleftBrow(:,:,3),Z(:,:,3),boxSrc,posDest);

boxSrc = [2 size(ZrightBrow,2)-2 2 size(ZrightBrow,1)-2];
posDest = [rightBrowRect(1,1)+1 rightBrowRect(1,3)+1];
%figure, imshow(uint8(Z));
Z(rightBrowRect(1,3)+2:rightBrowRect(1,4)-1,rightBrowRect(1,1)+2:rightBrowRect(1,2)-1,1) = ...
    my_poissonSolverSeamlessCloning(ZrightBrow(:,:,1),Z(:,:,1),boxSrc,posDest);
%figure, imshow(uint8(Z));
Z(rightBrowRect(1,3)+2:rightBrowRect(1,4)-1,rightBrowRect(1,1)+2:rightBrowRect(1,2)-1,2) = ...
    my_poissonSolverSeamlessCloning(ZrightBrow(:,:,2),Z(:,:,2),boxSrc,posDest);
%figure, imshow(uint8(Z));
Z(rightBrowRect(1,3)+2:rightBrowRect(1,4)-1,rightBrowRect(1,1)+2:rightBrowRect(1,2)-1,3) = ...
    my_poissonSolverSeamlessCloning(ZrightBrow(:,:,3),Z(:,:,3),boxSrc,posDest);
%figure, imshow(uint8(Z));

%{
meanValue = double(ycbcr2rgb(uint8(meanValue)));

leftEarRect(1,1) = floor(faceProfile(1,1))-10;
if leftEarRect(1,1) <= 0
    leftEarRect(1,1) = 1;
end
leftEarRect(1,2) = floor(faceProfile(1,1));
leftEarRect(1,3) = leftBrowRect(1,3);
leftEarRect(1,4) = mouthRect(1,4);
leftEarCount = 0;
leftEarMeanValue = zeros(1,1,3);
for i = leftEarRect(1,3):leftEarRect(1,4)
    for j = leftEarRect(1,1):leftEarRect(1,2)
        [in on] = inpolygon(j,i,xv,yv);
        if (in == 0 || on == 0 ) && ...
                (X(i,j,1) >= minY && X(i,j,1) <= maxY && ...
                X(i,j,2) >= minCb && X(i,j,2) <= maxCb && ...
                X(i,j,3) >= minCr && X(i,j,3) <= maxCr || ...
                clustered(i,j,4) == faceCluster)
            leftEarCount = leftEarCount + 1;
            leftEarMeanValue(1,1,:) = leftEarMeanValue(1,1,:) + Z(i,j,:);
        end
    end
end
leftEarMeanValue(1,1,:) = round(leftEarMeanValue(1,1,:) / leftEarCount)
for i = leftEarRect(1,3):leftEarRect(1,4)
    for j = leftEarRect(1,1):leftEarRect(1,2)
        [in on] = inpolygon(j,i,xv,yv);
        if (in == 0 || on == 0 ) && ...
                (X(i,j,1) >= minY && X(i,j,1) <= maxY && ...
                X(i,j,2) >= minCb && X(i,j,2) <= maxCb && ...
                X(i,j,3) >= minCr && X(i,j,3) <= maxCr || ...
                clustered(i,j,4) == faceCluster)
            Z(i,j,:) = Z(i,j,:) - (leftEarMeanValue(1,1,:) - meanValue(1,1,:));
        end
    end
end

rightEarRect(1,1) = ceil(faceProfile(1,17));
rightEarRect(1,2) = ceil(faceProfile(1,17))+10;
if rightEarRect(1,2) > xlength
    rightEarRect(1,2) = xlength;
end
rightEarRect(1,3) = rightBrowRect(1,3);
rightEarRect(1,4) = mouthRect(1,4);
rightEarCount = 0;
rightEarMeanValue = zeros(1,1,3);
for i = rightEarRect(1,3):rightEarRect(1,4)
    for j = rightEarRect(1,1):rightEarRect(1,2)
        [in on] = inpolygon(j,i,xv,yv);
        if (in == 0 || on == 0 ) && ...
                (X(i,j,1) >= minY && X(i,j,1) <= maxY && ...
                X(i,j,2) >= minCb && X(i,j,2) <= maxCb && ...
                X(i,j,3) >= minCr && X(i,j,3) <= maxCr || ...
                clustered(i,j,4) == faceCluster)
            rightEarCount = rightEarCount + 1;
            rightEarMeanValue(1,1,:) = rightEarMeanValue(1,1,:) + Z(i,j,:);
        end
    end
end
rightEarMeanValue(1,1,:) = round(rightEarMeanValue(1,1,:) / rightEarCount)
for i = rightEarRect(1,3):rightEarRect(1,4)
    for j = rightEarRect(1,1):rightEarRect(1,2)
        [in on] = inpolygon(j,i,xv,yv);
        if (in == 0 || on == 0 ) && ...
                (X(i,j,1) >= minY && X(i,j,1) <= maxY && ...
                X(i,j,2) >= minCb && X(i,j,2) <= maxCb && ...
                X(i,j,3) >= minCr && X(i,j,3) <= maxCr || ...
                clustered(i,j,4) == faceCluster)
            Z(i,j,:) = Z(i,j,:) - (rightEarMeanValue(1,1,:) - meanValue(1,1,:));
        end
    end
end
%}

disp('Start post-processing(NL-means filter).');
toc

% Post-Processing(NL-means filter)
sigma_NLmeans = 10;
Zf_NLmeans(:,:,1) = NLmeansfilter((Z(:,:,1)),5,2,sigma_NLmeans);
Zf_NLmeans(:,:,2) = NLmeansfilter((Z(:,:,2)),5,2,sigma_NLmeans);
Zf_NLmeans(:,:,3) = NLmeansfilter((Z(:,:,3)),5,2,sigma_NLmeans);
%Zf_NLmeans(:,:,2:3) = Z(:,:,2:3);
%sigma_NLmeans = 20;
%Zf_NLmeans(1:eyeMinH,:,1) = NLmeansfilter(Z(1:eyeMinH,:,1),5,2,sigma_NLmeans);

%{
Zf_NLmeans(mouthRect(1,3)+2:mouthRect(1,4)-1,mouthRect(1,1)+2:mouthRect(1,2)-1,:) = ...
    ZmouthSave;
Zf_NLmeans(noseRect(1,3)+2:noseRect(1,4)-1,noseRect(1,1)+2:noseRect(1,2)-1,:) = ...
    ZnoseSave;
Zf_NLmeans(leftEyeRect(1,3)+2:leftEyeRect(1,4)-1,leftEyeRect(1,1)+2:leftEyeRect(1,2)-1,:) = ...
    ZleftEyeSave;
Zf_NLmeans(rightEyeRect(1,3)+2:rightEyeRect(1,4)-1,rightEyeRect(1,1)+2:rightEyeRect(1,2)-1,:) = ...
    ZrightEyeSave;
%}

Z_ori_NLmeans(:,:,1) = NLmeansfilter((Z_ori(:,:,1)),5,2,sigma_NLmeans);
Z_ori_NLmeans(:,:,2:3) = Z_ori(:,:,2:3);

disp('Start post-processing(guided filter).');
toc

% Post-Processing(guided filter by synthesis, edge-preserving smoothing)
r = 2; % try r=2, 4, or 8
eps = 0.1^2; % try eps=0.1^2, 0.2^2, 0.4^2
Zh_synthesisGuided(:,:,1) = guidedfilter(double(Z(:,:,1))/255, double(Z(:,:,1))/255, r, eps);
Zh_synthesisGuided(:,:,1) = Zh_synthesisGuided(:,:,1) * 255;
Zh_synthesisGuided(:,:,2) = guidedfilter(double(Z(:,:,2))/255, double(Z(:,:,2))/255, r, eps);
Zh_synthesisGuided(:,:,2) = Zh_synthesisGuided(:,:,2) * 255;
Zh_synthesisGuided(:,:,3) = guidedfilter(double(Z(:,:,3))/255, double(Z(:,:,3))/255, r, eps);
Zh_synthesisGuided(:,:,3) = Zh_synthesisGuided(:,:,3) * 255;
%Zh_synthesisGuided(:,:,2:3) = Z(:,:,2:3);

%{
Zh_synthesisGuided(rightEyeRect(1,3)+2:rightEyeRect(1,4)-1,rightEyeRect(1,1)+2:rightEyeRect(1,2)-1,:) = ...
    ZrightEyeSave;
Zh_synthesisGuided(leftEyeRect(1,3)+2:leftEyeRect(1,4)-1,leftEyeRect(1,1)+2:leftEyeRect(1,2)-1,:) = ...
    ZleftEyeSave;
Zh_synthesisGuided(noseRect(1,3)+2:noseRect(1,4)-1,noseRect(1,1)+2:noseRect(1,2)-1,:) = ...
    ZnoseSave;
Zh_synthesisGuided(mouthRect(1,3)+2:mouthRect(1,4)-1,mouthRect(1,1)+2:mouthRect(1,2)-1,:) = ...
    ZmouthSave;
%}

Z_ori_synthesisGuided(:,:,1) = guidedfilter(double(Z_ori(:,:,1))/255, double(Z_ori(:,:,1))/255, r, eps);
Z_ori_synthesisGuided(:,:,1) = Z_ori_synthesisGuided(:,:,1) * 255;
Z_ori_synthesisGuided(:,:,2:3) = Z_ori(:,:,2:3);

% Show/Save and compare photo(X) and cartoon synthesis(Z)
h = figure;
set(h,'name','Comparison photo, cartoon drawn by artist and cartoon synthesis');
subplot(3,3,1),imshow(ycbcr2rgb(uint8(X))),title('Photo');
subplot(3,3,4),imshow(ycbcr2rgb(uint8(Z_ori))),title('Original cartoon synthesis');
subplot(3,3,5),imshow(ycbcr2rgb(uint8(Z_ori_NLmeans))),title('NL-means filtered');
subplot(3,3,6),imshow(ycbcr2rgb(uint8(Z_ori_synthesisGuided))),title('Guided filtered by synthesis');
subplot(3,3,7),imshow(uint8(Z)),title('Replaced cartoon synthesis');
subplot(3,3,8),imshow(uint8(Zf_NLmeans)),title('NL-means filtered');
subplot(3,3,9),imshow(uint8(Zh_synthesisGuided)),title('Guided filtered by synthesis');
truesize;
% Save
im = frame2im(getframe(h));
imwrite(im,['./synthesis_lasso/',name,'_photoCartoonComparison',ext],ext(2:length(ext)));

% Show/Save final cartoon synthesis
% Show
%figure,imshow(ycbcr2rgb(uint8(Z))),title('Cartoon synthsis');
% Save
imwrite(ycbcr2rgb(uint8(X)),['./synthesis_lasso/',name,'_cutPhoto',ext],ext(2:length(ext)));
imwrite(uint8(Z),['./synthesis_lasso/',name,'_cartoonSynthsis_ori',ext],ext(2:length(ext)));
imwrite(uint8(Zf_NLmeans),['./synthesis_lasso/',name,'_cartoonSynthsis_NLmeans',ext],ext(2:length(ext)));

disp('Finish synthesizing cartoon.');
toc