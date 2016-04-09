function demo2_synthesize(photoPath)
% demo2_synthesize Synthesizes a face cartoon from a face photo.
% 
%IN:
%   photoPath - Path of the face photo

tic

addpath(genpath(pwd));

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

[clustered,C] = pixel_cluster(X,3);

X = rgb2ycbcr(X);

% Get name and extension of the input photo
[~, name, ext] = fileparts(photoPath);

% Compute PHOG over all patches the input photo
D = my_phog(X(:,:,1),bin,Angle,L,patches);

% Construct test set
testSet = zeros(size(phogs{1}{1,1},1),patchRowCount*patchColCount);
for i = 1:patchRowCount
    for j = 1:patchColCount
        testSet(:,(i-1)*patchColCount+j) = D{i,j}; 
    end
end

disp('Start searching nearest neighbors.');
toc

neighborsCount = 100;
[result, dists] = flann_search(index, testSet, neighborsCount, parameters);
flann_free_index(index);

disp('Start synthesizing cartoon.');
toc

% Synthesize all cartoon patches
Zp = cell(patchRowCount, patchColCount);
Xp = cell(patchRowCount, patchColCount);
W1 = zeros(neighborsCount,1);
W2 = zeros(neighborsCount,1);
W3 = zeros(neighborsCount,1);
W4 = zeros(neighborsCount,1);
W = zeros(neighborsCount,1);
patchIndexes = zeros(neighborsCount,1);
for m = 1:patchRowCount
	for n = 1:patchColCount
		Zp{m,n} = 0;
		Xp{m,n} = X(patches{m,n}(1,1):patches{m,n}(2,1),patches{m,n}(3,1):patches{m,n}(4,1),1);
        
        patchIndexes = result(:,(m-1)*patchColCount+n);
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
                pause
            end
            
			%Iq_photo = photos{k}(patches{i,j}(1,1):patches{i,j}(2,1),patches{i,j}(3,1):patches{i,j}(4,1),1);
			Iq_cartoon = cartoons{k}(patches{i,j}(1,1):patches{i,j}(2,1),patches{i,j}(3,1):patches{i,j}(4,1),1);
					
            %W1(count) = w1(Xp{m,n},Iq_photo,sigma1);
            W1(count) = 1;

            W2(count) = w2(D{m,n},phogs{k}{i,j},sigma2);

            if m==1
                Zpt = [];
            else
                Zpt = Zp{m-1,n};
            end

            if n==1
                Zpl = [];
            else
                Zpl = Zp{m,n-1};
            end

            W3(count) = w3(Zpt,Zpl,Iq_cartoon,patchSize,overlappingDis,sigma3);

            W4(count) = w4([patches{m,n}(1,1);patches{m,n}(3,1)],[patches{i,j}(1,1);patches{i,j}(3,1)],sigma4);

            W(count) = W1(count)*W2(count)*W3(count)*W4(count);

            Zp{m,n} = Zp{m,n} + W(count)*(Iq_cartoon-mean(mean(Iq_cartoon)));
        end

		Zp{m,n} = Zp{m,n} / sum(W);
		Zp{m,n} = Zp{m,n} + mean(mean(Xp{m,n}));

        for ii = 1:patchSize
            for jj = 1:patchSize
                if Zp{m,n}(ii,jj)>255
					Zp{m,n}(ii,jj) = 255;
                end
                if Zp{m,n}(ii,jj)<0
					Zp{m,n}(ii,jj) = 0;
                end
            end
        end
       
        %fprintf('Finish synthesizing patch{%d,%d}.\n',m,n);
        %toc
	end
end

% Show/Save and compare photo patches and synthesized cartoon patches
for m = 1:patchRowCount
	for n = 1:patchColCount
		A(patchSize*(m-1)+1:patchSize*m,patchSize*(n-1)+1:patchSize*n) = Xp{m,n};
		B(patchSize*(m-1)+1:patchSize*m,patchSize*(n-1)+1:patchSize*n) = Zp{m,n};
	end
end
% Show
h = figure;
set(h,'name','Comparison between photo patches and synthesized cartoon patches');
subplot(1,2,1),imshow(uint8(A)),title('Patches of photo');
subplot(1,2,2),imshow(uint8(B)),title('Patches of synthesized cartoon');
truesize;
% Save
im = frame2im(getframe(h));
imwrite(im,['./synthesis_weight_function/',name,'_patchesComparison',ext],ext(2:length(ext)));

disp('Start putting together all patches.');
toc

% Put together all patches to synthesize cartoon Z
Z = zeros(ylength,xlength,3);

overlapCount = zeros(ylength,xlength);
for m = 1:patchRowCount
	for n = 1:patchColCount
		Z(patches{m,n}(1,1):patches{m,n}(2,1),patches{m,n}(3,1):patches{m,n}(4,1),1) = ...
            Z(patches{m,n}(1,1):patches{m,n}(2,1),patches{m,n}(3,1):patches{m,n}(4,1),1) + Zp{m,n};
        overlapCount(patches{m,n}(1,1):patches{m,n}(2,1),patches{m,n}(3,1):patches{m,n}(4,1)) = ...
            overlapCount(patches{m,n}(1,1):patches{m,n}(2,1),patches{m,n}(3,1):patches{m,n}(4,1)) + 1;        
	end
end

for i = 1:ylength
	for j = 1:xlength
        if overlapCount(i,j) ~= 0
    		Z(i,j,1) = Z(i,j,1)/overlapCount(i,j);    
        end
	end
end

% Color-Preserving Colorization
Z(:,:,2:3) = X(:,:,2:3);

Z_ori = Z;

% Replace background
%{
background = rgb2ycbcr(clustered(1,1,1:3));
for i = 1:ylength-2
    for j = 1:xlength
        if clustered(i,j,4) == clustered(1,1,4)
            Z(i,j,:) = background;
        end
    end
end
%}

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
faceMiddleH = floor(min(yv));
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
        if inpolygon(j,i,xv,yv) == 1 && faceCluster == clustered(i,j,4)
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

faceMinH = faceMiddleH-(faceMaxH-faceMiddleH);
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
            if inpolygon(j,i,xv,yv) == 1 && ...
                    X(i,j,1) >= minY && X(i,j,1) <= maxY && ...
                    X(i,j,2) >= minCb && X(i,j,2) <= maxCb && ...
                    X(i,j,3) >= minCr && X(i,j,3) <= maxCr
                Z(i,j,:) = meanValue(1,1,:);
                %imMask(i,j) = 1;
            end
        else
            if X(i,j,1) >= minY && X(i,j,1) <= maxY && ...
                    X(i,j,2) >= minCb && X(i,j,2) <= maxCb && ... 
                    X(i,j,3) >= minCr && X(i,j,3) <= maxCr
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
    my_poissonSolverSeamlessCloning(Zmouth(:,:,2),Z(:,:,3),...
    [2 size(Zmouth,2)-2 2 size(Zmouth,1)-2],[mouthRect(1,1)+1 mouthRect(1,3)+1]);
Z(mouthRect(1,3)+2:mouthRect(1,4)-1,mouthRect(1,1)+2:mouthRect(1,2)-1,3) = imRegion;

% Change nose
Z(noseRect(1,3)+2:noseRect(1,4)-1,noseRect(1,1)+2:noseRect(1,2)-1,1) = ...
    my_poissonSolverSeamlessCloning(Znose(:,:,1),Z(:,:,1),...
    [2 size(Znose,2)-2 2 size(Znose,1)-2],[noseRect(1,1)+1 noseRect(1,3)+1]);
Z(noseRect(1,3)+2:noseRect(1,4)-1,noseRect(1,1)+2:noseRect(1,2)-1,2) = ...
    my_poissonSolverSeamlessCloning(Znose(:,:,2),Z(:,:,2),...
    [2 size(Znose,2)-2 2 size(Znose,1)-2],[noseRect(1,1)+1 noseRect(1,3)+1]);
Z(noseRect(1,3)+2:noseRect(1,4)-1,noseRect(1,1)+2:noseRect(1,2)-1,3) = ...
    my_poissonSolverSeamlessCloning(Znose(:,:,2),Z(:,:,3),...
    [2 size(Znose,2)-2 2 size(Znose,1)-2],[noseRect(1,1)+1 noseRect(1,3)+1]);

% Change eyes
Z(leftEyeRect(1,3)+2:leftEyeRect(1,4)-1,leftEyeRect(1,1)+2:leftEyeRect(1,2)-1,1) = ...
    my_poissonSolverSeamlessCloning(ZleftEye(:,:,1),Z(:,:,1),...
    [2 size(ZleftEye,2)-2 2 size(ZleftEye,1)-2],[leftEyeRect(1,1)+1 leftEyeRect(1,3)+1]);
Z(leftEyeRect(1,3)+2:leftEyeRect(1,4)-1,leftEyeRect(1,1)+2:leftEyeRect(1,2)-1,2) = ...
    my_poissonSolverSeamlessCloning(ZleftEye(:,:,2),Z(:,:,2),...
    [2 size(ZleftEye,2)-2 2 size(ZleftEye,1)-2],[leftEyeRect(1,1)+1 leftEyeRect(1,3)+1]);
Z(leftEyeRect(1,3)+2:leftEyeRect(1,4)-1,leftEyeRect(1,1)+2:leftEyeRect(1,2)-1,3) = ...
    my_poissonSolverSeamlessCloning(ZleftEye(:,:,2),Z(:,:,3),...
    [2 size(ZleftEye,2)-2 2 size(ZleftEye,1)-2],[leftEyeRect(1,1)+1 leftEyeRect(1,3)+1]);

Z(rightEyeRect(1,3)+2:rightEyeRect(1,4)-1,rightEyeRect(1,1)+2:rightEyeRect(1,2)-1,1) = ...
    my_poissonSolverSeamlessCloning(ZrightEye(:,:,1),Z(:,:,1),...
    [2 size(ZrightEye,2)-2 2 size(ZrightEye,1)-2],[rightEyeRect(1,1)+1 rightEyeRect(1,3)+1]);
Z(rightEyeRect(1,3)+2:rightEyeRect(1,4)-1,rightEyeRect(1,1)+2:rightEyeRect(1,2)-1,2) = ...
    my_poissonSolverSeamlessCloning(ZrightEye(:,:,2),Z(:,:,2),...
    [2 size(ZrightEye,2)-2 2 size(ZrightEye,1)-2],[rightEyeRect(1,1)+1 rightEyeRect(1,3)+1]);
Z(rightEyeRect(1,3)+2:rightEyeRect(1,4)-1,rightEyeRect(1,1)+2:rightEyeRect(1,2)-1,3) = ...
    my_poissonSolverSeamlessCloning(ZrightEye(:,:,2),Z(:,:,3),...
    [2 size(ZrightEye,2)-2 2 size(ZrightEye,1)-2],[rightEyeRect(1,1)+1 rightEyeRect(1,3)+1]);

% Change brows
Z(leftBrowRect(1,3)+2:leftBrowRect(1,4)-1,leftBrowRect(1,1)+2:leftBrowRect(1,2)-1,1) = ...
    my_poissonSolverSeamlessCloning(ZleftBrow(:,:,1),Z(:,:,1),...
    [2 size(ZleftBrow,2)-2 2 size(ZleftBrow,1)-2],[leftBrowRect(1,1)+1 leftBrowRect(1,3)+1]);
Z(leftBrowRect(1,3)+2:leftBrowRect(1,4)-1,leftBrowRect(1,1)+2:leftBrowRect(1,2)-1,2) = ...
    my_poissonSolverSeamlessCloning(ZleftBrow(:,:,2),Z(:,:,2),...
    [2 size(ZleftBrow,2)-2 2 size(ZleftBrow,1)-2],[leftBrowRect(1,1)+1 leftBrowRect(1,3)+1]);
Z(leftBrowRect(1,3)+2:leftBrowRect(1,4)-1,leftBrowRect(1,1)+2:leftBrowRect(1,2)-1,3) = ...
    my_poissonSolverSeamlessCloning(ZleftBrow(:,:,2),Z(:,:,3),...
    [2 size(ZleftBrow,2)-2 2 size(ZleftBrow,1)-2],[leftBrowRect(1,1)+1 leftBrowRect(1,3)+1]);

Z(rightBrowRect(1,3)+2:rightBrowRect(1,4)-1,rightBrowRect(1,1)+2:rightBrowRect(1,2)-1,1) = ...
    my_poissonSolverSeamlessCloning(ZrightBrow(:,:,1),Z(:,:,1),...
    [2 size(ZrightBrow,2)-2 2 size(ZrightBrow,1)-2],[rightBrowRect(1,1)+1 rightBrowRect(1,3)+1]);
Z(rightBrowRect(1,3)+2:rightBrowRect(1,4)-1,rightBrowRect(1,1)+2:rightBrowRect(1,2)-1,2) = ...
    my_poissonSolverSeamlessCloning(ZrightBrow(:,:,2),Z(:,:,2),...
    [2 size(ZrightBrow,2)-2 2 size(ZrightBrow,1)-2],[rightBrowRect(1,1)+1 rightBrowRect(1,3)+1]);
Z(rightBrowRect(1,3)+2:rightBrowRect(1,4)-1,rightBrowRect(1,1)+2:rightBrowRect(1,2)-1,3) = ...
    my_poissonSolverSeamlessCloning(ZrightBrow(:,:,2),Z(:,:,3),...
    [2 size(ZrightBrow,2)-2 2 size(ZrightBrow,1)-2],[rightBrowRect(1,1)+1 rightBrowRect(1,3)+1]);

disp('Start post-processing(NL-means filter).');
toc

% Post-Processing(NL-means filter)
sigma_NLmeans = 20;
Zf_NLmeans(:,:,1) = NLmeansfilter((Z(:,:,1)),5,2,sigma_NLmeans);
Zf_NLmeans(:,:,2:3) = Z(:,:,2:3);

Z_ori_NLmeans(:,:,1) = NLmeansfilter((Z_ori(:,:,1)),5,2,sigma_NLmeans);
Z_ori_NLmeans(:,:,2:3) = Z_ori(:,:,2:3);

disp('Start post-processing(guided filter).');
toc

% Post-Processing(guided filter by synthesis, edge-preserving smoothing)
r = 2; % try r=2, 4, or 8
eps = 0.1^2; % try eps=0.1^2, 0.2^2, 0.4^2
Zh_synthesisGuided(:,:,1) = guidedfilter(double(Z(:,:,1))/255, double(Z(:,:,1))/255, r, eps);
Zh_synthesisGuided(:,:,1) = Zh_synthesisGuided(:,:,1) * 255;
Zh_synthesisGuided(:,:,2:3) = Z(:,:,2:3);

Z_ori_synthesisGuided(:,:,1) = guidedfilter(double(Z_ori(:,:,1))/255, double(Z_ori(:,:,1))/255, r, eps);
Z_ori_synthesisGuided(:,:,1) = Z_ori_synthesisGuided(:,:,1) * 255;
Z_ori_synthesisGuided(:,:,2:3) = Z_ori(:,:,2:3);

clustered(mouthRect(1,3):mouthRect(1,4),mouthRect(1,1):mouthRect(1,2),1:3) = ycbcr2rgb(uint8(Zmouth));
clustered(noseRect(1,3):noseRect(1,4),noseRect(1,1):noseRect(1,2),1:3) = ycbcr2rgb(uint8(Znose));
clustered(leftEyeRect(1,3):leftEyeRect(1,4),leftEyeRect(1,1):leftEyeRect(1,2),1:3) = ycbcr2rgb(uint8(ZleftEye));
clustered(rightEyeRect(1,3):rightEyeRect(1,4),rightEyeRect(1,1):rightEyeRect(1,2),1:3) = ycbcr2rgb(uint8(ZrightEye));

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
imwrite(im,['./synthesis_weight_function/',name,'_photoCartoonComparison',ext],ext(2:length(ext)));

% Show/Save final cartoon synthesis
% Show
%figure,imshow(ycbcr2rgb(uint8(Z))),title('Cartoon synthsis');
% Save
imwrite(uint8(Zf_NLmeans),['./synthesis_weight_function/',name,'_cartoonSynthsis',ext],ext(2:length(ext)));

disp('Finish synthesizing cartoon.');
toc