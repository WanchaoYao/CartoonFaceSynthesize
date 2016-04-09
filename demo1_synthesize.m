function demo1_synthesize(photoNo)
% demo1_synthesize Synthesizes a odd number face photo in the training set.
% 
%IN:
%   photoNo - Number(must be odd) of the face photo in the training set

if rem(photoNo,2) ~= 1
    disp('Fail in synthesizing cartoon.');
    disp('Number of the face photo in the training set must be odd.');
    return;
end

tic

addpath(genpath(pwd));

% Read input photo(odd number)
X = rgb2ycbcr(imread(['./photo-1-100_cut/',num2str(photoNo),'.jpg.jpg']));

% Load data and index
load './demo1.mat'
index = flann_load_index('./demo1_index.index', trainingSet);

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
            
			Iq_photo = photos{k}(patches{i,j}(1,1):patches{i,j}(2,1),patches{i,j}(3,1):patches{i,j}(4,1),1);
			Iq_cartoon = cartoons{k}(patches{i,j}(1,1):patches{i,j}(2,1),patches{i,j}(3,1):patches{i,j}(4,1),1);
					
            W1(count) = w1(Xp{m,n},Iq_photo,sigma1);

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
imwrite(im,['./synthesis/',num2str(photoNo),'_patchesComparison.jpg'],'jpg');

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
size(Z)
size(X)
Z(:,:,2:3) = X(:,:,2:3);

disp('Start post-processing(NL-means filter).');
toc

% Post-Processing(NL-means filter)
sigma_NLmeans = 10;
Zf_NLmeans(:,:,1) = NLmeansfilter((Z(:,:,1)),21,7,sigma_NLmeans);
Zf_NLmeans(:,:,2:3) = Z(:,:,2:3);

disp('Start post-processing(guided filter).');
toc

% Post-Processing(guided filter by synthesis, edge-preserving smoothing)
r = 4; % try r=2, 4, or 8
eps = 0.1^2; % try eps=0.1^2, 0.2^2, 0.4^2
Zh_synthesisGuided(:,:,1) = guidedfilter(double(Z(:,:,1))/255, double(Z(:,:,1))/255, r, eps);
Zh_synthesisGuided(:,:,1) = Zh_synthesisGuided(:,:,1) * 255;
Zh_synthesisGuided(:,:,2:3) = Z(:,:,2:3);

% Post-Processing(guided filter by photo, edge-preserving smoothing)
r = 4; % try r=2, 4, or 8
eps = 0.1^2; % try eps=0.1^2, 0.2^2, 0.4^2
Zh_photoGuided(:,:,1) = guidedfilter(double(X(:,:,1))/255, double(Z(:,:,1))/255, r, eps);
Zh_photoGuided(:,:,1) = Zh_photoGuided(:,:,1) * 255;
Zh_photoGuided(:,:,2:3) = Z(:,:,2:3);

% Show/Save and compare photo(X),cartoon drawn by artist(Y) and cartoon synthesis(Z)
Y = imread(['./cartoon-1-100_cut/',num2str(photoNo),'.jpg.jpg']);
% Show
h = figure;
set(h,'name','Comparison photo, cartoon drawn by artist and cartoon synthesis');
subplot(2,3,1),imshow(ycbcr2rgb(uint8(X))),title('Photo');
subplot(2,3,2),imshow(Y),title('Cartoon drawn by artist');
subplot(2,3,3),imshow(ycbcr2rgb(uint8(Z))),title('Original cartoon synthesis');
subplot(2,3,4),imshow(ycbcr2rgb(uint8(Zf_NLmeans))),title('NL-means filtered');
subplot(2,3,5),imshow(ycbcr2rgb(uint8(Zh_synthesisGuided))),title('Guided filtered by synthesis');
subplot(2,3,6),imshow(ycbcr2rgb(uint8(Zh_photoGuided))),title('Guided filtered by photo');
truesize;
% Save
im = frame2im(getframe(h));
imwrite(im,['./synthesis/',num2str(photoNo),'_photoCartoonComparison.jpg'],'jpg');

% Show/Save final cartoon synthesis
% Show
%figure,imshow(ycbcr2rgb(uint8(Z))),title('Cartoon synthsis');
% Save
imwrite(ycbcr2rgb(uint8(Z)),['./synthesis/',num2str(photoNo),'_cartoonSynthsis.jpg'],'jpg');

disp('Finish synthesizing cartoon.');
toc