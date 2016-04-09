function demo2_prepare()

tic

addpath(genpath(pwd));

xlength = 96;
ylength = 128;
pairCount = 1000;

patchSize = 16;
overlappingDis = 11;
spacing = patchSize-overlappingDis;

sigma1 = 60.0;
sigma2 = 0.02;
sigma3 = 10.0;
sigma4 = 50.0;

% Opening the matlabpool
if matlabpool('size')==0
	matlabpool open 4
end

% Read training cartoons
cartoons = cell(pairCount,1);
for i = 1:pairCount
    if i <= 100
        cartoons{i} = double(rgb2ycbcr(imread(['./cartoon-1-100_cut/',num2str(i),'.jpg.jpg'])));
    else
        cartoons{i} = double(rgb2ycbcr(imread(['./additional_cartoon_cut/',num2str(i),'.jpg'])));
    end
end

disp('Finish reading cartoons.');
toc

% Calculate coordinate of all patches
patchRowCount = floor((ylength-overlappingDis)/(patchSize-overlappingDis));
patchColCount = floor((xlength-overlappingDis)/(patchSize-overlappingDis));
patches = cell(patchRowCount, patchColCount);
parfor i = 1:patchRowCount
	for j = 1:patchColCount
		% patches{i,j} = [ytop;ybottom;xleft;xright]
		patches{i,j} = [spacing*(i-1)+1;spacing*(i-1)+patchSize;spacing*(j-1)+1;spacing*(j-1)+patchSize];
	end
end

% Compute PHOG over all patches of cartoons and input photo
bin = 8;
Angle = 360;
L = 3;
phogs = cell(pairCount,1);
parfor k = 1:pairCount
	phogs{k} = my_phog(cartoons{k}(:,:,1),bin,Angle,L,patches);
    k
end

matlabpool close

% Construct trainings set
trainingSet = zeros(size(phogs{1}{1,1},1),pairCount*patchRowCount*patchColCount);
for k = 1:pairCount
    for i = 1:patchRowCount
        for j = 1:patchColCount
            trainingSet(:,(k-1)*patchRowCount*patchColCount+(i-1)*patchColCount+j) = phogs{k}{i,j};
        end
    end
end

disp('Start building index.');
toc

%build_params.algorithm = 'autotuned';
%build_params.target_precision = 0.90;
%build_params.build_weight = 1;
%build_params.memory_weight = 1;
[index,parameters] = flann_build_index(trainingSet, struct('algorithm','kdtree','trees',1,'checks',2048,'cb_index',0));
parameters

save -v7.3 './demo2.mat'

flann_save_index(index, './demo2_index.index');
flann_free_index(index);

disp('Finish preparation.');
toc
