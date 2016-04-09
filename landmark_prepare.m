addpath(genpath(pwd));

xlength = 96;
ylength = 128;
pairCount = 100;

% Read training cartoons
myCutCartoons = cell(pairCount,1);
faceProfiles = cell(pairCount,1);
leftEyes = cell(pairCount,1);
rightEyes = cell(pairCount,1);
leftBrows = cell(pairCount,1);
rightBrows = cell(pairCount,1);
noses = cell(pairCount,1);
mouths = cell(pairCount,1);
for i = 1:pairCount
    photoPath = ['./cartoon-1-100/',num2str(i),'.jpg'];
	[X faceProfile leftEye rightEye leftBrow rightBrow nose mouth] = pre_process(photoPath, false);
    
    %figure, imshow(X)
    
    myCutCartoons{i} = X;
    faceProfiles{i} = faceProfile;
    leftEyes{i} = leftEye;
    rightEyes{i} = rightEye;
    leftBrows{i} = leftBrow;
    rightBrows{i} = rightBrow;
    noses{i} = nose;
    mouths{i} = mouth;
end

save './landmark.mat' myCutCartoons faceProfiles leftEyes rightEyes leftBrows rightBrows noses mouths;