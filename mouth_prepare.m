addpath(genpath(pwd));

load './landmark.mat'

bin = 8;
Angle = 360;
L = 3;
pairCount = 100;

% Save mouths' images
for i = 1:pairCount
    % Left -> right
    mouthRect(1,1) = min(mouths{i}(1,:));
    mouthRect(1,2) = max(mouths{i}(1,:));
    % Top -> bottom
    mouthRect(1,3) = min(mouths{i}(2,:));
    mouthRect(1,4) = max(mouths{i}(2,:));
    
    middlePoints{i} = [round((mouths{i}(1,4)+mouths{i}(1,10))/2) round((mouths{i}(2,4)+mouths{i}(2,10))/2)];
    % Decrease top and increase bottom
    %height = mouthRect(1,4) - mouthRect(1,3);
    %mouthRect(1,3) = mouthRect(1,3) - 0.25*height;
    %mouthRect(1,4) = 0.75*height + mouthRect(1,4);
    %mouthRect = round(mouthRect);
    
    %mouthImages{i} = myCutCartoons{i}(mouthRect(1,3):mouthRect(1,4),mouthRect(1,1):mouthRect(1,2));
    %figure(1), imshow(myCutCartoons{i}),truesize;
    %figure(2), imshow(mouthImages{i}),truesize;
    
    if i == 1
        globalMouthRect = mouthRect(1,:);
    else
        globalMouthRect(1,1) = min(globalMouthRect(1,1), mouthRect(1,1));
        globalMouthRect(1,2) = max(globalMouthRect(1,2), mouthRect(1,2));
        globalMouthRect(1,3) = min(globalMouthRect(1,3), mouthRect(1,3));
        globalMouthRect(1,4) = max(globalMouthRect(1,4), mouthRect(1,4));
    end        
    
    %pause
end

globalMouthRect = round(globalMouthRect(1,:));
globalMouthRectWidth = globalMouthRect(1,2) - globalMouthRect(1,1);
globalMouthRectHeight = globalMouthRect(1,4) - globalMouthRect(1,3);

for i = 1:pairCount
    i
    % Left -> right
    mouthRect(1,1) = middlePoints{i}(1,1) - round(globalMouthRectWidth/2);
    mouthRect(1,2) = middlePoints{i}(1,1) + round(globalMouthRectWidth/2);
    % Top -> bottom
    mouthRect(1,3) = middlePoints{i}(1,2) - round(globalMouthRectHeight/2);
    mouthRect(1,4) = middlePoints{i}(1,2) + round(globalMouthRectHeight/2);
    % Decrease top and increase bottom
    %height = mouthRect(1,4) - mouthRect(1,3);
    %mouthRect(1,3) = mouthRect(1,3) - 0.25*height;
    %mouthRect(1,4) = 0.75*height + mouthRect(1,4);
    %mouthRect = round(mouthRect);
    mouthImages{i} = myCutCartoons{i}(mouthRect(1,3):mouthRect(1,4),mouthRect(1,1):mouthRect(1,2),:);
    mouthROI(1,1) = mouthRect(1,3);
    mouthROI(2,1) = mouthRect(1,4);
    mouthROI(3,1) = mouthRect(1,1);
    mouthROI(4,1) = mouthRect(1,2);
    mouthPhogs(:,i) = anna_phog(myCutCartoons{i}(:,:,1),bin,Angle,L,mouthROI);
    i
    %figure(1), imshow(myCutCartoons{i}),truesize;
    %figure(2), imshow(mouthImages{i}),truesize;       
    
    %pause
end

[mouthIndex,mouthParameters] = flann_build_index(mouthPhogs, struct('algorithm','kdtree','trees',1,'checks',2048,'cb_index',0));

flann_save_index(mouthIndex, './mouth_index.index');
flann_free_index(mouthIndex);

save './mouth.mat' globalMouthRectWidth globalMouthRectHeight mouthImages mouthPhogs mouthParameters