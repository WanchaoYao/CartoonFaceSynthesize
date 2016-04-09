addpath(genpath(pwd));

load './landmark.mat'

bin = 8;
Angle = 360;
L = 3;
pairCount = 100;

% Save noses' images
for i = 1:pairCount
    % Left -> right
    noseRect(1,1) = min(noses{i}(1,:));
    noseRect(1,2) = max(noses{i}(1,:));
    noseRect(1,1) = noseRect(1,1)-(max(noses{i}(1,:))-min(noses{i}(1,:)))*0.5;
    noseRect(1,2) = noseRect(1,2)+(max(noses{i}(1,:))-min(noses{i}(1,:)))*0.5;
    % Top -> bottom
    noseRect(1,3) = (noses{i}(2,7)+noses{i}(2,8))/2;
    noseRect(1,4) = max(noses{i}(2,:));
    
    middlePoints{i} = [round((noseRect(1,1)+noseRect(1,2))/2) round((noses{i}(2,7)+noseRect(1,4))/2)];
    % Decrease top and increase bottom
    %height = noseRect(1,4) - noseRect(1,3);
    %noseRect(1,3) = noseRect(1,3) - 0.25*height;
    %noseRect(1,4) = 0.75*height + noseRect(1,4);
    %noseRect = round(noseRect);
    
    %noseImages{i} = myCutCartoons{i}(noseRect(1,3):noseRect(1,4),noseRect(1,1):noseRect(1,2));
    %figure(1), imshow(myCutCartoons{i}),truesize;
    %figure(2), imshow(noseImages{i}),truesize;
    
    if i == 1
        globalNoseRect = noseRect(1,:);
    else
        globalNoseRect(1,1) = min(globalNoseRect(1,1), noseRect(1,1));
        globalNoseRect(1,2) = max(globalNoseRect(1,2), noseRect(1,2));
        globalNoseRect(1,3) = min(globalNoseRect(1,3), noseRect(1,3));
        globalNoseRect(1,4) = max(globalNoseRect(1,4), noseRect(1,4));
    end        
    
    %pause
end

globalNoseRect = round(globalNoseRect(1,:))
globalNoseRectWidth = globalNoseRect(1,2) - globalNoseRect(1,1);
globalNoseRectHeight = globalNoseRect(1,4) - globalNoseRect(1,3);

for i = 1:pairCount
    i
    % Left -> right
    noseRect(1,1) = middlePoints{i}(1,1) - round(globalNoseRectWidth/2);
    noseRect(1,2) = middlePoints{i}(1,1) + round(globalNoseRectWidth/2);
    % Top -> bottom
    noseRect(1,3) = middlePoints{i}(1,2) - round(globalNoseRectHeight/2);
    noseRect(1,4) = middlePoints{i}(1,2) + round(globalNoseRectHeight/2);
    % Decrease top and increase bottom
    %height = noseRect(1,4) - noseRect(1,3);
    %noseRect(1,3) = noseRect(1,3) - 0.25*height;
    %noseRect(1,4) = 0.75*height + noseRect(1,4);
    %noseRect = round(noseRect);
    noseImages{i} = myCutCartoons{i}(noseRect(1,3):noseRect(1,4),noseRect(1,1):noseRect(1,2),:);
    noseROI(1,1) = noseRect(1,3);
    noseROI(2,1) = noseRect(1,4);
    noseROI(3,1) = noseRect(1,1);
    noseROI(4,1) = noseRect(1,2);
    nosePhogs(:,i) = anna_phog(myCutCartoons{i}(:,:,1),bin,Angle,L,noseROI);
    
    %i
    %figure(1), imshow(myCutCartoons{i}),truesize;
    %figure(2), imshow(noseImages{i}),truesize;      
    %pause
end

[noseIndex,noseParameters] = flann_build_index(nosePhogs, struct('algorithm','kdtree','trees',1,'checks',2048,'cb_index',0));

flann_save_index(noseIndex, './nose_index.index');
flann_free_index(noseIndex);

save './nose.mat' globalNoseRectWidth globalNoseRectHeight noseImages nosePhogs noseParameters