addpath(genpath(pwd));

load './landmark.mat'

bin = 8;
Angle = 360;
L = 3;
pairCount = 100;

for i = 1:pairCount
    % Left brows
    % Left -> right
    leftBrowRect(1,1) = min(leftBrows{i}(1,:));
    leftBrowRect(1,2) = max(leftBrows{i}(1,:));
    % Top -> bottom
    leftBrowRect(1,3) = min(leftBrows{i}(2,:));
    leftBrowRect(1,4) = max(leftBrows{i}(2,:));
    
    leftMiddlePoints{i} = [round((leftBrowRect(1,1)+leftBrowRect(1,2))/2) round((max(leftBrows{i}(2,:))+min(leftBrows{i}(2,:)))/2)];
    
    if i == 1
        globalLeftBrowRect = leftBrowRect(1,:);
    else
        globalLeftBrowRect(1,1) = min(globalLeftBrowRect(1,1), leftBrowRect(1,1));
        globalLeftBrowRect(1,2) = max(globalLeftBrowRect(1,2), leftBrowRect(1,2));
        globalLeftBrowRect(1,3) = min(globalLeftBrowRect(1,3), leftBrowRect(1,3));
        globalLeftBrowRect(1,4) = max(globalLeftBrowRect(1,4), leftBrowRect(1,4));
    end
    
    % Right brows
    % Right -> right
    rightBrowRect(1,1) = min(rightBrows{i}(1,:));
    rightBrowRect(1,2) = max(rightBrows{i}(1,:));
    % Top -> bottom
    rightBrowRect(1,3) = min(rightBrows{i}(2,:));
    rightBrowRect(1,4) = max(rightBrows{i}(2,:));
    
    rightMiddlePoints{i} = [round((rightBrowRect(1,1)+rightBrowRect(1,2))/2) round((max(rightBrows{i}(2,:))+min(rightBrows{i}(2,:)))/2)];
    
    if i == 1
        globalRightBrowRect = rightBrowRect(1,:);
    else
        globalRightBrowRect(1,1) = min(globalRightBrowRect(1,1), rightBrowRect(1,1));
        globalRightBrowRect(1,2) = max(globalRightBrowRect(1,2), rightBrowRect(1,2));
        globalRightBrowRect(1,3) = min(globalRightBrowRect(1,3), rightBrowRect(1,3));
        globalRightBrowRect(1,4) = max(globalRightBrowRect(1,4), rightBrowRect(1,4));
    end
    
    %pause
end

% Left
globalLeftBrowRect = round(globalLeftBrowRect(1,:));
globalLeftBrowRectWidth =  globalLeftBrowRect(1,2) - globalLeftBrowRect(1,1);
globalLeftBrowRectHeight = globalLeftBrowRect(1,4) - globalLeftBrowRect(1,3);

% Right
globalRightBrowRect = round(globalRightBrowRect(1,:));
globalRightBrowRectWidth = globalRightBrowRect(1,2) - globalRightBrowRect(1,1);
globalRightBrowRectHeight = globalRightBrowRect(1,4) - globalRightBrowRect(1,3);

for i = 1:pairCount
    i
    
    % Left
    % Left -> right
    leftBrowRect(1,1) = leftMiddlePoints{i}(1,1) - round(globalLeftBrowRectWidth/2);
    leftBrowRect(1,2) = leftMiddlePoints{i}(1,1) + round(globalLeftBrowRectWidth/2);
    % Top -> bottom
    leftBrowRect(1,3) = leftMiddlePoints{i}(1,2) - round(globalLeftBrowRectHeight/2);
    leftBrowRect(1,4) = leftMiddlePoints{i}(1,2) + round(globalLeftBrowRectHeight/2);

    leftBrowImages{i} = myCutCartoons{i}(leftBrowRect(1,3):leftBrowRect(1,4),leftBrowRect(1,1):leftBrowRect(1,2),:);
    leftBrowROI(1,1) = leftBrowRect(1,3);
    leftBrowROI(2,1) = leftBrowRect(1,4);
    leftBrowROI(3,1) = leftBrowRect(1,1);
    leftBrowROI(4,1) = leftBrowRect(1,2);
    leftBrowPhogs(:,i) = anna_phog(myCutCartoons{i}(:,:,1),bin,Angle,L,leftBrowROI);
 
    % Right
    % Right -> right
    rightBrowRect(1,1) = rightMiddlePoints{i}(1,1) - round(globalRightBrowRectWidth/2);
    rightBrowRect(1,2) = rightMiddlePoints{i}(1,1) + round(globalRightBrowRectWidth/2);
    % Top -> bottom
    rightBrowRect(1,3) = rightMiddlePoints{i}(1,2) - round(globalRightBrowRectHeight/2);
    rightBrowRect(1,4) = rightMiddlePoints{i}(1,2) + round(globalRightBrowRectHeight/2);

    rightBrowImages{i} = myCutCartoons{i}(rightBrowRect(1,3):rightBrowRect(1,4),rightBrowRect(1,1):rightBrowRect(1,2),:);
    rightBrowROI(1,1) = rightBrowRect(1,3);
    rightBrowROI(2,1) = rightBrowRect(1,4);
    rightBrowROI(3,1) = rightBrowRect(1,1);
    rightBrowROI(4,1) = rightBrowRect(1,2);
    rightBrowPhogs(:,i) = anna_phog(myCutCartoons{i}(:,:,1),bin,Angle,L,rightBrowROI);

    %i
    %figure(1), imshow(myCutCartoons{i}),truesize;
    %figure(2), imshow(leftBrowImages{i}),truesize;    
    %figure(3), imshow(rightBrowImages{i}),truesize; 
    %pause
end

% Left
[leftBrowIndex,leftBrowParameters] = flann_build_index(leftBrowPhogs, struct('algorithm','kdtree','trees',1,'checks',2048,'cb_index',0));

flann_save_index(leftBrowIndex, './leftBrow_index.index');
flann_free_index(leftBrowIndex);

save './leftBrow.mat' globalLeftBrowRectWidth globalLeftBrowRectHeight leftBrowImages leftBrowPhogs leftBrowParameters

% Right
[rightBrowIndex,rightBrowParameters] = flann_build_index(rightBrowPhogs, struct('algorithm','kdtree','trees',1,'checks',2048,'cb_index',0));

flann_save_index(rightBrowIndex, './rightBrow_index.index');
flann_free_index(rightBrowIndex);

save './rightBrow.mat' globalRightBrowRectWidth globalRightBrowRectHeight rightBrowImages rightBrowPhogs rightBrowParameters