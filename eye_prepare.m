addpath(genpath(pwd));

load './landmark.mat'

bin = 8;
Angle = 360;
L = 3;
pairCount = 100;

for i = 1:pairCount
    % Left eyes
    % Left -> right
    leftEyeRect(1,1) = min(leftEyes{i}(1,:));
    leftEyeRect(1,2) = max(leftEyes{i}(1,:));
    leftEyeRect(1,1) = leftEyeRect(1,1)-(max(leftEyes{i}(1,:))-min(leftEyes{i}(1,:)))*0.5;
    leftEyeRect(1,2) = leftEyeRect(1,2)+(max(leftEyes{i}(1,:))-min(leftEyes{i}(1,:)))*0.5;
    % Top -> bottom
    leftEyeRect(1,3) = min(leftEyes{i}(2,:));
    leftEyeRect(1,4) = max(leftEyes{i}(2,:));
    leftEyeRect(1,3) = leftEyeRect(1,3)-(max(leftEyes{i}(2,:))-min(leftEyes{i}(2,:)))*0.5;
    leftEyeRect(1,4) = leftEyeRect(1,4)+(max(leftEyes{i}(2,:))-min(leftEyes{i}(2,:)))*0.75;
    
    leftMiddlePoints{i} = [round((leftEyeRect(1,1)+leftEyeRect(1,2))/2) round((max(leftEyes{i}(2,:))+min(leftEyes{i}(2,:)))/2)];
    
    if i == 1
        globalLeftEyeRect = leftEyeRect(1,:);
    else
        globalLeftEyeRect(1,1) = min(globalLeftEyeRect(1,1), leftEyeRect(1,1));
        globalLeftEyeRect(1,2) = max(globalLeftEyeRect(1,2), leftEyeRect(1,2));
        globalLeftEyeRect(1,3) = min(globalLeftEyeRect(1,3), leftEyeRect(1,3));
        globalLeftEyeRect(1,4) = max(globalLeftEyeRect(1,4), leftEyeRect(1,4));
    end
    
    % Right eyes
    % Right -> right
    rightEyeRect(1,1) = min(rightEyes{i}(1,:));
    rightEyeRect(1,2) = max(rightEyes{i}(1,:));
    rightEyeRect(1,1) = rightEyeRect(1,1)-(max(rightEyes{i}(1,:))-min(rightEyes{i}(1,:)))*0.25;
    rightEyeRect(1,2) = rightEyeRect(1,2)+(max(rightEyes{i}(1,:))-min(rightEyes{i}(1,:)))*0.25;
    % Top -> bottom
    rightEyeRect(1,3) = min(rightEyes{i}(2,:));
    rightEyeRect(1,4) = max(rightEyes{i}(2,:));
    rightEyeRect(1,3) = rightEyeRect(1,3)-(max(rightEyes{i}(2,:))-min(rightEyes{i}(2,:)))*0.25;
    rightEyeRect(1,4) = rightEyeRect(1,4)+(max(rightEyes{i}(2,:))-min(rightEyes{i}(2,:)))*0.75;
    
    rightMiddlePoints{i} = [round((rightEyeRect(1,1)+rightEyeRect(1,2))/2) round((max(rightEyes{i}(2,:))+min(rightEyes{i}(2,:)))/2)];
    
    if i == 1
        globalRightEyeRect = rightEyeRect(1,:);
    else
        globalRightEyeRect(1,1) = min(globalRightEyeRect(1,1), rightEyeRect(1,1));
        globalRightEyeRect(1,2) = max(globalRightEyeRect(1,2), rightEyeRect(1,2));
        globalRightEyeRect(1,3) = min(globalRightEyeRect(1,3), rightEyeRect(1,3));
        globalRightEyeRect(1,4) = max(globalRightEyeRect(1,4), rightEyeRect(1,4));
    end
    
    %pause
end

% Left
globalLeftEyeRect = round(globalLeftEyeRect(1,:))
globalLeftEyeRectWidth = globalLeftEyeRect(1,2) - globalLeftEyeRect(1,1);
globalLeftEyeRectHeight = globalLeftEyeRect(1,4) - globalLeftEyeRect(1,3);

% Right
globalRightEyeRect = round(globalRightEyeRect(1,:))
globalRightEyeRectWidth = globalRightEyeRect(1,2) - globalRightEyeRect(1,1);
globalRightEyeRectHeight = globalRightEyeRect(1,4) - globalRightEyeRect(1,3);

for i = 1:pairCount
    i
    
    % Left
    % Left -> right
    leftEyeRect(1,1) = leftMiddlePoints{i}(1,1) - round(globalLeftEyeRectWidth/2);
    leftEyeRect(1,2) = leftMiddlePoints{i}(1,1) + round(globalLeftEyeRectWidth/2);
    % Top -> bottom
    leftEyeRect(1,3) = leftMiddlePoints{i}(1,2) - round(globalLeftEyeRectHeight/2);
    leftEyeRect(1,4) = leftMiddlePoints{i}(1,2) + round(globalLeftEyeRectHeight/2);

    leftEyeImages{i} = myCutCartoons{i}(leftEyeRect(1,3):leftEyeRect(1,4),leftEyeRect(1,1):leftEyeRect(1,2),:);
    leftEyeROI(1,1) = leftEyeRect(1,3);
    leftEyeROI(2,1) = leftEyeRect(1,4);
    leftEyeROI(3,1) = leftEyeRect(1,1);
    leftEyeROI(4,1) = leftEyeRect(1,2);
    leftEyePhogs(:,i) = anna_phog(myCutCartoons{i}(:,:,1),bin,Angle,L,leftEyeROI);
 
    % Right
    % Right -> right
    rightEyeRect(1,1) = rightMiddlePoints{i}(1,1) - round(globalRightEyeRectWidth/2);
    rightEyeRect(1,2) = rightMiddlePoints{i}(1,1) + round(globalRightEyeRectWidth/2);
    % Top -> bottom
    rightEyeRect(1,3) = rightMiddlePoints{i}(1,2) - round(globalRightEyeRectHeight/2);
    rightEyeRect(1,4) = rightMiddlePoints{i}(1,2) + round(globalRightEyeRectHeight/2);

    rightEyeImages{i} = myCutCartoons{i}(rightEyeRect(1,3):rightEyeRect(1,4),rightEyeRect(1,1):rightEyeRect(1,2),:);
    rightEyeROI(1,1) = rightEyeRect(1,3);
    rightEyeROI(2,1) = rightEyeRect(1,4);
    rightEyeROI(3,1) = rightEyeRect(1,1);
    rightEyeROI(4,1) = rightEyeRect(1,2);
    rightEyePhogs(:,i) = anna_phog(myCutCartoons{i}(:,:,1),bin,Angle,L,rightEyeROI);
    
    %i
    figure(1), imshow(myCutCartoons{i}),truesize;
    figure(2), imshow(leftEyeImages{i}),truesize;    
    figure(3), imshow(rightEyeImages{i}),truesize; 
    pause
end

% Left
[leftEyeIndex,leftEyeParameters] = flann_build_index(leftEyePhogs, struct('algorithm','kdtree','trees',1,'checks',2048,'cb_index',0));

flann_save_index(leftEyeIndex, './leftEye_index.index');
flann_free_index(leftEyeIndex);

save './leftEye.mat' globalLeftEyeRectWidth globalLeftEyeRectHeight leftEyeImages leftEyePhogs leftEyeParameters

% Right
[rightEyeIndex,rightEyeParameters] = flann_build_index(rightEyePhogs, struct('algorithm','kdtree','trees',1,'checks',2048,'cb_index',0));

flann_save_index(rightEyeIndex, './rightEye_index.index');
flann_free_index(rightEyeIndex);

save './rightEye.mat' globalRightEyeRectWidth globalRightEyeRectHeight rightEyeImages rightEyePhogs rightEyeParameters