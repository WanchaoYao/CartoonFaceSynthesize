function [X faceProfile leftEye rightEye leftBrow rightBrow nose mouth angle] = pre_process(photoPath, showImage, saveImage)

addpath(genpath(pwd));

xlength = 96;
ylength = 128;

% Face alignment
[X faceProfile leftEye rightEye leftBrow rightBrow nose mouth angle] = my_face_detect(photoPath, showImage, saveImage);

if isempty(X)
    return;
end

% Pre-processing (crop input photo)
leftEyePos = [mean(leftEye(1,:)) mean(leftEye(2,:))];
rightEyePos = [mean(rightEye(1,:)) mean(rightEye(2,:))];
eyeSpacing = rightEyePos(1,1) - leftEyePos(1,1);
rect(1,1) = leftEyePos(1,1) - eyeSpacing;
if rect(1,1) < 0
    photoPath;
    rect(1,1) = 0;
end
rect(1,2) = leftEyePos(1,2) - 2 * eyeSpacing;
if rect(1,2) < 0
    photoPath;
    rect(1,2) = 0;
end
rect(1,3) = 3*eyeSpacing;
if (rect(1,3)+rect(1,1))  > size(X,2)
    photoPath;
    rect(1,3) = size(X,2) - rect(1,1);
end
rect(1,4) = 4*eyeSpacing;
if (rect(1,4)+rect(1,2))  > size(X,1)
    photoPath;
    rect(1,4) = size(X,1) - rect(1,2);
end

oriYlength = rect(1,4);
oriXlength = rect(1,3);

X = imcrop(X, rect);
X = imresize(X, [ylength xlength]);

faceProfile(1,:) = (faceProfile(1,:) - rect(1,1)) * (xlength / oriXlength);
faceProfile(2,:) = (faceProfile(2,:) - rect(1,2)) * (ylength / oriYlength);

leftEye(1,:) = (leftEye(1,:) - rect(1,1)) * (xlength / oriXlength);
leftEye(2,:) = (leftEye(2,:) - rect(1,2)) * (ylength / oriYlength);

rightEye(1,:) = (rightEye(1,:) - rect(1,1)) * (xlength / oriXlength);
rightEye(2,:) = (rightEye(2,:) - rect(1,2)) * (ylength / oriYlength);

leftBrow(1,:) = (leftBrow(1,:) - rect(1,1)) * (xlength / oriXlength);
leftBrow(2,:) = (leftBrow(2,:) - rect(1,2)) * (ylength / oriYlength);

rightBrow(1,:) = (rightBrow(1,:) - rect(1,1)) * (xlength / oriXlength);
rightBrow(2,:) = (rightBrow(2,:) - rect(1,2)) * (ylength / oriYlength);

nose(1,:) = (nose(1,:) - rect(1,1)) * (xlength / oriXlength);
nose(2,:) = (nose(2,:) - rect(1,2)) * (ylength / oriYlength);

mouth(1,:) = (mouth(1,:) - rect(1,1)) * (xlength / oriXlength);
mouth(2,:) = (mouth(2,:) - rect(1,2)) * (ylength / oriYlength);

if showImage
    figure,imshow(X);
end