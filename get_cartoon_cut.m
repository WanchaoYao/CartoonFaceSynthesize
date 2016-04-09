clear all;
clc;

ims = dir('./additional_cartoon/*.jpg');
for i = 1:length(ims),
    close all;
    [X faceProfile leftEye rightEye leftBrow rightBrow nose mouth angle] = ...
        pre_process(['./additional_cartoon/' ims(i).name], true, true);
    
    if isempty(X)
        continue;
    end
    
    imwrite((uint8(X)),['./additional_cartoon_cut/',ims(i).name]);
    
    %disp('press any key to continue');
    %pause;
end