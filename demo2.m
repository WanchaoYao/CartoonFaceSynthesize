clear all;
clc;

for imNumber = 45:45
    close all;
    imNumber
    photoPath = ['./additional_photo/' num2str(uint8(imNumber)) '.jpeg'];
    demo2_lasso(photoPath);
    %demo2_synthesize(photoPath);
    %pause
end
%{
pause;
for imNumber = 1:45
    close all;
    imNumber
    photoPath = ['./additional_photo/' num2str(uint8(imNumber)) '.jpeg'];
    demo2_lasso(photoPath);
    demo2_synthesize(photoPath);
end
%}

