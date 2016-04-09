function [im faceProfile leftEye rightEye leftBrow rightBrow nose mouth angle] = my_face_detect(photoPath, showImage, saveImage)

% load and visualize model
% Pre-trained model with 146 parts. Works best for faces larger than 80*80
% load face_p146_small.mat

% % Pre-trained model with 99 parts. Works best for faces larger than 150*150
% load face_p99.mat

% % Pre-trained model with 1050 parts. Give best performance on localization, but very slow
load multipie_independent.mat

%{
disp('Model visualization');
visualizemodel(model,1:13);
disp('press any key to continue');
pause;
%}

% 5 levels for each octave
model.interval = 5;
% set up the threshold
model.thresh = min(-0.65, model.thresh);

% define the mapping from view-specific mixture id to viewpoint
if length(model.components)==13 
    posemap = 90:-15:-90;
elseif length(model.components)==18
    posemap = [90:-15:15 0 0 0 0 0 0 -15:-15:-90];
else
    error('Can not recognize this model');
end

im = imread(photoPath);
%%%%% TMP %%%%%
%im = imresize(im, 0.5);
%clf; imagesc(im); axis image; axis off; drawnow;

bs = detect(im, model, model.thresh);
bs = clipboxes(im, bs);
bs = nms_face(bs,0.3);

% show highest scoring one
%{
if size(bs) == 0
    im = imresize(im, 2);
    if size(im,1) <= 400 && size(im,2) <= 400
        bs = detect(im, model, model.thresh);
        bs = clipboxes(im, bs);
        bs = nms_face(bs,0.3);
    end
end
%}

for i = 1:3
    i;
    if size(bs) == 0
        if size(im,1) >= 600 && size(im,2) >= 600
            break
        end
        im = imresize(im, 2);
        bs = detect(im, model, model.thresh);
        bs = clipboxes(im, bs);
        bs = nms_face(bs,0.3);
    else
        break;
    end
end

if size(bs) == 0
    disp([photoPath ' - No face! ! ! ! ! ! ! ! ! ! !']);
    im = [];
    faceProfile = [];
    leftEye = [];
    rightEye = [];
    leftBrow = [];
    rightBrow = [];
    nose = [];
    mouth = [];
    angle = [];
    return;
end

[faceProfile leftEye rightEye leftBrow rightBrow nose mouth] = showboxes(im, bs(1),posemap);
if size(bs) == 0
    disp([photoPath ' - No face! ! ! ! ! ! ! ! ! ! !']);
    im = [];
    return;
end

b = bs(1);
angle = posemap(b.c);

if showImage
    h = figure;
    title('Highest scoring detection');
    imagesc(im);
    hold on;
    axis image;
    axis off;
    
    partsize = b.xy(1,3)-b.xy(1,1)+1;
    tx = (min(b.xy(:,1)) + max(b.xy(:,3)))/2;
    ty = min(b.xy(:,2)) - partsize/2;
    text(tx,ty, num2str(angle),'fontsize',18,'color','c');
    
    for i = 1:length(leftEye)
        plot(leftEye(1,i),leftEye(2,i),'g.','markersize',15);
        plot(rightEye(1,i),rightEye(2,i),'b.','markersize',15);
    end

    for i = 1:length(leftBrow)
        plot(leftBrow(1,i),leftBrow(2,i),'c.','markersize',15);
        plot(rightBrow(1,i),rightBrow(2,i),'m.','markersize',15);
    end

    for i = 1:length(nose)
        plot(nose(1,i),nose(2,i),'y.','markersize',15);
        %pause
    end

    for i = 1:length(mouth)
        plot(mouth(1,i),mouth(2,i),'k.','markersize',15);
    end

    x = faceProfile(1,1):1:faceProfile(1,17);
    for i = 1:length(faceProfile)
        plot(faceProfile(1,i),faceProfile(2,i),'r.','markersize',15);
    end

    func = polyfit(faceProfile(1,:), faceProfile(2,:), 2);  %进行拟合，c为2次拟合后的系数
    result = polyval(func, x, 1);  %拟合后，每一个横坐标对应的值即为d
    plot(x, result, 'w', 'markersize', 10);   %拟合后的曲线

    values = spcrv(faceProfile,2);
    plot(values(1,:),values(2,:), 'g');

    bottom = faceProfile(1,4):1:faceProfile(1,14);
    funcBottom = polyfit(faceProfile(1,4:14), faceProfile(2,4:14), 2);  %进行拟合，c为2次拟合后的系数
    resultBottom = polyval(funcBottom, bottom, 1);  %拟合后，每一个横坐标对应的值即为d
    plot(bottom, resultBottom, 'b', 'markersize', 10);  %拟合后的曲线    

    drawnow;

    % show all
    %figure,showboxes(im, bs,posemap),title('All detections above the threshold');

    %figure, imshow(im(:,:,1));
    
    if saveImage
        figureShot = frame2im(getframe(h));
        % Get name and extension of the input photo
        [~, name, ext] = fileparts(photoPath);
        saveas(h,['./additional_cartoon_face_alignment/',name],'jpg');
    end
end
