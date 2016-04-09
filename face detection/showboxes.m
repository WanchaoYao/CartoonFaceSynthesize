function [faceProfile leftEye rightEye leftBrow rightBrow nose mouth] = showboxes(im, boxes, posemap)

% showboxes(im, boxes)
% Draw boxes on top of image.

b = boxes;
%partsize = b.xy(1,3)-b.xy(1,1)+1;
%tx = (min(b.xy(:,1)) + max(b.xy(:,3)))/2;
%ty = min(b.xy(:,2)) - partsize/2;
%text(tx,ty, num2str(posemap(b.c)),'fontsize',18,'color','c');

if size(b.xy,1) ~= 68
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

for i = size(b.xy,1):-1:1;
    x1 = b.xy(i,1);
    y1 = b.xy(i,2);
    x2 = b.xy(i,3);
    y2 = b.xy(i,4);
    %line([x1 x1 x2 x2 x1]', [y1 y2 y2 y1 y1]', 'color', 'b', 'linewidth', 1);

    %i
    %(x1+x2)/2,(y1+y2)/2

    if i <= 68 && i >=52
        tmpFaceProfile(1,68 - i + 1) = (x1+x2)/2;
        tmpFaceProfile(2,68 - i + 1) = (y1+y2)/2;
        %plot((x1+x2)/2,(y1+y2)/2,'r.','markersize',15);
    elseif i <= 51 && i >=32
        mouth(1, i-31) = (x1+x2)/2;
        mouth(2, i-31) = (y1+y2)/2;
        %plot((x1+x2)/2,(y1+y2)/2,'g.','markersize',15);
    elseif i <= 31 && i >=27
        rightBrow(1, i-26) = (x1+x2)/2;
        rightBrow(2, i-26) = (y1+y2)/2;
        %plot((x1+x2)/2,(y1+y2)/2,'b.','markersize',15);
    elseif i <= 26 && i >=21
        rightEye(1, i-20) = (x1+x2)/2;
        rightEye(2, i-20) = (y1+y2)/2;
        %plot((x1+x2)/2,(y1+y2)/2,'c.','markersize',15);
    elseif i <= 20 && i >=16
        leftBrow(1, i-15) = (x1+x2)/2;
        leftBrow(2, i-15) = (y1+y2)/2;
        %plot((x1+x2)/2,(y1+y2)/2,'m.','markersize',15);
    elseif i <= 15 && i >=10
        leftEye(1, i-9) = (x1+x2)/2;
        leftEye(2, i-9) = (y1+y2)/2;
        %plot((x1+x2)/2,(y1+y2)/2,'y.','markersize',15);
    elseif i <= 9 && i >=1
        nose(1, i) = (x1+x2)/2;
        nose(2, i) = (y1+y2)/2;
        %plot((x1+x2)/2,(y1+y2)/2,'k.','markersize',15);
    end

    %pause;
end

for j = 9:16,
    faceProfile(:,j-8) = tmpFaceProfile(:,j);
end
faceProfile(:,9) = tmpFaceProfile(:,17);
for j = 8:-1:1,
    faceProfile(:,18-j) = tmpFaceProfile(:,j);
end
