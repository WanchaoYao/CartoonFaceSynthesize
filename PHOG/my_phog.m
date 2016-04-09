function p = my_phog(I,bin,angle,L,patches)
% my_PHOG Computes Pyramid Histogram of Oriented Gradient over all patches.
% 
% Given an image I, phog computes the Pyramid Histogram of Oriented Gradients
% over L pyramid levels and over all patches
%
%IN:
%   I - Image of size MxN (Color or Gray)
%   bin - Number of bins on the histogram 
%   angle - 180 or 360
%   L - number of pyramid levels
%   patches - All patches (ytop,ybottom,xleft,xright)
%
%OUT:
%   p - pyramid histogram of oriented gradients(PHOG) over all patches

if size(I,3) == 3
    G = rgb2gray(I);
else
    G = I;
end
bh = [];
bv = [];

if sum(sum(G))>100
    [~,threshold] = edge(G,'canny',[],1.0);
    k = 0.50;
    E = edge(G,'canny',k*threshold,1.0);

    [GradientX,GradientY] = gradient(double(G));
    GradientYY = gradient(GradientY);
    Gr = sqrt((GradientX.*GradientX)+(GradientY.*GradientY));
            
    index = GradientX == 0;
    GradientX(index) = 1e-5;
            
    YX = GradientY./GradientX;
    if angle == 180, A = ((atan(YX)+(pi/2))*180)/pi; end
    if angle == 360, A = ((atan2(GradientY,GradientX)+pi)*180)/pi; end
                                
    [bh bv] = anna_binMatrix(A,E,Gr,angle,bin);
else
    bh = zeros(size(I,1),size(I,2));
    bv = zeros(size(I,1),size(I,2));
end

[patchRowCount patchColCount] = size(patches);
for i = 1:patchRowCount
    for j = 1:patchColCount
        % roi - Region Of Interest (ytop,ybottom,xleft,xright)
        roi = patches{i,j};
        bh_roi = bh(roi(1,1):roi(2,1),roi(3,1):roi(4,1));
        bv_roi = bv(roi(1,1):roi(2,1),roi(3,1):roi(4,1));
        p{i,j} = anna_phogDescriptor(bh_roi,bv_roi,L,bin);
    end
end