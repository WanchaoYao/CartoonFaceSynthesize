function [clustered,C] = pixel_cluster(im,K)

if size(im,3) == 3
    pixs = double(reshape(im, size(im,1)*size(im,2), 3));
    [IDX,C] = kmeans(pixs,K,'distance','cityblock','emptyaction','singleton',...
        'start','uniform','Replicates',5);

    pixsCluster = reshape(IDX,size(im,1),size(im,2));

    clustered = zeros(size(im,1), size(im,2), 4);
    for i = 1:size(im,1)
        for j = 1:size(im,2)
            clustered(i,j,1:3) = C(pixsCluster(i,j),:);
            clustered(i,j,4) = pixsCluster(i,j);
        end
    end

    clustered = uint8(clustered);

    %figure, imshow(im);
    figure, imshow(clustered(:,:,1:3));
end

if size(im,3) == 1
    pixs = double(reshape(im, size(im,1)*size(im,2), 1));
    [IDX,C] = kmeans(pixs,K,'distance','cityblock','emptyaction','singleton',...
        'start','uniform','Replicates',10);

    pixsCluster = reshape(IDX,size(im,1),size(im,2));

    clustered = zeros(size(im,1), size(im,2), 2);
    for i = 1:size(im,1)
        for j = 1:size(im,2)
            clustered(i,j,1) = C(pixsCluster(i,j),:);
            clustered(i,j,2) = pixsCluster(i,j);
        end
    end

    clustered = uint8(clustered);

    %figure, imshow(im);
    figure, imshow(clustered(:,:,1));
end