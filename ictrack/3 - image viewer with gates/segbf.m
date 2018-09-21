function imout = segbf(im, xnow, ynow)

blur = 3;
crange = 5;
canlo = 0.2;
canhi = 0.5;
bw = 3; 

X = size(im,2);
Y=  size(im,1);
regsize = 50;

xmin = max(1, xnow-regsize);
ymin = max(1, ynow-regsize);
xmax = min(X, xnow+regsize);
ymax = min(Y, ynow+regsize);
ynew = ynow - ymin;
xnew = xnow - xmin;
imsub = im(ymin:ymax, xmin:xmax);  % do image calculations based on subimage

imout = zeros(size(im));    % create an output matrix

log = fspecial('log',10, 2);
imsub2 = imfilter(imsub,log,'Replicate');   % laplacian of gaussian
thresh = max(prctile(imsub2(:),80),1);   % threshold for edge, arbitrarily set
imb = (imsub2 > thresh);
imb2 = imopen(imb, [1 1; 1 1]);  % remove small isolated bits of crap
imb3 = imclose(imb2,strel('disk',1));
imb4 = bwmorph(imb3,'thin','Inf');
imb5 = imclearbordern(imb4,5);   % remove objects close to the boundary;   
imdisk = zeros(size(imb5)); imdisk(ynew, xnew) = 1; imdisk = imdilate(imdisk,strel('disk',4)); % make
imb5 = (imb5 & (~imdisk));     % remove objects in the neighborhood for the center 

iml = bwlabel(imb5);   % find objects that are super extended
lens = regionprops(iml,'Area');
long = find([lens.Area] > 20);   % find long line segments
imb5 = ismember(iml,long);    

% only link edges if region is not closed --
linklen = 0;
maxlink = 15;
imb6 = zeros(size(imb5));    % the default image is empty if image cannot be closed

for i = 1:maxlink   % try increasing link lengths
    imbtest = edgelink(imb5, i);
    imbtestfill = imfill(imbtest,'holes');   % see if the point is now in a closed region
    %figure; subplot(2,1,1); imshow(imbtest)
    %subplot(2,1,2); imshow(imbtestfill); title(['edgelink length = ' num2str(i)]);
    if imbtestfill(ynew,xnew);  % region is closed
        imb6 = imfill(imbtest,[ynew, xnew]);   % fill only the region specified
        break;  % done with the loop!
    end
end 
imb7 = imopen(imb6,strel('disk',6));    % now remove debris from the that is not part of the cell

figure(10);
subplot(3,3,1); imshow(imsub,[]);
subplot(3,3,2); imshow(imsub2,[]);
subplot(3,3,3); imshow(imb);
subplot(3,3,4); imshow(imb2);
subplot(3,3,5); imshow(imb3);
subplot(3,3,6); imshow(imb4);
subplot(3,3,7); hold off; imshow(imb5); hold on; plot(xnew,ynew,'rx');
subplot(3,3,8); imshow(imb6); title(['linklen=' num2str(i)]);
subplot(3,3,9); imshow(imb7);
%subplot(4,3,10); imshow(imb8);
figure(1);
imout(ymin:ymax, xmin:xmax) = imb7;

% need to reconvert this into 


%figure(4); imshow(imb2);
%a = waitforbuttonpress
        % calculate image histogram limits based on subregion centered
        
        %[xp, dr, dl] = findxpw(imsub);
        %imedge = edge(im,'canny',[canlo canhi]);
        %% perform edge detection using canny filter, then correct
        %% artifacts
        %
        %imedge = imclose(imedge,[1 1]); imedge = imclose(imedge,[1;1]);  % bridge 4-connected regions
        %imedge = imfill(imedge,'holes');  % fill all closed regions
        %imedge = imopen(imedge, [1 1; 1 1]);  % disconnect objects that
        %are separated by one object
