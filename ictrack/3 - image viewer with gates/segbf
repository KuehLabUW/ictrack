function imout = segbf(im, xnow, ynow)

blur = 3;
crange = 5;
canlo = 0.2;
canhi = 0.5;
bw = 3; 

X = size(im,2);
Y=  size(im,1);
regsize = 30;

xmin = max(1, xnow-regsize);
ymin = max(1, ynow-regsize);
xmax = min(X, xnow+regsize);
ymax = min(Y, ynow+regsize);
ynew = ynow - ymin;
xnew = xnow - xmin;
imsub = im(ymin:ymax, xmin:xmax);  % do image calculations based on subimage

log = fspecial('log',10, 3);
imsub2 = imfilter(imsub,log,'Replicate');   % laplacian of gaussian
thresh = max(prctile(imsub2(:),90),1);   % threshold for edge, arbitrarily set
imb = (imsub2 > thresh);
imb2 = imopen(imb, [1 1; 1 1]);  % remove small isolated bits of crap
imb3 = edgelink(imb2, 5);
imb4 = imfill(imb3,'holes');
% now remove other objects
iml = bwlabel(imb4);
imb5 = ismember(iml, iml(ynew,xnew));
imb6 = imopen(imb5,strel('disk',8));

figure(10);
subplot(4,2,1); imshow(imsub,[]);
subplot(4,2,2); imshow(imsub2,[]);
subplot(4,2,3); imshow(imb);
subplot(4,2,4); imshow(imb2);
subplot(4,2,5); imshow(imb3);
subplot(4,2,6); imshow(imb4);
subplot(4,2,7); imshow(imb5);
subplot(4,2,8); imshow(imb6);
set(gcf,'HandleVisibility', 'off');
imout = imb6;
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
return