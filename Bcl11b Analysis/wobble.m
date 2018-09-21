function [w] = wobble(i9,im)
props = regionprops(i9,im,'Centroid'); % Find centroids of each object
bw2 = bwperim(i9); % Take only perimeter pixels of each object
props2 = regionprops(bw2,'PixelList'); % List pixel locations on perimeter for each object
centroids = [props.Centroid]'; % Vector of centroid values
mcentroids = zeros(); % Blank matrix for centroid values
k = 1; % Row index for centroid matrix
for j = 1:2:length(centroids)
    mcentroids(k,1) = centroids(j); % X value of centroid in first column of matrix
    mcentroids(k,2) = centroids(j+1); % Y value of centroid in second column of matrix
    k = k+1; % Increment row index
end
num_obj = length(centroids)/2; % Number of objects
pixel_perimeter = zeros(); % Blank matrix for pixel perimeter locations
num_pix = 0; % Initial pixel perimeter number
w = zeros();
for i = 1:num_obj
    d = zeros(); % Blank vector for distances from perimeter pixels to centroid
    pixel_perimeter = props2(i).PixelList; % Pixel perimeter locations for object i
    num_pix = length(pixel_perimeter(:,1)); % Number of pixels in perimeter of object i
    for n = 1:num_pix
        d(n) = pdist([mcentroids(i,:); pixel_perimeter(n,:)]); % Vector of distances from perimeter pixel to centroid
    end
    [poly3fit, gof] = fit([1:num_pix]',d','poly3'); % Fit 3rd degree polynomial to distance values
    w(i) = gof.rmse/num_pix; % Wobble is root mean square error of curve/number of pixels in perimeter
    %figure(1000+i)
    %plot(poly3fit,1:num_pix,d)
end
end
