close all; imtool close all; clear all;

wtotal = 376; % Total number of well positions
well = round(rand()*wtotal)% Generate random well position from total number of wells 
ch = 3; % Fluorescence channel for analysis, 3 = CFP
htype = 1; % Type of histogram generated: 1 = area of objects, 2 = eccentricity of objects, 3 = wobble of objects, 4 = ratio of objects

for j = 1:40:177 % For images 1,41,81,121, and 161
    cd(['/data/kueh/060117 - Bcl11b-YFP-RFP-CFP ANALYZED/pos ' num2str(well)]) % Enter folder for well position
    
    if j == 1
        im = load(['imgf_000' num2str(j) '.mat']); % Load image structure at j index
        img = double(im.images(ch).im); % Convert image to double
        cd('/data/kueh/060117 - Bcl11b-YFP-RFP-CFP code/') % Enter folder with cellseg function
        [seg] = cellseg060117(img,j); % Segment image
        [B,imb] = bwboundaries(seg,8,'noholes'); % Label objects
        imtool(imb) % Display segmented image
        props = regionprops(imb,'Area','Perimeter','MajorAxisLength','MinorAxisLength'); % Find area and major and minor axis length of each object
        if htype == 1
            figure(j+1)
            hist([props.Area],100) % Plot histogram of object areas
        elseif htype == 2
            figure(j+2)
            hist([props.MajorAxisLength]./[props.MinorAxisLength],100) % Plot histogram of object eccentricities
        elseif htype == 3
            figure(j+3)
            hist(w,100) % Plot histogram of object wobble
        elseif htype == 4
            figure(j+4)
            hist([props.Perimeter].^(2)./[props.Area],100) % Plot histogram of object perimeter squared/area ratio
        end
        
    elseif j == 41 || j == 81
        im = load(['imgf_00' num2str(j) '.mat']); % Load image structure at j index
        img = double(im.images(ch).im); % Convert image to double
        cd('/data/kueh/060117 - Bcl11b-YFP-RFP-CFP code/') % Enter folder with cellseg function
        [seg] = cellseg060117(img,j); % Segment image
        [B,imb] = bwboundaries(seg,8,'noholes'); % Label objects
        imtool(imb) % Display segmented image
        props = regionprops(imb,'Area','Perimeter','MajorAxisLength','MinorAxisLength'); % Find area and major and minor axis length of each object
        if htype == 1
            figure(j+1)
            hist([props.Area],100) % Plot histogram of object areas
        elseif htype == 2
            figure(j+2)
            hist([props.MajorAxisLength]./[props.MinorAxisLength],100) % Plot histogram of object eccentricities
        elseif htype == 3
            figure(j+3)
            hist(w,100) % Plot histogram of object wobble
        elseif htype == 4
            figure(j+4)
            hist([props.Perimeter].^(2)./[props.Area],100) % Plot histogram of object perimeter squared/area ratio
        end
        
    elseif j == 121 || j == 161
        im = load(['imgf_0' num2str(j) '.mat']); % Load image structure at j index
        img = double(im.images(ch).im);
        cd('/data/kueh/060117 - Bcl11b-YFP-RFP-CFP code/') % Enter folder with cellseg function
        [seg] = cellseg060117(img,j); % Segment image
        [B,imb] = bwboundaries(seg,8,'noholes'); % Label objects
        imtool(imb) % Display segmented image
        props = regionprops(imb,'Area','Perimeter','MajorAxisLength','MinorAxisLength'); % Find area and major and minor axis length of each object
        if htype == 1
            figure(j+1)
            hist([props.Area],100) % Plot histogram of object areas
        elseif htype == 2
            figure(j+2)
            hist([props.MajorAxisLength]./[props.MinorAxisLength],100) % Plot histogram of object eccentricities
        elseif htype == 3
            figure(j+3)
            hist(w,100) % Plot histogram of object wobble
        elseif htype == 4
            figure(j+4)
            hist([props.Perimeter].^(2)./[props.Area],100) % Plot histogram of object perimeter squared/area ratio
        end
    end
end