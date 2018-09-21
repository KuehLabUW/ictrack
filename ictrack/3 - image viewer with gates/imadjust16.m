function imout = imadjust16(im,c);
    cmin = c(1);
    cmax = c(2);
    imout = (im - cmin)*(65535/(cmax-cmin));