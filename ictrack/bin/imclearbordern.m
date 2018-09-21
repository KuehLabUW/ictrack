function imout = imclearbordern(im, n);
if (n>1)    
    im(1:n,:) = 0;
    im(:,1:n) = 0;
    im((end-n+1):end, :) = 0;
    im(:, (end-n+1):end) = 0;
end

imout = im;