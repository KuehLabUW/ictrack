function imout = myunsharp(im, um, alpha);

        imout = im - alpha.*imfilter(im,um);