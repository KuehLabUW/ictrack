function imc = myimadjust(im, xlo, xhi)

imc = (im - xlo)./(xhi - xlo);
imc(find(imc < 0)) = 0;
imc(find(imc > 1)) = 1;