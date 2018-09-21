function nstr = num2strn(n,d)
% takes an integer, then pads it so that it contains d digits
zerostr = '0';
nstr = num2str(round(n));
if (length(nstr) > d)
    error('integer has more characters than the total amount to be padded');
else
    padnum = (d - length(nstr));
    nstr = [repmat(zerostr,1,padnum) nstr];
end
