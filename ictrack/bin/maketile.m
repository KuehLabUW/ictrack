function imtile = maketile(imall, Xo, Yo, method);

Y = size(imall,1);
X = size(imall,2);
M = size(imall,3);
N = size(imall,4);

% method - gives the projection method, maximum intensity?  or other
% methods?


inx = NaN.*zeros(M*(Y-Yo)+Yo, N*(X-Xo)+Xo,M,N);
imtile = zeros(M*(Y-Yo)+Yo, N*(X-Xo)+Xo);

s = whos('imall');
s = s.class;

switch s    
    case 'logical'
        imtile = logical(imtile);
        if strmatch(method,'avg')
            method = 'maxp';   % cannot use averaging for binary images
        end
    case 'uint8'
        imtile = uint8(imtile);
    case 'uint16'
        imtile = uint16(imtile);
end

if strmatch(method,'avg')
    num = imtile;
end

for m= 1:M
    for n = 1:N
        xmin = (n-1)*(X-Xo)+1;
        ymin = (m-1)*(Y-Yo)+1;                
        xmax = xmin+X-1;
        ymax = ymin+Y-1;
        
        switch method
            case 'maxp'
                imtile(ymin:ymax, xmin:xmax) = max(imall(:,:,m,n), imtile(ymin:ymax,xmin:xmax)); 
            case 'minp'
                imtile(ymin:ymax, xmin:xmax) = min(imall(:,:,m,n), imtile(ymin:ymax,xmin:xmax)); 
            case 'avg'
                imtile(ymin:ymax, xmin:xmax) = imall(:,:,m,n) + imtile(ymin:ymax,xmin:xmax);
                num(ymin:ymax, xmin:xmax) = num(ymin:ymax, xmin:xmax) +1;
        end
    end
end

if strmatch(method,'avg')
    imtile = imtile./num
end
        

