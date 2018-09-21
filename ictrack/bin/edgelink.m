function imlink = edgelink(imb, nb)
%UNTITLED IMLINK = EDGELINK(IMB, NB) takes a binary image imb and links
%pixels a distance nb or less away that are not already connected

if (nb == 0)
    imlink = imb;   % do no linking of the neighborhood size is zero
end

I = size(imb,1);
J = size(imb,2);
imlink = imb;   % output image

p = 0:(0.3/nb):1;   % line parameterization

for i = 1:I
    
    for j = 1:J
        
        if (imb(i,j) == 0) % black pixel no need to link
            continue;
        end
        
        imin = max(1,i-nb);
        imax = min(I,i+nb);
        jmin = max(1,j-nb);
        jmax = min(J,j+nb);
        
        imbsub = imb(imin:imax, jmin:jmax);
        imlsub = bwlabel(imbsub);
        
        curr_ind  = imlsub(i-imin+1,j-jmin+1);    % index of current object
        not_connected = ~ismember(imlsub, [0 curr_ind]);   % find objects that are not connected
        
        [inc, jnc] = find(not_connected);  
        inc = inc-1+imin;   % convert back to regular indices 
        jnc = jnc-1+jmin; 
        
        for k = 1:length(inc)  % fill in the 
            ilist = round(i + (inc(k)-i).*p);
            jlist = round(j + (jnc(k)-j).*p);
            indlist = sub2ind([I J], ilist, jlist);
            imlink(indlist) = 1;
        end
    end
end
           
end

