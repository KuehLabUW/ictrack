function im = readmmim(indir,t,cname,z,n,m);
% reads an individual micromanager file from a multi-dimensional
% acquisition

dirname = [indir '/Pos_' num2strn(n-1,3) '_' num2strn(m-1,3)];
filename = ['img_' num2strn(t-1,9) '_' cname '_' num2strn(z-1,3) '.tif'];
im = imread([dirname '/' filename]);


end

