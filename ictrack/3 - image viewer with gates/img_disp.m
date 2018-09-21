close all;
clear all;



load('acq.mat');    % load the acquisition parameters
load('img_0001.mat');  % timelapse data
% load the acquisition variables

C = acq.C;
T = acq.T;
M = acq.M;
N = acq.N;
Z = acq.Z;
X = 348;
Y = 260;
Xo = 65;   % overlap in X
Yo = 65;   % overlap in Y

%% initialize viewing window
figure(1);
axis equal;
axis ij;
axis([0 N*(X-Xo)+Xo 0 M*(Y-Yo)+Yo]);
colormap gray
hold on;

%% set the contrast limits
C(1).minc = 150;
C(1).maxc = 6000;
C(2).minc = 500;
C(2).maxc = 3000;
C(3).minc = 10000;
C(3).maxc = 60000;

for c = 1:length(C)
    C(c).cR = C(c).cR/255
    C(c).cG = C(c).cG/255
    C(c).cB = C(c).cB/255
end

% for tomorrow - incorporate this into the gui!
% then maybe add a feature for multiple channels...
% toggling, etc.......
% flatten background also!  important....
% update the pview structure.

tic
for m = 1:M
    for n = 1:N
        imrgb = uint16(zeros(Y,X,3));
        for c = 1:length(C);            
            im = imadjust16(images(c).im(:,:,m,n,1), [C(c).minc C(c).maxc]);
            imrgb = imrgb + cat(3, im*C(c).cR, im*C(c).cG, im*C(c).cB);
        end 
        tileh(m,n) = imagesc((n-1)*(X-Xo), (m-1)*(Y-Yo), imrgb);
    end
end
toc;     
 

