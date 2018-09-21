% lineage information for PU-T 10 IL-3 10/26
clear all;

lingraph(1).nodes = 1:11;
lingraph(2).nodes = 12:22;
lingraph(3).nodes = 23:33;
lingraph(4).nodes = 34:48;
lingraph(5).nodes = 49:59;
lingraph(6).nodes = 60;
lingraph(7).nodes = 61:75;

% script to update the tracks.mat file for 
infile = 'D:\Rothenberg Lab\local data\101910 - putd123 ANALYZED\tracks4';
load(infile);
outfile = [infile ' addlin 3'];
    
for i = 1:length(lingraph)    
    nodes = lingraph(i).nodes;
    for j = 1:length(nodes)
        tracks(nodes(j)).lin = i;
    end
end
save([outfile '.mat'],'tracks');
