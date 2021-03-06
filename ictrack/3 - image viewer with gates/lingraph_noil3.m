% lineage information for PU-T 10 IL-3 10/26
% lingraph(1).nodes = 1:11;
% lingraph(2).nodes = 12:22;
% lingraph(3).nodes = 23:33;
% lingraph(4).nodes = 34:48;
% lingraph(5).nodes = 49:59;
% lingraph(6).nodes = 60;
% lingraph(7).nodes = 61:75;

% script to update the tracks.mat file for 
infile = 'D:\Rothenberg Lab\local data\101910 - putd123 ANALYZED\no il3 tracks';
load(infile);
outfile = [infile 'addlin'];
    
% lineage information for PU-T 0.08 IL3 10/26
lingraph(1).nodes = 1:3;
lingraph(2).nodes = 4:6;
lingraph(3).nodes = 7:17;
lingraph(4).nodes = 18:19;
lingraph(5).nodes = 20:22;
lingraph(6).nodes = 23:29;
lingraph(7).nodes = 30:40;
lingraph(8).nodes = 41:53;
lingraph(9).nodes = 54;

for i = 1:length(lingraph)    
    nodes = lingraph(i).nodes;
    for j = 1:length(nodes)
        tracks(nodes(j)).lin = i;
    end
end
save([outfile '.mat'],'tracks');
