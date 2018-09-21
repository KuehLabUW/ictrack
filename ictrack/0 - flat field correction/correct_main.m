indir= '~/060315 - Bcl11b-YFP imaging/060815-beadscal';

base = 'beads3';   % the base name for the files

channel = 'w2YFP';   % the name of the channel
out = '060315-YFP.mat'
correct(indir,base,channel,out);

channel = 'w3RFP';
out = '060315-RFP.mat'
correct(indir,base,channel,out);