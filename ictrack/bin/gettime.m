function t = gettime(acq, segchannel)

tseg = acq.C(segchannel).tlist;
tcommon = find(tseg);   % the list of all planes where segmentation channel is present
t = (acq.tr - acq.tr(1)).*24;   % time in hours
t = t(tcommon);     % this is the time in hours