% for completion of figure 2d


%[gfp_all, complete_all, tlength_all, meantime_all] = pact('schgate63_all',[0.8 0.8 0.8]);
[gfp_pro, complete_pro, tlength_pro, meantime_pro] = pact('schgate63_pro',[0.3 0.3 0.3]);
[gfp_b, complete_b, tlength_b, meantime_b] = pact('schgate63_b',[1 0 0]);
[gfp_mac, complete_mac, tlength_mac, meantime_mac] = pact('schgate63_mac',[0 0 1]);

figure(1); axis([0 100 -500 5000])


% histogram of promoter activities, fig 2d
figure(5);
bins = 40;
x = [-500:(5500/bins):5000];
n = hist(gfp_pro, x);
bar(x,n,1)
set(gca,'XLim',[-500 5000]);




% close all;
% figure(10); hist(gfp_pro(find(complete_pro)), 20); %set(gca,'Xlim',[3 9])
% 
% figure(11); hist(gfp_mac(find(complete_mac)), 20 ,'r'); %set(gca,'Xlim',[3 9])
% 
% 
% mean(gfp_pro(tlength(complete_pro)))
% 
% mean(gfp_mac(tlength(complete_mac)))
