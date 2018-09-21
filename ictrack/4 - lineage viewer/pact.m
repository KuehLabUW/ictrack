function [gfpslope, complete, tlength, meantime] = pact(linname,color);

load(linname)
gfpslope = [];
rfpslope = [];
start = [];
finish = [];
complete = [];
startfig = 101;

L = length(lin);

minlen = 30;   % minimum length for a trajectory

for i = 1:L   % loop through all the lineages              
    crand = rand(1,3);   % random color;    
    trs = lin(i).trs;    
    
    for j = 1:length(trs)   % loop through all the cell lineages within a track
        
        if (trs(j).ingate)
            linestyle = '-';
        else            
            linestyle = '--';
        end
        
        t = thr(trs(j).ts);   % time in hours        
        if isempty(t)
            continue
        end        
        gfp = [trs(j).data.gfp];
        rfp = [trs(j).data.rfp];
        
        if (length(t) < minlen)
            continue;  % skip trajectory if too short
        end
        
        [fgfp, gof, output] = fit(t', gfp', 'poly1','robust','On');  % perform the robust fitting
        %frfp = fit(t', rfp', 'poly1','robust','On');  % perform robust fitting             
        % can we find a way to exclude the ones where the iteration limit
        % was reached?
        
        if (output.exitflag == 0)
            continue
        end
                
        gfp_val = fgfp.p1;
        conf = confint(fgfp);
        gfp_min     = conf(1,1);
        gfp_max = conf(2,1);
        
        figure(1); hold on;  % plot the gfp slope        
        axis([0 100 -500 8000]);
        plot([t(1) t(end)], [gfp_val gfp_val], 'Color',color,'LineStyle',linestyle,'LineWidth',0.5);
        
        % plot the         
        %text(t(1),fgfp.p1,num2str(trs(j).tr));
        
        if (j > 1)   % descendant, plot starting point
            plot(t(1), fgfp.p1, 'Marker', '<','Color', color,'MarkerFaceColor',color,'MarkerSize',5);
        end       
        if (~isempty(trs(j).children))
            plot(t(end), fgfp.p1, 'Marker', '>','Color', color,'MarkerFaceColor',color,'MarkerSize',5);
        end
        
        % is the trajectory complete?
        
        if ((j>1) && (~isempty(trs(j).children)))
            complete = [complete 1];
        else 
            complete = [complete 0];
        end
                        
%         figure(2);  % plot the rfp slope        
%         plot([t(1) t(end)], [frfp.p1 frfp.p1], 'LineStyle', '-','Color', crand);
%         if (j > 1)   % descendant, plot starting point
%             plot(t(1), frfp.p1, 'Marker', 'o','Color',crand);
%         elseif (~isempty(trs(j).children))
%             plot(t(end), frfp.p1, 'Marker', 'o','Color',crand);
%         end
%         hold on;
%         
%         figure(startfig); hold on; % plot the gfp slope                
%         plot([t(1) t(end)], [gfp_val gfp_val], 'LineStyle', '-','Color', crand);
%         plot([t(1) t(end)], [gfp_min gfp_min], 'LineStyle', ':','Color', crand);
%         plot([t(1) t(end)], [gfp_max gfp_max], 'LineStyle', ':','Color', crand);
%         
        

        %if (j > 1)   % descendant, plot starting point
        %    plot(t(1), fgfp.p1, 'Marker', 'o','Color',color);
        %end
        %if (~isempty(trs(j).children))
        %    plot(t(end), fgfp.p1, 'Marker', 'x','Color', color);
        %end
        axis([0 100 -500 8000]);        
        
%        if (fgfp.p1 > 4000)
%            keyboard
%        end
        
        gfpslope = [gfpslope fgfp.p1];
        start = [start t(1)];
        finish = [finish t(end)];
        
    end
    startfig = startfig+1;
    
    
end

figure; title('GFP promoter activities');
%figure; title('RFP promoter activities');

%figure(200); hist(gfpslope, 30);
meantime = (start + finish)./2;
tlength = (finish - start);   % length of trajectory
%figure(201); plot(meantime, gfpslope, 'kx');

%slopes = [];
%identity = [];
% % here we are classifying the promoter activities according to the mean
% % time of the trajectory.  Don't need to do it here
% for i = 1:length(meantime)  
%     if (meantime(i) < 10)
%         slopes = [slopes gfpslope(i)];
%         identity = [identity 0];
%     elseif ((meantime(i) > 20)&&(meantime(i) < 30))
%         slopes = [slopes gfpslope(i)];
%         identity = [identity 1];
%     elseif (meantime(i) > 60)
%         slopes = [slopes gfpslope(i)];
%         identity = [identity 2];
%     end
% end
% 
% early ; mid ; late
% 
% figure(3);
%     
    

        
            
        
        
        
        
        
            
       
        
        
        