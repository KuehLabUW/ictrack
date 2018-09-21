function plotlin(thrs,lintr,field,color,scale)   % plots a lineage tree out based on information from a lineage structure

        % thrs -- the vector of time in hours       
    for i = 1:length(lintr)                
        % plots the parent trajectory
        tr = lintr(i).tr;   % display the track number?
        gen = lintr(i).gen;
        delta = lintr(i).delta;
        ts = lintr(i).ts;
        children = lintr(i).children;        
        deltas = delta.*ones(size(ts));
        edata = lintr(i).edata;
        
        if (length(ts)<1)
            continue
        end
        hold on;
        if (nargin > 2)   % plot additional fluorescence intensity data additional arguments are there            
            % check to see whether             
            if (isfield(lintr(i).data, field));                        
                out = [lintr(i).data.(field)];                
                area(thrs(ts), delta+(out./scale),delta,'EdgeColor','k','FaceColor',color);        
            end
        end        
        plot(thrs(ts), deltas,'k','LineWidth',1.5);                                      
        text(thrs(ts(1)), delta-0.05, num2str(tr));    % comment or
        %uncomment if necessary.  this is the track number
            % don't plot text, for
        % the mean time
        
        if isfield(edata,'death');
            %fprintf('field exists');            
            if ~isempty(edata.death);
%                 plot(thrs(ts(end)), delta, 'rx','MarkerSize', 12);
                  % death parameter comment or uncomment if necessary
            end
        end
        
        % plots connectors to children
        if isempty(children)  % plot vertical lines for cell division events            
            continue
        end                        
        for i = 1:length(children)            
            c = children(i);            
            plot([thrs(ts(end)) thrs(lintr(c).ts(1))], [delta, lintr(c).delta;],'k', 'LineWidth', 1.5);
        end                
    end
hold off;
        
        
            