function plotlinflat(thrs,lintr,field,color,scale)   % plots a lineage tree out based on information from a lineage structure
 
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
                plot(thrs(ts), (out./scale),'Color',color,'LineWidth',1);        
            end
        end
        
        if isfield(edata,'death');
            %fprintf('field exists');            
            if ~isempty(edata.death);
                plot(thrs(ts(end)), (out(end)./scale), 'rx','MarkerSize', 12);
            end
        end                
    end
hold off;
        
        
            