classdef Schnitz3 % < handle
    
    properties
        timepoints;   % contains the objects for each timepoint
        tracks;       % contains information for individual tracks        
    end
    
    properties (Dependent = true)
        T;            % the total number of timepoints
        Tr;           % the largest track number in the whole object        
        tracklist;    % the list of tracks that exist   
        tr_approved;  % list of all approved tracks
    end
    
    methods        
        %% constructor function
        function SCH = Schnitz3(timepoints, varargin)
            % load objects for each time point                                    
            SCH.timepoints = timepoints;                        
            if (nargin == 1)   % need to generate tracks structure
                Tr = 0;
                for i = 1:length(timepoints)   % loop through all timepoints
                    if isempty(timepoints(i).obj)
                        continue
                    end
                    trnos = [timepoints(i).obj.trno];
                    Tr = max(Tr,max(trnos));   % get the largest track number                     
                end
                
                % if the data subfield does not exist for the object
                if ~isfield(timepoints(1).obj,'data');
                    fprintf('need to add data field to timepoints structure.\n');
                    for i = 1:length(timepoints)
                        obj = timepoints(i).obj;
                        if (~isempty(obj))
                            SCH.timepoints(i).obj(1).data = [];
                        end
                    end
                end
                
                % generate the tracks structure, with color information
                tracks.approved = 0;
                tracks.parent = [];
                tracks.children = [];
                tracks.data = [];
                tracks = repmat(tracks,1,Tr);
                for i = 1:Tr
                    tracks(i).rgb = rand(1,3)*0.9 + 0.1;;  % generate random colors for the tracks
                end
            else        % tracks structure already present
                tracks = varargin{1};
            end
            % store data in the object
            
            SCH.tracks = tracks;            
        end
        
        %% dependent variable get functions
        function val = get.Tr(SCH)     % get function for the total number of tracks
            val = length(SCH.tracks);
        end       
        
        function val = get.T(SCH)     % get function for the total number of timepoints
            val = length(SCH.timepoints);
        end       
         
        function val = get.tracklist(SCH)   % get list of tracks that exist in the structure 
            timepoints = [SCH.timepoints];
            obj = [timepoints.obj];
            tracklist = [obj.trno];   % list of the all the tracks, with overlapping indices
            val = unique(tracklist);
            val(find(val==0)) = [];
        end                
        
        %% get functions for timepoints
        function val = tr(SCH, t, varargin);   % get the track number for object with given time point T and index IND            \            
            if (nargin == 3)                   % return all tracks for a given time point if no index specified
                ind = varargin{1};
                val = SCH.timepoints(t).obj(ind).trno;
            else
                obj = SCH.timepoints(t).obj;
                val = [obj.trno];
            end
        end                
       
        function val = obj(SCH, t, varargin);   % get objects with given time point T (and index IND if specified)
            if (nargin == 3)   % index specified             
                ind = varargin{1};                
                val = SCH.timepoints(t).obj(ind);
            else    % no index specified
                val = SCH.timepoints(t).obj;
            end
        end                
        
        function val = objdata(SCH, t, ind);        % get data correspondind to a given object                              
            val = SCH.timepoints(t).obj(ind).data;
        end  
        
        %% get functions for tracks
        function val = rgb(SCH, tr);     % get the rgb information for a given track number
            val = SCH.tracks(tr).rgb;
        end
        
        function val = children(SCH, tr);
            val = SCH.tracks(tr).children;
        end
        
        function val = approved(SCH, tr);            
            if (nargin == 2);                   
                if (tr)
                    if (tr>0)
                        val = SCH.tracks(tr).approved;   % query only one track
                    else
                        val = 0;
                    end
                    
                else
                    val = 0;
                end
            else
                val = [SCH.tracks.approved];     % query all tracks
            end
        end
        
        function [ts, os] = trackobj(SCH, tr);    % find objects associated with single track                                    
            iall = [];
            tall = [];            
            timepoints = SCH.timepoints;            
            for i = 1:length(timepoints);                
                L = length(timepoints(i).obj);                
                if (L>0)   % objects exist at this timepoint
                    iall = [iall 1:L];   % create long vector with all object index numbers
                    tall = [tall i.*ones(1,L)];   % create long vector with all time points
                end
            end            
            oall = [timepoints.obj];
            trall = [oall.trno];            
            inds = find(trall == tr);            
            ts = tall(inds);
            os = iall(inds);                  
        end               
           
        function out = istrack(SCH, tr)        % returns 1 if trno is a track, otherwise returns zero            
            if (isempty(intersect(SCH.tracklist, tr)))
                out = 0;
            else
                out = 1;
            end
        end
        
        function val = trdata(SCH, tr);           % get the event data for the track
            val = SCH.tracks(tr).data;
        end
        
        
        %% get functions for lineage information          
        function root = getlinroot(SCH, tr)   % find the root of a given track        
            root = tr;
            parent = SCH.tracks(root).parent;
            while (~isempty(parent))
                root = parent;
                parent = SCH.tracks(root).parent;
            end
        end            
            
        function lintr = getlintr(SCH, tr)            
            % do pre-order DFS algorithm            
            node = SCH.getlinroot(tr);   % starting node                                    
            lintr = [];            
            lintr = preorder(lintr, node, 1, 0, SCH);                                    
            
            function lintr = preorder(lintr, node, gen, delta, SCH)                                               
                if (gen == 1)
                    i = 1;
                else
                    i = length(lintr)+1;
                end                
                % obtain track data for all objects in this track
                [ts, inds] = SCH.trackobj(node);
                data = [];                
                for t = 1:length(ts)
                    data = [data SCH.objdata(ts(t),inds(t))];
                end                
                % obtain the event data for this track
                edata = SCH.trdata(node);                
                % children
                children = SCH.children(node);      
                                         
                % store lineage information
                lintr(i).tr = node;
                lintr(i).gen = gen;
                lintr(i).delta = delta;                
                lintr(i).ts = ts;
                lintr(i).data = data;
                lintr(i).edata = edata;
                lintr(i).children = children; % update the lineage structure                                                
                
                % error checking to prevent self-references
                node
                children
                
                
                if (~isempty(children))  % perform recursive search of binary tree
                    
                    if (node == children(1))                        
                        beep
                        fprintf('self-reference')
                        keyboard
                    end
                    
                    lintr = preorder(lintr, children(1), gen+1, delta-(0.5^gen), SCH);                    
                    if (length(children)>1)  % other daughter exists       
                        if (node == children(2))
                            beep
                            fprintf('self-reference')
                            keyboard
                        end                        
                        lintr = preorder(lintr, children(2), gen+1, delta+(0.5^gen), SCH);
                    end
                end
            end                        
            % re-number the children according to the lin tracks             
            trs = [lintr.tr];            
            for i = 1:length(lintr)                
                c = lintr(i).children;
                if isempty(c)
                    continue
                end
                for j = 1:length(c)
                    lintr(i).children(j) = find(trs==c(j));
                end
            end            
        end
                      
        %% set functions for timepoints
                
        function SCH = delobj(SCH, t, ind)   % delete object
            SCH.timepoints(t).obj(ind) = [];
        end
            
        function SCH = addobj(SCH, t, objnew)   % add new objects to the current timepoint
            SCH.timepoints(t).obj = [SCH.timepoints(t).obj objnew];  % append the newly generated objects                    
        end
        
        function SCH = addobjdata(SCH, t, ind, data);    % set data for a particular object
            SCH.timepoints(t).obj(ind).data = data;
        end
        
        %% set functions for tracks
        
        function SCH = approvetrack(SCH, tr, varargin)            
            if (nargin == 2)                
                SCH.tracks(tr).approved = 1;
            elseif (nargin == 3)
                flag = varargin{1};
                if isempty(intersect([1 0], flag));
                    error('Approval flag must be 1 or 0');
                end
                SCH.tracks(tr).approved = flag;
            end
        end
        
        function SCH = approveall(SCH)   % approve all lineages that have parents or children            
            for i = 1:length(SCH.tracks)                
                if (~isempty(SCH.tracks(i).children))|(~isempty(SCH.tracks(i).parent))                    
                    if (SCH.tracks(i).approved)
                        fprintf(['Track no. ' num2str(i) ' already approved.\n']);
                    else
                        SCH.tracks(i).approved = 1;
                        fprintf(['Approving track no ' num2str(i) '.\n']);
                    end
                end
            end
        end
        
        function SCH = addtotrack(SCH, tr, t, ind)   % add an unassigned object to a track            
                                   
            if (tr == 0)   % if trno is zero, then generate a new track
                SCH = SCH.newtrack();
                tr = SCH.Tr;
            end             
            tracks = SCH.tr(t);   % list of tracks at this given timepoint                       
            if (~tracks(ind))   % object not yet assigned to a track                
                ind_old = find(tracks==tr);   % index of object already assigned to the track at the timepoint, if exists                                                 
                if (~isempty(ind_old))
                   SCH = SCH.removefromtrack(tr, t, ind_old);
                end
                SCH.timepoints(t).obj(ind).trno = tr;                                            
            else
                error('Object already belongs to a track.')
            end
        end
        
        function SCH = removefromtrack(SCH,tr, t, ind)   % remove an object from a track
            if (SCH.timepoints(t).obj(ind).trno==tr)
                SCH.timepoints(t).obj(ind).trno = 0;
            else
                error('Object does not belong to specified track.');
            end
        end
        
        function SCH = newtrack(SCH);             % generate a new track
            trnew = SCH.Tr+1;           
            SCH.tracks(trnew).approved = 0;
            SCH.tracks(trnew).parent = [];
            SCH.tracks(trnew).children = [];                                    
            SCH.tracks(trnew).rgb = rand(1,3)*0.9 + 0.1;  % generate random colors for the tracks            
        end
        
        function SCH = deletetrack(SCH, tr);             % remove entire track             
            if (SCH.istrack(tr))
                [ts, inds] = SCH.trackobj(tr);
                for i = 1:length(ts)
                    SCH.timepoints(ts(i)).obj(inds(i)).trno = 0;
                end            
            else
                error('track does not exist.');
            end
        end
                            
        function SCH = splittrack(SCH, tr, t, varargin);                        
            % function SPLITTRACK(TRNO, T, VARARGIN) splits an existing
            % track into two.  Objects from all timepoints < t are assigned
            % to the first track, objects from all timepoints >= t are
            % assigned to the second track.
            % By default, the existing track number is assigned to objects
            % with timepoints < t; the passing of an additional argument
            % will cause the existing track number to be assigned to
            % objects with timepoints >= t.                                    
            if (~SCH.istrack(tr))
                error('Specified track does not exist');               
            end            
            [ts, inds] = SCH.trackobj(tr);  % specified             
            
            ind1 = find(ts < t);   % indices of all objects before T
            ind2 = find(ts >= t);  % indices of all objects on or after T
            
            % re-assign objects to new track            
            if ((nargin > 3)&(~isempty(ind1)))   % assign the earlier timepoints to the new track                                
                SCH = SCH.newtrack();   % create new track            
                for i = ind1
                    SCH.timepoints(ts(i)).obj(inds(i)).trno = SCH.Tr;
                end                
                parent = SCH.tracks(tr).parent;                                     
                if (~isempty(parent))                
                    SCH.tracks(SCH.Tr).parent = parent;
                    children = SCH.tracks(parent).children;                    
                    SCH.tracks(parent).children(find(children == tr)) = SCH.Tr;
                end                
                SCH.tracks(tr).parent = [];                  
            elseif ((nargin == 3)&(~isempty(ind2)))   % assign the later timepoints to the new track                                
                % re-assign any children to the new track, orphanize old track
                SCH = SCH.newtrack();
                for i = ind2                
                    SCH.timepoints(ts(i)).obj(inds(i)).trno = SCH.Tr;   % assign all the trajectories to the new track
                end                
                children = SCH.tracks(tr).children;    
                if (~isempty(children))
                    SCH.tracks(end).children = children;
                    for  j = 1:length(children)
                        SCH.tracks(children(j)).parent = SCH.Tr;
                    end
                end
                SCH.tracks(tr).children = [];                
            else
                %fprintf('nothing to reassign...');
            end
        end                
        
        function SCH = jointrack(SCH, tr1, t1, tr2, t2);            
            % join tracks tr1 and tr2, assigning the new track number tr1
            % to both the union of both tracks
            % the difference in times determine the temporal sequence
            if (t1 == t2)
                error('Cannot have two coincidental timepoints...\n');
            elseif (t1 < t2)  % tr1 precedes tr2
                SCH = SCH.splittrack(tr1, t1+1);   % disconnect tr1 from all subsequent data points
                SCH = SCH.splittrack(tr2, t2);     % disconnect tr2 from all previous time points                                
                children = SCH.tracks(SCH.Tr).children;
                if (~isempty(children))
                    SCH.tracks(tr1).children = children;    % reassign children of new track to Tr1
                    for  j = 1:length(children)
                        SCH.tracks(children(j)).parent = tr1;
                    end
                end
                SCH.tracks(SCH.Tr).children = [];
            else
                SCH = SCH.splittrack(tr1, t1, 'earlier');  
                SCH = SCH.splittrack(tr2, t2+1, 'earlier');
                
                parent = SCH.tracks(SCH.Tr).parent;                                                     
                if (~isempty(parent))
                    SCH.tracks(tr1).parent = parent;
                    children = SCH.tracks(parent).children;                    
                    SCH.tracks(parent).children(find(children == SCH.Tr)) = tr1;
                end                
                SCH.tracks(SCH.Tr).parent = [];                                  
            end
            Tr = SCH.Tr;
            [ts, inds] = SCH.trackobj(Tr);  % obtain indices of all objects to be assigned to tr1
            if (~isempty(ts))
                for i = 1:length(ts)
                    SCH.timepoints(ts(i)).obj(inds(i)).trno = tr1;   % reassign new objects to tr1
                end
            end            
            
        end                            
        
        function SCH = delchild(SCH, tr)            
            cs = SCH.tracks(tr).children;            
            if isempty(cs)
                return % nothing to delete
            end
            
            for i = 1:length(cs)
                SCH.tracks(cs(i)).parent = [];
            end
            SCH.tracks(tr).children = [];
        end
            
        function [SCH, added] = addchild(SCH, tr1, tr2)   % tr1 for parent, tr2 for child            
            children = SCH.tracks(tr1).children ;
            parent = SCH.tracks(tr2).parent;                        
            if (length(children) == 2)
                fprintf('Two children already exist for current track\n');
                added = 0;
                return
            end
            %elseif ~isempty(parent)
            %    fprintf('Selected child already parented\n');
            %    added = 0;
            %    return                
            %end            
            SCH.tracks(tr1).children = [children tr2];
            SCH.tracks(tr2).parent = tr1;
            added = 1;
        end              
        
        function SCH = addtrdata(SCH, tr, data);    % set data for a given track
            SCH.tracks(tr).data = data;
        end
        
        %% set functions for lineages        
        function SCH = checklins(SCH)            
            for i = 1:length(SCH.tracks)
                p = SCH.tracks(i).parent;             
                if ~isempty(p)
                    c = SCH.tracks(p).children;
                    if isempty(intersect(c,i));   % parent does not recognize this child
                        SCH.tracks(i).parent = [];     % remove the parent from the child
                        fprintf(['De-parented track number ' num2str(i) '\n']);
                    end
                end
            end
        end
           
    end   % methods
end   % classdef
