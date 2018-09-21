classdef Fcs
    
    properties
        data        % contains the raw data points
        fieldinfo   % contains field names and display ranges
        gates       % contains information about specific gates
    end
    
    %properties (Dependent = true)
        %fnames  % cell array containing names of all the fields
        %gnames  % cell array containing names of all the gates
    %end
    
    methods
        %% constructor function
        function FCS = Fcs(timepoints, thr, varargin)
            % FCS(timepoints, thr)          uses timepoints info, timestamps
            % FCS(timepoints, thr, tracks)     uses above and track information
            % FCS(timepoints, thr, tracks, fieldinfo, gates)  gate info
            % FCS(timepoints, thr, [], fieldinfo, gates) gates, no tracks
            
            %% process optional arguments
            if (nargin == 2)
                withgate = 0; 
                withtrack = 0;
            elseif (nargin == 3)
                fprintf('withtracks\n');
                withtrack = 1;                
                tracks = varargin{1};
                withgate = 0;
            elseif (nargin == 5)
                withgate = 1;
                fieldinfo = varargin{2};
                gates = varargin{3};                
                if isempty(varargin{1})
                    withtrack = 0;
                else                    
                    withtrack = 1;
                    tracks = varargin{1};
                end
            else
                error('Unknown input argument structure...');
            end            
            
            %% populate data structure
            dtest = timepoints(1).obj(1).data;            
            isyfp = isfield(dtest,'YFP');
            isgfp = isfield(dtest,'gfp');
            isrfp = isfield(dtest,'RFP');
            iscfp = isfield(dtest,'cfp');
            iscy5 = isfield(dtest,'cy5');
            isA488 = isfield(dtest,'A488');
            isapc = isfield(dtest,'apc');

            
            
            data = [];            
            for t = 1:length(timepoints)    
                
                obj = timepoints(t).obj;  
                
                t
                
                if isempty(obj)   % skip frame if no data here
                    continue;
                end
                
                data1 = cat(2,obj.data);
                
                delobj = [];      % array for objects that have been deleted
                for i = 1:length(data1)                    
                    
                    data1(i).t = thr(t);   % this is the time in hours
                    data1(i).tint = t;      % this is the time in frames
                    data1(i).ind = i;
                    
                    if isyfp
                        data1(i).lyfp = log10(data1(i).YFP);
                    end                    
                    
                    if isapc
                        data1(i).lapc = log10(data1(i).apc);
                    end
                    
                    if isA488
                        data1(i).lA488 = log10(data1(i).A488);
                    end
                    
                    if isgfp
                        data1(i).lgfp = log10(data1(i).gfp);
                    end
                    
                    if isrfp
                        data1(i).lrfp = log10(data1(i).RFP);
                    end
                    
                    if iscfp
                        data1(i).lcfp = log10(data1(i).cfp);
                    end
                    
                    if iscy5
                        data1(i).lcy5 = log10(data1(i).cy5);
                    end
                    
                    %% enter track information if present
                    if (withtrack)                        
                        trno = obj(i).trno;
                        if (trno == -1)
                            delobj = [delobj i];
                            continue;
                        end                        
                        if (trno)
                            data1(i).trno = trno;
                            data1(i).approved = tracks(trno).approved;
                        else
                            data1(i).trno = 0;
                            data1(i).approved = 0;
                        end
                    end                    
                end
                data1(delobj) = [];   % delete the deleted objects                                                
                if ~isempty(data1)
                    data = [data data1];                
                end
            end
            
            FCS.data = data;
            
            %% initialize the field info structure            
            if (~withgate)    % initialize new gates and fields                                
                bins = 100;     % this is the default number of bins
                names = fieldnames(data);
                for i = 1:length(names)
                    fieldinfo(i).name = names{i};   % the name
                    data1 = [data.(names{i})];   % get the list of all the data points for this field 
                    dmin = min(data1);    % find the minimum
                    dmax = max(data1);    % find the maximum
                    fieldinfo(i).range = dmin:(dmax-dmin)./bins:dmax;     % the display range                    
                    if strmatch(names{i},'t')
                        fieldinfo(i).range = thr;
                    end
                end
                % changet this later %%%%%%%%%%%%%%%%%%%%%%%%%%%
                %% initialize the root gate
                gates(1).name = 'root';
                gates(1).ingate = ones(length(data), 1);  % every cell is in the root gate
                gates(1).parent = 0;
                gates(1).children = 0;    
                FCS.fieldinfo = fieldinfo;
            else              % gates already exist, find all objects within gate                
                FCS.fieldinfo = fieldinfo;
                for i = 2:length(gates)   % skip the parent gate, but look for all objects within the defined child gates   
                    parent = gates(i).parent;                    
                    bw = gates(i).bw;                    
                    if ~isempty(bw)
                        Xf = gates(i).Xf;
                        Yf = gates(i).Yf;
                        gates(i).ingate = (gates(parent).ingate).*(FCS.findingate(bw, Xf, Yf));
                    end
                    % re-assign gates                    
                end
            end            
            %% save all fcs class properties             
            
            FCS.gates = gates;            
            %% recompute gates if not yet specified
            
        end                
        %%%%% get functions        
        function names = fnames(FCS, varargin)
            names = fieldnames(FCS.data);            
            if (nargin > 1)   % pick specific fieldname
                names = names{varargin{1}};
            end            
        end
        
        function names = gnames(FCS, varargin)
            names = {FCS.gates.name};
            if (nargin > 1)   % pick specific gatename
                names = names{varargin{1}};
            end            
        end       
        
        function range = getrange(FCS, fieldnum)             
            range = FCS.fieldinfo(fieldnum).range;
        end
        
        function data = getdata(FCS, fieldnum, varargin)    % get data from a specific field                        
            name = FCS.fieldinfo(fieldnum).name;
            data = FCS.data;            
            data = [data.(name)]';
            if (nargin == 3)   % specific gate is specified 
                gatenum = varargin{1};
                inds = find(FCS.gates(gatenum).ingate);  % indices of cells that are in the gate                
                data = data(inds);
            end            
        end
        
        function in = ingate(FCS, gatenum)   % returns the vector containing the cells that are outside or inside a specified gate
            in = FCS.gates(gatenum).ingate;
        end
        
        function ingate = findingate(FCS, bw, Xf, Yf);  % determines whether cell is in gate or not            
            Xdata = FCS.getdata(Xf);
            Ydata = FCS.getdata(Yf);            
            
            xr = FCS.getrange(Xf);
            xmin = xr(1);
            xmax = xr(end);
            xB = length(xr);
            
            yr = FCS.getrange(Yf);
            ymin = yr(1);
            ymax = yr(end);
            yB = length(yr);
            
            Xbnum = floor((xB-1).*(Xdata - xmin)./(xmax-xmin))+1;   % bin numbers for each data point
            Ybnum = floor((yB-1).*(Ydata - ymin)./(ymax-ymin))+1;
            Xbnum(find(Xbnum<1)) = 1;
            Xbnum(find(Xbnum>(xB-1))) = xB-1;
            Ybnum(find(Ybnum<1)) = 1;
            Ybnum(find(Ybnum>(yB-1))) = yB-1;
            
            %% find cells within the gate
            mask = find(bw);
            linYX = sub2ind([yB-1,xB-1],Ybnum,Xbnum);
            ingate = ismember(linYX,mask);                                                            
        end
            
        function b = isgate(FCS, gatenum, Xf, Yf);       % verify whether the X and Y axis correspond to that define gate
            if (gatenum == 1)
                b = 0;
            else
                Xfg = FCS.gates(gatenum).Xf;
                Yfg = FCS.gates(gatenum).Yf;        
                if (Xfg == Xf)&(Yfg == Yf)   % correct plots are displayed
                    b = 1;
                else
                    b = 0;
                end
            end
        end
        
        %%%%%% set functions        
        function FCS = setrange(FCS, f, varargin)            
            if (nargin == 2)   % no option arguments, use dialog to get range
                range = FCS.getrange(f);            
                dmin = min(range);
                dmax = max(range);
                bins = length(range)-1;
                r = inputdlg(['Enter min (' num2str(dmin) '), max (' num2str(dmax) '), and bins (' num2str(bins)  ') for ' FCS.fnames(f) ':']);
                r = r{1};
                r = eval(r);
                range = r(1):(r(2)-r(1))./r(3):r(2);
            else   % range automatically passed to the 
                range = varargin{1};
            end
            FCS.fieldinfo(f).range = range;
        end
                   
        function FCS = newgate(FCS, name, bw, xi, yi, Xf, Yf, parent);       
            % creates new gate based on the provided information            
            %% load relevant information
            gates = FCS.gates;
                        
            %% return the X and Y bin numbers for all objects                        
            ingate = (FCS.ingate(parent)).*(FCS.findingate(bw, Xf, Yf));   % this is the indices of all objects that fall inside this gate and also only in the parent gate
            
            %% name of the gate
            gatename = [gates(parent).name ':' name];  % gatename, with parent identity appended            
            %% save gate information
            G = length(gates)+1;            
            gates(G).name = gatename;
            gates(G).ingate = ingate;
            gates(G).parent = parent;
            gates(parent).children = G;    % the parent gate now has a child
            gates(G).children = 0;
            gates(G).Xf = Xf;
            gates(G).Yf = Yf;
            gates(G).xi = xi;
            gates(G).yi = yi; 
            gates(G).bw = bw;
            FCS.gates = gates;                        
        end
        
        function FCS = newfield(FCS, name, dataf);                        
            data = FCS.data;
            
            if (length(dataf) ~= length(FCS.data))
                error('Existing and new data does not match...')
            end
                        
            for i = 1:length(dataf)
                data(i).(name) = dataf(i);    % add new field
            end
            
            bins = 100;     % this is the default number of bins
            dmin = min(dataf);    % find the minimum
            dmax = max(dataf);    % find the maximum            
            
            FCS.data = data;
            FCS.fieldinfo(end+1).name = name;            
            FCS.fieldinfo(end).range = dmin:(dmax-dmin)./bins:dmax;     % the display range                                
        end
    end
end