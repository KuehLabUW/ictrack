function [z channel] = readmmheader(inputdir)

headerfile = [inputdir '/Pos_000_000/Acqusition.xml'];
fid = fopen(headerfile);

% regular expressions for capturing header information
pos_expr = 'Pos_(\d\d\d)_(\d\d\d)';  % regular expression
zbottom_expr =  '<entry key=\"acqZbottom\" value=\"(.+)\"/>';
ztop_expr = '<entry key=\"acqZtop\" value=\"(.+)\"/>';
zstep_expr = '<entry key=\"acqZstep\" value=\"(.+)\"/>';
cnum_expr = '<entry key=\"acqNumchannels\" value=\"(.+)\"/>';
cin_expr    =  '<entry key=\"acqChannelName\d\" value=\"(.+)\"/>';  % channel name
expose_expr =  '<entry key=\"acqChannelExp\d\" value=\"(.+)\"/>';   % exposure time
doz_expr    =  '<entry key=\"acqChannelDoZStack\d\" value=\"(.+)\"/>';  % do z-stack or not
zoffset_expr = '<entry key=\"acqChannelZOffset\d\" value=\"(.+)\"/>';   % z offset of channel
cR_expr = '<entry key=\"acqChannelColorR\d\" value=\"(.+)\"/>';    % red color setting (can be used to determine whether channel is fluorescent or not
cG_expr = '<entry key=\"acqChannelColorG\d\" value=\"(.+)\"/>';    % green color setting
cB_expr = '<entry key=\"acqChannelColorB\d\" value=\"(.+)\"/>';    % blue color setting
tskip_expr = '<entry key=\"acqSkip\d\" value=\"(.+)\"/>';

% read z stack information
z.bottom = str2num(fregexp(fid, zbottom_expr));
z.top = str2num(fregexp(fid, ztop_expr));
z.step = str2num(fregexp(fid, zstep_expr));
cnums = str2num(fregexp(fid, cnum_expr));

for i = 1:cnums   % cycle through all the channels
    channel(i).name = fregexp(fid, cin_expr);
    channel(i).expose = str2num(fregexp(fid, expose_expr));
    
    doz = fregexp(fid, doz_expr);
    if strmatch(doz,'false')
        channel(i).doz = 0;
    elseif strmatch(doz,'true')
        channel(i).doz = 1;
    end
    
    channel(i).zoffset = str2num(fregexp(fid, zoffset_expr));
    channel(i).cR = str2num(fregexp(fid, cR_expr));
    channel(i).cG = str2num(fregexp(fid, cG_expr));
    channel(i).cB = str2num(fregexp(fid, cB_expr));
    channel(i).tskip = str2num(fregexp(fid, tskip_expr));
end
fclose(fid);

function token = fregexp(fid, expr)
tline = 0;
tok = [];
while(isempty(tok)&(tline ~= -1))
    tline = fgetl(fid);
    tok = regexp(tline, expr, 'tokens');
end

if (tline == -1)
    error(['end of file reached, string ' expr ' not found\n']);
end

token = tok{1}{1};


