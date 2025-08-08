function samplelist = concatenateAndCleanBinaries(dpfile, dptarget, varargin)

%dpfile = {'', ''}; % fill this with the paths you want 
%dptarget = 'finalbin.bin';

if nargin <3
    Ncomps = 3;
else
    Ncomps = varargin{1};
end

[targetfolder, ~] = fileparts(dptarget);

fidtarget    = fopen(dptarget, 'W');
NchanTOT     = 385;
batchsamples = 0.5*30000;
ichuse       = 1:384;
Nuse         = numel(ichuse);
samplelist   = zeros(numel(dpfile), 1); % samplelist is important to tell us how to split files later
scaleproc    = 200;
%-----------------------------------------------------------------------------
for ii = 1:numel(dpfile)
    [fpath, ~] = fileparts(dpfile{ii});
    metapath   = dir(fullfile(fpath, '*.ap.meta'));
    res(ii) = SGLXMetaToCoords(fullfile(metapath.folder,metapath.name));
end

if numel(dpfile) > 1
    assert(isequal(res(:).xcoords) & isequal(res(:).ycoords), "Shanks don't match!!!")
end
%-----------------------------------------------------------------------------
[fpath, ~] = fileparts(dpfile{1});
metapath   = dir(fullfile(fpath, '*.ap.meta'));
res        = SGLXMetaToCoords(fullfile(metapath.folder,metapath.name), targetfolder);
%-----------------------------------------------------------------------------
[foldertop, ~] = fileparts(targetfolder);
%-----------------------------------------------------------------------------
ops.fs = 30000; 
ops.fshigh = 300;

%-----------------------------------------------------------------------------
tic;
for ifile = 1:numel(dpfile)
    fprintf('Going through file %d/%d...\n', ifile, numel(dpfile));
    s = dir(dpfile{ifile});         

    bytes =  s.bytes; % size in bytes of raw binary
    nTimepoints = floor(bytes/NchanTOT/2); % number of total timepoints
    Nbatch      = ceil(nTimepoints/batchsamples);
    
    synchpulsename  = sprintf('%02d_synch_pulse.dat', ifile);
    fidsynch        = fopen(fullfile(foldertop, synchpulsename), 'W');

    fid = fopen(dpfile{ifile},'r');  msg = [];
    for ibatch = 1:Nbatch
        
        sampsread = min(batchsamples, nTimepoints - (ibatch-1)*batchsamples);
        fseek(fid, (ibatch-1)*batchsamples*NchanTOT*2, 'bof'); % fseek to batch start in raw file
        dat = fread(fid, [NchanTOT sampsread], '*int16');

        datsynch = dat(end, :);
        dataRAW = gpuArray(dat(ichuse,:));
        dataRAW = dataRAW';
        dataRAW = single(dataRAW)/scaleproc;

        dataRAW = gpufilter(dataRAW, ops);

        % SVD model for noise
        [aa,bb,cc] = svd(dataRAW,'econ');
        toremove   = aa(:,1:Ncomps)*bb(1:Ncomps,1:Ncomps)*cc(:,1:Ncomps)';
        dataDEN    = dataRAW - toremove;

        datcpu  = gather(int16(dataDEN*scaleproc)); % convert to int16, and gather on the CPU side

        fwrite(fidtarget, datcpu', 'int16'); % write main chans
        fwrite(fidsynch, datsynch', 'int16'); % write synch file

        %--------------------------------------------------------------------------
        fprintf(repmat('\b', 1, numel(msg)));
        msg = sprintf('Batch %d/%d. Time %d s \n', ibatch, Nbatch, round(toc));
        fprintf(msg);
        %--------------------------------------------------------------------------
    end
    fclose(fid);
    samplelist(ifile) = nTimepoints;
end
fclose(fidtarget);
%--------------------------------------------------------------------------
% save samplelist
save(fullfile(targetfolder, 'samplelist.mat'), 'samplelist','dpfile');
% save event markers
stimtimes = [0;cumsum(samplelist(1:end-1))];
fidOut= fopen([targetfolder,filesep,'eventmarkers.txt'], 'W');
fprintf(fidOut,'%d\r\n',stimtimes);
fclose(fidOut);
% save event marker names
fidOut= fopen([targetfolder,filesep,'eventmarkernames.txt'], 'W');
[~, stimnames, ~] = fileparts(dpfile');
fprintf(fidOut,'%s\r\n',stimnames{:});
fclose(fidOut);
%--------------------------------------------------------------------------
end


