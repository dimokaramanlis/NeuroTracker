function [experiment] = loadExperimentNPKilosort(experiment, namestoread)
%LOADEXPERIMENTKILOSORT Loads the raw data for all stimuli in a MEA experiment
%
% This function loads everthing about on expriment folder, this include all
% parameters of stimulus, goodcells and clusters, frametimings, spiketimings
% and binned spikes.

ksdir = fullfile(experiment.originalpath,'ks_sorted');
sesslistpath   = 'D:\DATA_folder\Sessions\session_list.txt';
%==========================================================================
fprintf('Loading Kilosort data... '); tic;
spike_clusters  = readNPY(fullfile(ksdir, 'spike_clusters.npy'));
spike_times     = readNPY(fullfile(ksdir, 'spike_times.npy'));
spike_templates = readNPY(fullfile(ksdir, 'spike_templates.npy'));
templates       = readNPY(fullfile(ksdir, 'templates.npy'));
chanmap         = readNPY(fullfile(ksdir, 'channel_map.npy')) + 1;
fprintf('Done! Took %2.2fs\n',toc);
%==========================================================================
% reduce spikes to include only the good ones
fprintf('Getting good cells...\n'); tic;

[cids, cvals] = readClusterGroupsCSV(fullfile(ksdir, 'cluster_group.tsv'));
goodCellIds = cids(cvals == 2); Ncells = numel(goodCellIds);
goodCellIds = cast(goodCellIds,'like', spike_clusters);

fprintf('Found %d\n', Ncells);
stori = [double(spike_times), double(spike_clusters)];

% read ratings
clustratings = NaN(Ncells, 1);
qualfile = fullfile(ksdir, 'cluster_quality.tsv');
if exist(qualfile,'file')
    [qids, quals] = readClusterGroupsCSV(qualfile);
    for ientry = 1:numel(qids)
        clustratings(goodCellIds == qids(ientry)) = quals(ientry);
    end
end

% read comments
clustcomments = cell(Ncells, 1);
clustcomments(:) = {''};
commfile = fullfile(ksdir, 'cluster_comment.tsv');
if exist(qualfile,'file')
    [comids, comms] = readClusterCommentsCSV(commfile);
    for ientry = 1:numel(comids)
        if nnz(goodCellIds == comids(ientry))>0
            clustcomments{goodCellIds == comids(ientry)} = comms{ientry};
        end
    end
end
stuse = ismembc(spike_clusters,goodCellIds);
spike_times     = spike_times     (stuse);
spike_templates = spike_templates (stuse);
spike_clusters  = spike_clusters  (stuse);

[clustids, ~, new_clusters] = unique(spike_clusters);

[spike_times, avg_temps] = fixMisalignedUnits(spike_times, new_clusters,...
    spike_templates, templates, clustcomments);

st_cellarray = accumarray(new_clusters, spike_times,[], @(x) {x});

fprintf('Done! Took %2.2fs\n',toc);
%==========================================================================
% defining the function for separating spike vectors into clusters.
disp('Assigning templates and cell Ids...'); tic;

tempchannels = squeeze(min(avg_temps,[],2));
[~, imin] = min(tempchannels,[],2);
channasignments = chanmap(imin);

[channasignments, isort] = sort(channasignments,'ascend');

classignments = ones(Ncells, 1, 'single');
allds = diff(channasignments);
for icell = 2:Ncells
    if allds(icell-1) == 0
        classignments(icell) = classignments(icell-1) + 1;
    end
end

% set names
experiment.clusters = [clustids(isort) channasignments classignments clustratings(isort)];
experiment.comments = clustcomments(isort);

% re-sort data
st_cellarray = st_cellarray(isort);
sc_cellarray = cellfun(@(x,idx) idx*ones(numel(x),1, class(spike_times)),...
    st_cellarray, num2cell(1:Ncells)','un',false);
st = [cell2mat(st_cellarray) cell2mat(sc_cellarray)];
st = sortrows(st, 1);


fprintf('Done! Took %2.2fs\n',toc);
%==========================================================================
% load stimulus spliting info
bfile         = load(fullfile(ksdir,'samplelist.mat'));
experiment.stimsamples = bfile.samplelist;
experiment.stimpaths   = namestoread;
experiment.sessionids  = namestoread;
%==========================================================================
disp('Calculating electrical images and autocorrelations...'); tic;
%--------------------------------------------------------------------------
experiment.fs = 30000;
analyzeElecImageNPKilosort(experiment, st);
% plotElecImageNPKilosort()
%==========================================================================
disp('Saving stimulus parameters, pulses, rasters and binning...'); tic;
%--------------------------------------------------------------------------
targetPath = [experiment.originalFolder,'\session_data'];
if not(exist(targetPath, 'dir')); mkdir(targetPath);end
%--------------------------------------------------------------------------
Nrecordings    = numel(experiment.stimsamples);

msg = [];
csamples = [0; cumsum(experiment.stimsamples)];
csamples = cast(csamples, 'like', st);

for istim = 1 : Nrecordings
    %--------------------------------------------------------------------------
    stimdata.originalbin = experiment.stimpaths{istim};
    %--------------------------------------------------------------------------
    % Get rasters
    stimSt = st((spike_times>=csamples(istim)+1) & (spike_times<csamples(istim+1)+1),:);
    stimSt(:, 1) = stimSt(:, 1) - csamples(istim) + 1; % +1 for python
    stimdata.spiketimes = double(stimSt);
    %--------------------------------------------------------------------------
    % Prepare frametimings
    [aa, bb]    = fileparts(experiment.stimpaths{istim});
    binpath     = dir(fullfile(aa,'..','*.bin'));
    metapath    = dir(fullfile(aa,'..','*.meta'));
    metapathrec = dir(fullfile(experiment.stimpaths{istim},'..','*.ap.meta'));
    metadata    = readINI(fullfile(metapath.folder, metapath.name));
    metadatarec = readINI(fullfile(metapathrec.folder, metapathrec.name));
    stimdata.analogfs =  metadata.niSampRate;
    stimdata.apdatafs =  metadatarec.imSampRate;
    %--------------------------------------------------------------------------
    stimdata.sessid = extractSessIdNeuropixels(stimdata.originalbin, metadata);
    experiment.sessionids{istim} = stimdata.sessid;
    %--------------------------------------------------------------------------
    fid = fopen(fullfile(binpath.folder, binpath.name), 'r');
    data = fread(fid, [metadata.nSavedChans Inf], 'int16');
    fclose(fid);
    camerapulses = extractSynchTimings(single(data(2, :)'), metadata.niSampRate, 3);
    screenpulses = extractSynchTimings(single(data(3, :)'), metadata.niSampRate, 50);
    if contains(metapath.folder,'VisualStim')
        stimtype = 'visualstim';
    else
        stimtype = 'behavior';
    end
    stimdata.stimtype = stimtype;

    camerapulses.fonsets  = round(camerapulses.fonsets*stimdata.apdatafs/stimdata.analogfs);
    camerapulses.foffsets = round(camerapulses.fonsets*stimdata.apdatafs/stimdata.analogfs);

    screenpulses.fonsets  = round(screenpulses.fonsets*stimdata.apdatafs/stimdata.analogfs);
    screenpulses.foffsets = round(screenpulses.fonsets*stimdata.apdatafs/stimdata.analogfs);
    
    stimdata.camerapulses = camerapulses;
    stimdata.screenpulses = screenpulses;

% 
%     sons = screenTimingsFromPulses(screenpulses, stimtype);
%     for ii =1:numel(sons)
%         stimdata.stimonsets{ii} = round(sons{ii}*stimdata.apdatafs/stimdata.analogfs);
%     end

    %Save data
    sname = fullfile( targetPath, [num2str(istim) '_raw_data.mat']);
    save(sname,'-struct', 'stimdata');
    %--------------------------------------------------------------------------
    fprintf(repmat('\b', 1, numel(msg)));
    msg = sprintf('Recording %d/%d. Time elapsed %2.2f s...\n', istim, Nrecordings, toc);
    fprintf(msg);
    %--------------------------------------------------------------------------
end
disp('Done with loading!')




%     
%     stimdata = struct();
%     %--------------------------------------------------------------------------
%     %Get stimulus params
%     stimtext = fullfile(experiment.originalFolder, 'stimuli', stimtextpaths{istim});
%     stimdata.stimPara = readStimulusParameters(stimtext,...
%         experiment.projector.screen(2), experiment.projector.screen(1)); 
%     stimNames{istim} = stimdata.stimPara.stimulus;
%     %--------------------------------------------------------------------------
%     %Get frametimes
%     frametimetext = fullfile(experiment.originalFolder, 'frametimes', frametimestextpaths{istim});
%     stimuliFrames = load(frametimetext, 'fonsets','foffsets');
% 
%     stimdata.fonsets  = stimuliFrames.fonsets  + delsamples; 
%     stimdata.foffsets = stimuliFrames.foffsets + delsamples; 
%     %--------------------------------------------------------------------------
%     %Get rasters
%     stimSt = st((spike_times>=csamples(istim)+1) & (spike_times<csamples(istim+1)+1),:);
%     stimSt(:, 1) = stimSt(:, 1) - csamples(istim) + 1; % +1 for python
%     stimdata.spiketimes = double(stimSt);
%     %--------------------------------------------------------------------------
%     %Do initial binning
%     stimdata.spikesbin = loadBinnedSpikesKS(...
%         stimSt, Ncells, stimdata.fonsets, stimdata.foffsets, stimdata.stimPara);
%     %--------------------------------------------------------------------------
%     %Create folder for saving
%     folderName = fullfile(targetPath, [num2str(istim) '_' stimdata.stimPara.stimulus]);
%     if not(exist(folderName, 'dir')); mkdir(folderName);end
%     %Save data
%     sname = fullfile( folderName, [num2str(istim) '_raw_data.mat']);
%     save(sname,'-struct', 'stimdata');
%    %--------------------------------------------------------------------------
%     fprintf(repmat('\b', 1, numel(msg)));
%     msg = sprintf('Stimulus %d/%d. Time elapsed %2.2f s...\n', istim,Nrecordings,toc);
%     fprintf(msg);
% end
% experiment.stimnames = stimNames;
%==========================================================================


end
