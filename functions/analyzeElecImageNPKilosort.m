function res = analyzeElecImageNPKilosort(experiment, st)
%ANALYZEELECIMAGE Summary of this function goes here
%   Detailed explanation goes here
%--------------------------------------------------------------------------
targetPath=[experiment.originalFolder,'\electrical_data'];
if not(exist(targetPath, 'dir')); mkdir(targetPath); end
saveName = 'elecimage_analysis.mat';
%--------------------------------------------------------------------------
Ncells = size(experiment.clusters,1);
ksdir = fullfile(experiment.originalpath, 'ks_sorted');
binpath = fullfile(ksdir, 'alldata.dat');

chanMap = 'C:\Users\karamanl\Documents\GitHub\KiloSortMEA\configFiles\neuropixPhase3B1_kilosortChanMap.mat'; 
ccmap   = load(chanMap);


Tmin = ceil(-1.5e-3*experiment.fs); Tmax = floor(2.5*1e-3*experiment.fs);
dt = Tmin:Tmax;
Nt = numel(dt);
NchanTOT = numel(ccmap.chanMap);


stimSums   = zeros([Ncells, NchanTOT, Nt], 'single');
stimSpikes = zeros([Ncells, 1], 'single');

Nskip     = 10;
batchsize = 10 * experiment.fs;
Nbatches  = ceil(sum(experiment.stimsamples)/batchsize);

spktimes = st(:,1);
spkids   = st(:,2);
cspksall = accumarray(spkids, spktimes, [Ncells 1], @(x) {x});

%--------------------------------------------------------------------------
loadnew = true;
if exist(fullfile(targetPath,saveName), 'file')
    res = load(fullfile(targetPath,saveName));
    if size(res.stimTemplatesMean, 1) == Ncells
        loadnew = false;
    end
end
%--------------------------------------------------------------------------
if loadnew
    %--------------------------------------------------------------------------
    fprintf('Extracting electrical images... '); tic;
    fid = fopen(binpath, 'r');

    msg = [];
    for ibatch = 1:Nskip:Nbatches

        
        mint = batchsize * (ibatch - 1);
        maxt = batchsize * ibatch;

        offset = 2 * NchanTOT * batchsize * (ibatch - 1);
        fseek(fid, offset, 'bof');
        dat = fread(fid, [NchanTOT batchsize], '*int16');

        for icell = 1:Ncells
            cellspikes = double(cspksall{icell});
            cellspikes((cellspikes+Tmin)<(mint+1) | (cellspikes+Tmax)>maxt) = [];
            cellspikes  = cellspikes - mint;
            if isempty(cellspikes), continue; end
            Nspikes = numel(cellspikes);

            spseek =  dt' + cellspikes';
            spkwvfrms = single(dat(:, spseek));
            spkwvfrms = reshape(spkwvfrms, NchanTOT, numel(dt), Nspikes);
            stimSums(icell, :, :) = squeeze(stimSums(icell, :, :)) + sum(spkwvfrms, 3);
            stimSpikes(icell)     = stimSpikes(icell)  + Nspikes;
        end
        %--------------------------------------------------------------------------
        fprintf(repmat('\b', 1, numel(msg)));
        msg = sprintf('Batch %d/%d. Time elapsed %2.2f s...\n', ibatch,Nbatches,toc);
        fprintf(msg);
        %--------------------------------------------------------------------------
    end
    fclose(fid);
    %--------------------------------------------------------------------------
    res.stimTemplatesMean   = stimSums./stimSpikes;
    res.stimSpikes          = stimSpikes;
    res.templateTimes = dt/experiment.fs;
    %--------------------------------------------------------------------------
    fprintf('Saving pre-analysis results... '); tic;
    save(fullfile(targetPath,saveName), '-struct', 'res')
    fprintf('Done! Took %2.2f s\n', toc);
    %--------------------------------------------------------------------------
end
%--------------------------------------------------------------------------

res.coords = [ccmap.xcoords ccmap.ycoords];
%--------------------------------------------------------------------------
% extract features for clustering later
[t2ptime, t2pval] = extractWaveformFeatures(res.stimTemplatesMean, dt/experiment.fs, 3);
res.t2ptime = t2ptime;
res.t2pval  = t2pval;

%%

dtcorr  = 5e-4;
Ncorr   =  1000e-3/dtcorr;
autoCorrs          = NaN(Ncells, Ncorr, 'single');
for icell = 1:Ncells
    currtrain = double(cspksall{icell})/experiment.fs;
    K = ccg(currtrain, currtrain, Ncorr, dtcorr);
    K(Ncorr+1) = 0;
    autoCorrs(icell, :) = single(K((Ncorr+1):end-1));
end

res.autoCorrs     = autoCorrs;
res.autoCorrsTime = ((0:Ncorr-1)+0.5) * dtcorr;
res.unitspikes    = cellfun(@numel, cspksall);
res.unitfr        = res.unitspikes * experiment.fs /sum(experiment.stimsamples);


[sptcv, sptlv]=spikeTrainCvLv(cspksall, experiment.fs);
res.spikecv  = sptcv;
res.spikelv  = sptlv;



% Xmat = [zscore(t2ptime), zscore(t2pval)];
% [idx, C] = kmeans(Xmat, 2, 'Replicates',400);
% cluscol = cbrewer('qual','Dark2',6);
% % subplot(1,4,1)
% % plot(res.templateTimes, squeeze(wvsclus(:,1,:)), 'Color',[0 0 0 0.2])
% subplot(1,4,2); cla;
% hold on;
% for ii = 1:size(C, 1)
%     plot(t2ptime(idx==ii,1)*1e3, t2pval(idx==ii,1),'o','Color', cluscol(ii,:))
% end
% % subplot(1,4,3);cla;
% % hold on;
% % for ii = 1:size(C, 1)
% %     plot(res.templateTimes, squeeze(wvsclus(idx==ii,1,:))', 'Color', [cluscol(ii,:) 0.3])
% % end
% autoccg = res.autoCorrs(:, 1:200);
% autoccg = autoccg./sum(autoccg, 2);
% subplot(1,4,4);cla;
% hold on;
% for ii = 1:size(C, 1)
% %     histogram(res.unitfr(idx==ii),linspace(0,50,40),'Normalization','probability',...
% %         'FaceColor',cluscol(ii,:))
%      %plot(res.autoCorrsTime(1:200), autoccg(idx==ii,:),'Color', [cluscol(ii,:) 0.2])
%      plot(res.autoCorrsTime(1:200), mean(autoccg(idx==ii,:)),'Color', cluscol(ii,:))
% end
%%
%----------------------------------------------------------------------
fprintf('Saving analysis results... '); tic;
save(fullfile(targetPath,saveName), '-struct', 'res')
fprintf('Done! Took %2.2f s\n', toc);
%--------------------------------------------------------------------------
end

