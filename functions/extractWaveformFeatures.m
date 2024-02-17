function [t2ptime, t2pval, wvsclus,axonval, trepol] = extractWaveformFeatures(stimTemplatesMean, dt, NcloseChan)
%EXTRACTWAVEFORMFEATURES Summary of this function goes here
%   Detailed explanation goes here

[Ncells, NchanTOT, Nt] = size(stimTemplatesMean);
dtfit = -2:2;

% find three top channels
minpt      = min(stimTemplatesMean, [], 3);
[~, isort] = sort(minpt, 2, 'ascend');
%--------------------------------------------------------------------------
% initialize collections
wvsclus = zeros([Ncells, NcloseChan, Nt], 'single');
t2ptime = zeros(Ncells, NcloseChan);
t2pval  = zeros(Ncells, NcloseChan);
trepol  = zeros(Ncells, NcloseChan);
axonval = zeros(Ncells, NcloseChan);
%--------------------------------------------------------------------------
for icell = 1:Ncells

    currtemp =  squeeze(stimTemplatesMean(icell, isort(icell, 1:NcloseChan), :));
    currtemp = currtemp./sqrt(sum(currtemp.^2,'all'));
    wvsclus(icell, :, :) = currtemp;
    [~, imintemp] = min(currtemp,[],2);
    ttroughpeak   = zeros(NcloseChan, 2);
    valtroughpeak = zeros(NcloseChan, 2);
    trepolarize   = zeros(NcloseChan,1);

    for ii = 1:size(wvsclus, 2)
        %--------------------------------------------------------------------------
        % first peak prominence
        tempinit = currtemp(ii, 1:imintemp(ii));
        normfac = max(abs(max(tempinit)), abs(min(tempinit)));
        axonval(icell, ii) = range(tempinit)./sqrt(sum(currtemp(ii,:).^2));
        %--------------------------------------------------------------------------
        % trough extraction
        pfit = polyfit(dt(imintemp(ii)+dtfit), currtemp(ii, imintemp(ii)+dtfit), 2);
        ttroughpeak(ii, 1)   = -pfit(2)/(2*pfit(1));
        valtroughpeak(ii, 1) = pfit(3)-pfit(2)^2/(4*pfit(1));
        %--------------------------------------------------------------------------
        % peak extraction
        [~, isortpeak] = sort(currtemp(ii,:), 2, 'descend');
        tmax = isortpeak(findfirst(isortpeak-imintemp(ii)>0,2));
        induse = tmax+dtfit;
        induse(induse>Nt) = [];
        pfit2 = polyfit(dt(induse), currtemp(ii, induse), 2);
        ttroughpeak(ii, 2)   = -pfit2(2)/(2*pfit2(1));
        valtroughpeak(ii, 2) = pfit2(3)-pfit2(2)^2/(4*pfit2(1));
        %--------------------------------------------------------------------------
        % repolarization
        midrep = currtemp(ii, tmax)/2;
        [~, iminr] = min(abs(currtemp(ii, tmax:end) -  midrep));
        induse = tmax+iminr+dtfit;
        induse(induse>Nt) = [];
        pfit3 = polyfit(dt(induse), currtemp(ii, induse), 1);
        tsearch = linspace(min(dt(induse)), max(dt(induse)));
        [~, iminr] = min(abs(polyval(pfit3,tsearch)-midrep));
        trepolarize(ii) = tsearch(iminr);
        %--------------------------------------------------------------------------
    end
    t2ptime(icell, :) = diff(ttroughpeak, [],2);
    valtroughpeak = abs(valtroughpeak);
    t2pval(icell, :)  = (valtroughpeak(:, 1)-valtroughpeak(:,2))./ (valtroughpeak(:, 1)+valtroughpeak(:,2));
    trepol(icell, :)  = trepolarize - ttroughpeak(:,2);
end

end

