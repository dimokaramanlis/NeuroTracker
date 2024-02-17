function resfeatures = getTemplateFeatures(edata, coords)
%GETTEMPLATEFEATURES Summary of this function goes here
%   Detailed explanation goes here

%--------------------------------------------------------------------------
% copy properties from edata
resfeatures = struct();
resfeatures.spikecv = edata.spikecv;
resfeatures.spikelv = edata.spikelv;
[t2ptime, t2pval, ~, axonval,trepol] = ...
    extractWaveformFeatures(edata.stimTemplatesMean, edata.templateTimes, 3);
resfeatures.t2pval  = t2pval;
resfeatures.t2ptime = t2ptime;
resfeatures.trepol  = trepol;
%--------------------------------------------------------------------------
% extract baseline firing rate from ACG
dacg              = median(diff(edata.autoCorrsTime));
allacg            = edata.autoCorrs./edata.unitspikes/dacg;
resfeatures.bslfr = mean(allacg(:,end-10:end), 2);
%--------------------------------------------------------------------------
tfs = median(diff(edata.templateTimes));
templatemat = edata.stimTemplatesMean;
dt     = -7:15;
dtsmall = -3:3;
maxRexp = 200;
Nsitesprop = 13;
Ncompsdenoise = 3;
dttrough    = -20:40;
[Nunits, Nchan, Nt] = size(templatemat);
templatemat         = templatemat - median(templatemat,[2 3]);

twts      = range(templatemat, 3);
%--------------------------------------------------------------------------
% initialize collections
centsall       = NaN(Nunits, 2);
avgwfmall      = NaN(Nunits, Nt);
decaysall      = NaN(Nunits, 1);
centsallgauss  = NaN(Nunits, 2);
sigmaall       = NaN(Nunits, 1);
tfitprofile    = NaN(Nunits, 27);
cfitprofile    = NaN(Nunits, 27);
ampsall        = NaN(Nunits, 1);

speedtop       = zeros(Nunits, 1);
speedbot       = zeros(Nunits, 1);
axonvalues     = zeros(Nunits, 1);

%--------------------------------------------------------------------------


for iunit = 1:Nunits

    utemp       = squeeze(templatemat(iunit,:,:));
    [~, imsite]  = max(twts(iunit,:));
    [~, itrough] = min(utemp(imsite,:));
%     cwts  = max(abs(utemp(:, dtcell)), [], 2);
%     cwts  = cwts - quantile(cwts,0.05);
    allds = sqrt(sum((coords - coords(imsite,:)).^2,2));

    cwts  = twts(iunit, :)';
    cwfit = cwts(allds<maxRexp);
    pfit  = polyfit(allds(allds<maxRexp), log(cwfit), 1);
    %---------------------------------------------------------------------
    % extract decay-based center and avg waveform

    xmax     = log(10)/abs(pfit(1));
    iptsuse  = allds<xmax;
    currcent = cwts(iptsuse)'*coords(iptsuse,:)/sum(cwts(iptsuse));

    alldsuse = sqrt(sum((coords(iptsuse,:) - currcent).^2,2));
    wvfwts   = 1 - alldsuse/xmax;
    wvfwts   = wvfwts.*(wvfwts>0);
    wvfwts   = wvfwts/sum(wvfwts);
    avgwfm   = wvfwts'*utemp(iptsuse, :);

    ampsall(iunit) = wvfwts' * twts(iptsuse);

    centsall(iunit, :) = currcent;
    avgwfmall(iunit,:) = avgwfm;
    decaysall(iunit)   = mean((cwts(imsite)-cwts(iptsuse))./allds(iptsuse), 'omitnan');
    %---------------------------------------------------------------------
    % extract Gaussian parameters

    [~, isortd] = sort(alldsuse,'descend');
    imgfit      = cwts(iptsuse);
    imgfit      = imgfit - mean(imgfit(isortd(1:6)));
    imgfit      = imgfit/max(imgfit);
    fitparams   = fitCircGaussForMEA(coords(iptsuse,:), double(imgfit));

    centsallgauss(iunit, :) = fitparams(1:2);
    sigmaall(     iunit, :) = fitparams(3);
    fullparams  = [fitparams(1:2) fitparams([3 3]) 0 1];
    cel = getEllipseFromParams(fullparams, 2, 50);
    sitesel = inpolygon(coords(:,1), coords(:,2), cel(1,:), cel(2,:));
    elsiteidx = find(sitesel);
    axvals    = NaN(size(elsiteidx));
    for isite = 1:nnz(sitesel)
        sctemp = utemp(elsiteidx(isite), :);
        [~, iminsc] = min(sctemp);
        normfac = abs(min(sctemp(1:iminsc)));
        axvals(isite) = range(sctemp(1:iminsc))/normfac - 1;
    end
    axonvalues(iunit) = cwts(elsiteidx)'*axvals/sum(cwts(elsiteidx));
    %---------------------------------------------------------------------
    % extract spike propagation
    sitesuse = imsite + (-Nsitesprop:Nsitesprop);
    sitesuse(sitesuse <     1) = [];
    sitesuse(sitesuse > Nchan) = [];
    indsave  = sitesuse;


    timesfit  = NaN(numel(sitesuse), 1);
    valsfit   = zeros(numel(sitesuse), 1);
    wvfsearch = utemp(sitesuse, :);
    wvfsearch = wvfsearch - median(wvfsearch(:, 1:5), 2);
    dtsearch  = itrough + dttrough;

    wvfsmall  = wvfsearch(:, dtsearch);
    
    Ntroughdist = 12;

    % svd model of waveform
    [aa, bb, cc] = svd(wvfsmall);
    svdmodel =  aa(:,1:Ncompsdenoise)*bb(1:Ncompsdenoise,1:Ncompsdenoise)*cc(:,1:Ncompsdenoise)';
%     clf;      
%     subplot(1,2,1)
%     plot(wvfsmall' + (1:numel(sitesuse))*10)
%     subplot(1,2,2)
%     plot(svdmodel' + (1:numel(sitesuse))*10)
%     
    for ipt = 1:numel(sitesuse)
        wvfmcurr = svdmodel(ipt,:);
        %wvfmcurr = wvfmcurr - mean(wvfmcurr([1:3 end-3:end]));
        % wvfmcurr = highpass(wvfmcurr,1000,30000);
        % wvfmcurr = wvfmcurr - movmedian(wvfmcurr,40);
        [minval,   imincurr] = min(wvfmcurr);
        [negvals, inegpeaks] = findpeaks(-wvfmcurr,'NPeaks',3,'SortStr','descend','MinPeakDistance',6,'MinPeakProminence',0.1);
        [posval,   ipospeak] = findpeaks(wvfmcurr,'NPeaks',1,'SortStr','descend');
        
        inegpeaks(inegpeaks < 6 | inegpeaks > numel(wvfmcurr)-6 ) = [];
        if ~isempty(inegpeaks)
            imincurr = inegpeaks(1);
        end
%         % replace minimum point 
%         if abs(minval)<posval
%             iafter      = (inegpeaks - ipospeak)>0;
%             if sum(iafter)==1
%                 imincurr =  inegpeaks(iafter);
%             else
% %                 [~, ifar] = max(abs((inegpeaks - ipospeak)));
% %                 iafter = true(size(iafter));
% %                 iafter(ifar) = 0;
%             end
%         end

        timesfit(ipt) = imincurr;
        valsfit(ipt)  = svdmodel(ipt,imincurr);

%         indstime = imincurr+dtsmall;
%         indstime(indstime < 1 | indstime>numel(wvfmcurr)) = [];
%         pfit = polyfit(indstime, wvfmcurr(indstime), 2);
%         timesfit(ipt) = -pfit(2)/(2*pfit(1));
%         valsfit(ipt)  = -pfit(2)^2/(4*pfit(1)) + pfit(3);
    end
%   
    %---------------------------------------------------------------------
    % triage propagation
%     line(timesfit, valsfit +(1:numel(sitesuse))'*10,...
%         'LineStyle','none','Color','k','Marker','o')

    cfit    = coords(sitesuse,2) - coords(imsite,2);
    cfitprofile(iunit, indsave) = cfit;
    fitpoly2=fit(cfit,timesfit,'poly2', 'Weights',range(wvfsmall,2));

    tpred = feval(fitpoly2, cfit);
    
    properrs = timesfit - tpred;
    irem     = abs(properrs) > mad(properrs, 1) * 1.4 * 6;
    timesfit(irem) = nan;
%     line(timesfit, valsfit +(1:numel(sitesuse))'*10,...
%         'LineStyle','none','Color','r','Marker','o')
    tfitprofile(iunit, indsave) = timesfit;
    %---------------------------------------------------------------------
    [cfit, ~, icun] = unique(cfit);
    tfit = accumarray(icun, timesfit,[],@nanmean);
    
    ispeedtop = 0;
    iuse = cfit>=0 & ~isnan(tfit);
    if nnz(iuse) > 4
        coefftop = polyfit(cfit(iuse), tfit(iuse), 1);
        ispeedtop = abs(coefftop(1) * tfs)*1e6;
    end
    ispeedbot = 0;
    iuse = cfit<=0 & ~isnan(tfit);
    if  nnz(iuse)  > 4
        coeffbot = polyfit(cfit(iuse), tfit(iuse), 1);
        ispeedbot = abs(coeffbot(1) * tfs)*1e6;
    end
    
    speedtop(iunit) = ispeedtop;
    speedbot(iunit) = ispeedbot;
    %---------------------------------------------------------------------
% 
end
resfeatures.axonval = axonvalues;

resfeatures.unitcents = centsall;
resfeatures.avgwforms = avgwfmall;
resfeatures.wfmdecays = decaysall;
resfeatures.gausscents = centsallgauss;
resfeatures.gausssigma = sigmaall;
resfeatures.tfitprofile = tfitprofile;
resfeatures.cfitprofile = cfitprofile;
resfeatures.speedtop    = speedtop;
resfeatures.speedbottom = speedbot;
resfeatures.unitamps    = ampsall;
% resfeatures.wfmamps   = max(abs(avgwfmall), [], 2);

end

