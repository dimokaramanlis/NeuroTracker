function matchmat = matchTemplatesPairwise(edata1, edata2)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

%==========================================================================
% extract features that we will use
resfeatures1 = getTemplateFeatures(edata1, edata1.coords);
resfeatures2 = getTemplateFeatures(edata2, edata2.coords);
Nunits1 = size(edata1.autoCorrs,1);
Nunits2 = size(edata2.autoCorrs,1);

%==========================================================================
% euclidean distance metric
Cdist = pdist2(resfeatures1.gausscents, resfeatures2.gausscents);
mind  = min(Cdist,[],'all');
Sdist = (100 - Cdist)./(100 - min(mind,[],'all'));
Sdist = Sdist .* (Cdist<100);
Sdist(Sdist>1) = 1;
%==========================================================================
speed1 = max([resfeatures1.speedbottom resfeatures1.speedtop], [], 2);
speed2 = max([resfeatures2.speedbottom resfeatures2.speedtop], [], 2);

feat1ori = [mean(resfeatures1.trepol, 2) mean(resfeatures1.t2ptime, 2) resfeatures1.axonval...
    speed1 resfeatures1.unitamps resfeatures1.bslfr resfeatures1.gausssigma edata1.spikecv edata1.spikelv];
feat2ori = [mean(resfeatures2.trepol, 2) mean(resfeatures2.t2ptime, 2) resfeatures2.axonval...
    speed2 resfeatures2.unitamps resfeatures2.bslfr resfeatures2.gausssigma edata2.spikecv edata2.spikelv];

% feat1ori = [mean(resfeatures1.trepol, 2) mean(resfeatures1.t2ptime, 2) resfeatures1.axonval...
%     speed1 resfeatures1.bslfr resfeatures1.gausssigma edata1.spikecv edata1.spikelv];
% feat2ori = [mean(resfeatures2.trepol, 2) mean(resfeatures2.t2ptime, 2) resfeatures2.axonval...
%     speed2 resfeatures2.bslfr resfeatures2.gausssigma edata2.spikecv edata2.spikelv];
%==========================================================================

feat2 = (feat2ori - min(feat1ori, [], 1))./quantile(feat1ori,0.99,1);
feat2(feat2>1) = 1;
feat2(feat2<0) = 0;
feat1 = (feat1ori - min(feat1ori, [], 1))./quantile(feat1ori,0.99,1);
feat1(feat1>1) = 1;
%==========================================================================
Srob    = robustcov([feat1ori;feat2ori]);
Snorm    = cov([feat1ori;feat2ori]);
distall = pdist2(feat1ori, feat2ori, 'mahalanobis',Snorm);
distallrob = pdist2(feat1ori, feat2ori, 'mahalanobis',Srob);

feattot = [feat1ori;feat2ori];
feattot = zscore(feattot);
bb = feattot*pcacov(Snorm./sqrt(diag(Snorm)));
distallcos = pdist2(bb(1:Nunits1,:), bb(Nunits1+1:end,:), 'cosine');
% asm = munkres(Cdistuse.*distallcos);
% asm = munkres(Cdistuse/100 +distallcos/max(distallcos,[],'all'));
%==========================================================================
Smult = Sdist;
Smult(Smult ==0) = 1e-6;
finalsim = Smult./distallcos;
[~,   womenpref] = sort(finalsim, 2, 'descend');
[~,   menpref] = sort(finalsim', 2, 'descend');

finalmatch = galeShapley(menpref,womenpref);
fcorr = nan(size(finalmatch,1), 1);
dfall = nan(size(finalmatch,1), 1);
for ii = 1:size(finalmatch,1)
    temp1 = squeeze(edata1.stimTemplatesMean(ii,:,:));
    temp2 = squeeze(edata2.stimTemplatesMean(finalmatch(ii,2),:,:));
%     temp1 = squeeze(edata2.stimTemplatesMean(ii,:,:));
%     temp2 = squeeze(edata1.stimTemplatesMean(finalmatch(ii,2),:,:));
    fcorr(ii) = maxSlidingCorr(temp1, temp2, 30);
    dfall(ii) = finalsim(ii, finalmatch(ii,2));
%     dfall(ii) = finalsim(finalmatch(ii,2),ii);
end


%%
featlist = {'unitamps','gausssigma','speedtop','speedbottom','t2ptime','trepol','bslfr','spikecv','spikelv',};
clf;
p = panel();
p.pack('v',2)
p(1).pack('h', numel(featlist))
for ii = 1:numel(featlist)
    f1 = resfeatures1.(featlist{ii})(finalmatch(:,1));
    f2 = resfeatures2.(featlist{ii})(finalmatch(:,2));
    maxf = max([f1;f2]);
    p(1,ii).select();
    scatter(f1, f2, 10, fcorr,'filled')
    line([0 maxf], [0 maxf])
    axis equal; xlim([0 maxf]); ylim([0 maxf])
    title(featlist{ii})
end

%%
iuse = 205;
clf;plotSpikeTemplate(squeeze(edata1.stimTemplatesMean(iuse,:,:))', edata1.coords, 'k')
plotSpikeTemplate(squeeze(edata2.stimTemplatesMean(finalmatch(iuse,2),:,:))', edata1.coords, 'r')
yc = resfeatures1.gausscents(iuse,2)';
ylim(yc + [- 200 200])
%%
%==========================================================================
end