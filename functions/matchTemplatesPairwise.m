function [matchmat, fcorr, finalsim, featsori] = matchTemplatesPairwise(edata1, edata2)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

%==========================================================================
Nunits = [size(edata1.autoCorrs,1) size(edata2.autoCorrs,1)];
coords = edata1.coords;
%==========================================================================
% extract features that we will use
fprintf('Extracting template features... ');tic;
for isess = 1:2
    if isess == 1
        edata = edata1;
    else
        edata = edata2;
    end
    tempfeatures(isess) = getTemplateFeatures(edata, coords);
end
fprintf('Done! Took %2.2f\n', toc);
%==========================================================================
% euclidean distance metric
Cdist = pdist2(tempfeatures(1).gausscents, tempfeatures(2).gausscents);
mind  = min(Cdist,[],'all');
Sdist = (100 - Cdist)./(100 - min(mind,[],'all'));
Sdist = Sdist .* (Cdist<100);
Sdist(Sdist>1) = 1;
%==========================================================================
% construct featmat
propspeed = max([cat(1, tempfeatures(:).speedbottom) cat(1, tempfeatures(:).speedtop)],[],2);
axonval   = cat(1, tempfeatures(:).axonval);
t2ptime   = median(cat(1, tempfeatures(:).t2ptime), 2);
trepol    = median(cat(1, tempfeatures(:).trepol), 2);
spikecv   = cat(1, tempfeatures(:).spikecv);
spikelv   = cat(1, tempfeatures(:).spikelv);
bslfrate  = cat(1, tempfeatures(:).bslfr);
gsigma    = cat(1, tempfeatures(:).gausssigma);
unitamps  = cat(1, tempfeatures(:).unitamps);

featsori = [trepol t2ptime axonval gsigma unitamps propspeed bslfrate spikecv spikelv];
featsquant  = (featsori - min(featsori, [], 1))./quantile(featsori,0.99,1);
featsquant(featsquant>1) = 1;
featsquant(featsquant<0) = 0;
%==========================================================================
% define distances and compute final similarity metric
Srob       = robustcov(featsori,'Method','olivehawkins');
Snorm      = cov(featsori);
distall    = pdist2(featsori(1:Nunits(1), :), featsori(Nunits(1)+1:end, :), 'mahalanobis', Snorm);
distallrob = pdist2(featsori(1:Nunits(1), :), featsori(Nunits(1)+1:end, :), 'mahalanobis', Srob);

Smult = Sdist;
Smult(Smult ==0) = 1e-3;
finalsim = Smult./distallrob;
%==========================================================================
% perform stable matching 
fprintf('Performing stable matching... ');tic;
if Nunits(1) > Nunits(2)
    [~,  prefProp] = sort(finalsim, 1, 'descend');
    prefProp = prefProp';
    [~,   prefAcc]  = sort(finalsim, 2, 'descend');
    toperm  = true;
else
    [~,  prefProp] = sort(finalsim, 2, 'descend');
    [~,   prefAcc]  = sort(finalsim, 1, 'descend');
    prefAcc = prefAcc';
    toperm  = false;
end
finalmatch = flexibleGS(prefProp, prefAcc);

if toperm
    finalmatch = finalmatch(:, [2 1]);
end
matchmat = finalmatch;
fprintf('Done! Took %2.2f\n', toc);
%==========================================================================
% compute metrics
fprintf('Computing evaluation metrics... ');tic;

fcorr = nan(size(finalmatch,1), 1);
dfall = nan(size(finalmatch,1), 1);
for ii = 1:size(finalmatch,1)
    temp1 = squeeze(edata1.stimTemplatesMean(finalmatch(ii,1),:,:));
    temp2 = squeeze(edata2.stimTemplatesMean(finalmatch(ii,2),:,:));
   

%     temp1 = squeeze(edata2.stimTemplatesMean(ii,:,:));
%     temp2 = squeeze(edata1.stimTemplatesMean(finalmatch(ii,2),:,:));
    fcorr(ii) = maxSlidingCorr( denoiseTemplate(temp1, 3),  denoiseTemplate(temp2, 3), 40);
    dfall(ii) = finalsim(finalmatch(ii,1), finalmatch(ii,2));
end
fprintf('Done! Took %2.2f\n', toc);
%%
%==========================================================================
end