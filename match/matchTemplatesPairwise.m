function [matchmat, fcorr, finalsim, featsori] = matchTemplatesPairwise(feats1, feats2)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
%==========================================================================
tempfeatures = cat(1, feats1, feats2);
Nunits       = [size(feats1.bslfr,1), size(feats2.bslfr,1)];
%==========================================================================
% euclidean distance metric
Cdist = pdist2(feats1.gausscents, feats2.gausscents);

mind  = min(Cdist,[],'all');
Sdist = (100 - Cdist)./(100 - min(mind,[],'all'));
Sdist = Sdist .* (Cdist<100);
Sdist(Sdist>1) = 1;
%==========================================================================
% wvformdist = 1 - pdist2(feats1.avgwforms, feats2.avgwforms,@(x1,x2) barcodedist(x1,x2,20));
%==========================================================================
% construct featmat
propspeed = max([cat(1, tempfeatures(:).speedbottom) cat(1, tempfeatures(:).speedtop)],[],2);
axonval   = cat(1, tempfeatures(:).axonval);
t2ptime   = cat(1, tempfeatures(:).t2ptime);
trepol    = cat(1, tempfeatures(:).trepol);
spikecv   = cat(1, tempfeatures(:).spikecv);
spikelv   = cat(1, tempfeatures(:).spikelv);
bslfrate  = cat(1, tempfeatures(:).bslfr);
gsigma    = cat(1, tempfeatures(:).gausssigma);
unitamps  = cat(1, tempfeatures(:).unitamps);
t2pval    = cat(1, tempfeatures(:).t2pval);
wvfdecays = cat(1, tempfeatures(:).wfmdecays);

featsori  = [trepol t2ptime axonval gsigma unitamps wvfdecays bslfrate spikecv spikelv t2pval];
medianmat = repmat(median(featsori, 'omitmissing'), [size(featsori, 1), 1]);
featsori(isnan(featsori)) = medianmat(isnan(featsori));
featsquant  = (featsori - min(featsori, [], 1))./quantile(featsori,0.99,1);
featsquant(featsquant>1) = 1;
featsquant(featsquant<0) = 0;
%==========================================================================
% define distances and compute final similarity metric
Srob       = robustcov(featsori,'Method','olivehawkins');
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
    % 
    % corrmat = nan(size(feats1.templates, 4), size(feats2.templates, 4));
    % for i1 = 1:size(feats1.templates, 4)
    %     temp1 = denoiseTemplate(squeeze(feats1.templates(finalmatch(ii,1),:,:, i1)), 3);
    %     for i2 = 1:size(feats2.templates, 4)
    %         temp2 = denoiseTemplate(squeeze(feats2.templates(finalmatch(ii,2),:,:, i2)), 3);
    %         corrmat(i1, i2) = maxSlidingCorr(temp1,  temp2, 40);
    %     end
    % end
    % fcorr(ii) = max(corrmat, [], 'all');

    temp1 = squeeze(feats1.templates(finalmatch(ii,1),:,:));
    temp2 = squeeze(feats2.templates(finalmatch(ii,2),:,:));
    fcorr(ii) = maxSlidingCorr( denoiseTemplate(temp1, 3),  denoiseTemplate(temp2, 3), 40);

    % temp1 = feats1.avgwforms(finalmatch(ii,1),:);
    % temp2 = feats2.avgwforms(finalmatch(ii,2),:);
    % fcorr(ii) = maxSlidingCorr(temp1, temp2, 40);

    dfall(ii) = finalsim(finalmatch(ii,1), finalmatch(ii,2));
end
fprintf('Done! Took %2.2f\n', toc);
%%
%==========================================================================
end