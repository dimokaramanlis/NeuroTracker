function finalsim = calculateUnitSimilarity(feats1, feats2)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

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
%==========================================================================
% define distances and compute final similarity metric
Srob       = robustcov(featsori,'Method','olivehawkins');
distallrob = pdist2(featsori(1:Nunits(1), :), featsori(Nunits(1)+1:end, :), 'mahalanobis', Srob);

Smult = Sdist;
Smult(Smult ==0) = 1e-3;
finalsim = Smult./distallrob;
%==========================================================================
end