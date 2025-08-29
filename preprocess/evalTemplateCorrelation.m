function [fcorr, dfall] = evalTemplateCorrelation(finalmatch, finalsim, templates1, templates2)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
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

    temp1 = squeeze(templates1(finalmatch(ii,1),:,:));
    temp2 = squeeze(templates2(finalmatch(ii,2),:,:));
    fcorr(ii) = maxSlidingCorr( denoiseTemplate(temp1, 3),  denoiseTemplate(temp2, 3), 40);

    % temp1 = feats1.avgwforms(finalmatch(ii,1),:);
    % temp2 = feats2.avgwforms(finalmatch(ii,2),:);
    % fcorr(ii) = maxSlidingCorr(temp1, temp2, 40);

    dfall(ii) = finalsim(finalmatch(ii,1), finalmatch(ii,2));
end
fprintf('Done! Took %2.2f\n', toc);

end