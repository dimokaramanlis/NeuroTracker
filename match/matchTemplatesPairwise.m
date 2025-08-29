function [matchmat] = matchTemplatesPairwise(finalsim)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
%==========================================================================
Nunits = size(finalsim);
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
end