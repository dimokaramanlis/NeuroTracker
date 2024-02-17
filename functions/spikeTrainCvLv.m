function [cv, lv] = spikeTrainCvLv(sptrains, fs)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here


spks    = double(CelltoMatUE2(sptrains))/fs;
allisis = diff(spks, [], 2);
cv = std(allisis, [], 2, 'omitnan')./mean(allisis,2,'omitnan');

lv = 3* sum(((allisis(:, 2:end) - allisis(:, 1:end-1)).^2)./...
    ((allisis(:, 2:end) + allisis(:, 1:end-1)).^2),2,'omitnan');
normfac = sum(~isnan(allisis(:, 2:end-1)), 2) - 1;
lv = lv ./normfac;



end