function dtemplate = denoiseTemplate(template, N)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

[aa,bb,cc] = svds(double(template), N);
dtemplate  = aa * bb * cc';
dtemplate  = single(dtemplate);
end