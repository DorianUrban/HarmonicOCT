function [corr] = tfxcorr(spec, offset)
%TFXCORR Summary of this function goes here
%   Detailed explanation goes here

[~, I] = max(max(abs(spec), [],  1));

% ref = floor(size(spec, 2) / 2) + offset;
ref = I + offset;
template = abs(spec(:, ref));
corr = normxcorr2(template, abs(spec));

end

