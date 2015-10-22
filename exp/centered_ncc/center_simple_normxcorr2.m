function [xcr] = center_simple_normxcorr2(samp, intr)
% function [xcr] = center_simple_normxcorr2(samp, intr)
%
% Center each data window once, then compute the normalized cross correlaton.
% %

mask_samp = samp ~= 0;
mean_samp = mean(samp(mask_samp));
samp(mask_samp) = samp(mask_samp)-mean_samp;

mask_intr = intr ~= 0;
mean_intr = mean(intr(mask_intr));
intr(mask_intr) = intr(mask_intr)-mean_intr;

xcr = normxcorr2(samp, intr);

