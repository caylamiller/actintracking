function fd = WblWithNoise(drange, sig1squared, sig2squared, fr1, srange, ...
    wblShape, wblScale)
% Calculates the probability density of d (measured displacement) from an
% underlying weibull distribution of true displacements
% Inputs:
% drange = list of measured displacements (vector) for mle, or range of
% measured displacements, for lsq
% sig1squared, sig2squared = best fit sig1^2 and sig2^2 from fixed cell fit
% fr1 = fractional weighting of sig1 from fixed cell fit
% srange = range of true displacements over which to integrate 
% wblShape = shape parameter
% wblScale = scale parameter

L = length(drange);
ds = mean(srange(2:end) - srange(1:end-1)) ;
fd = zeros(size(drange)); 

s1 = srange(1:end-1);   % left-hand edges of s bins
s2 = srange(2:end);     % right-hand edges of s bins

for i=1:L
    d = drange(i);
    % calculate p(d|S) for each bin, averaging left and right bin edge:
    pd_Sj = 0.5*fr1*(ncxDistPdf(d, sig1squared, s1) + ncxDistPdf(d, sig1squared, s2)) + ...
        0.5*(1-fr1)*(ncxDistPdf(d, sig2squared, s1) + ncxDistPdf(d, sig2squared, s2)); 
    % calculate p(S) for each bin, averaging left and right bin edge:
    pSj = 0.5*(wblpdf(s1,wblScale,wblShape) + wblpdf(s2,wblScale,wblShape));
    pSj = pSj / sum(pSj)/ ds; 
    % for each bin, calculate p(d|S)*p(S)*ds:
    f = ds .* pd_Sj .* pSj ;
    % sum up each bin:
    fd(i) = sum(f) ;
end



