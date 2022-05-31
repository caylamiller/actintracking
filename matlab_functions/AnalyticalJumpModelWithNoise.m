function fd = AnalyticalJumpModelWithNoise(dlist, sig1squared, sig2squared, ...
    fr1, srange, b, mu, t)
% Calculates the probability density of d (measured displacement) according
% to the jump model, with the correction for finite, discrete
% probability at true displacements of 0 ("0"), and with noise
% contributions from a mixture of two sigma values ("noisy 2") 
% Inputs:
% dlist = list of measured displacements (vector) for mle, or range of
% measured displacements, for lsq
% sig1squared, sig2squared = best fit sig1^2 and sig2^2 from fixed cell fit
% fr1 = fractional weighting of sig1 from fixed cell fit
% srange = range of true displacements over which to integrate (>0)
% b = 1/mean jump distance
% mu = 1/mean time between jumps
% t = observation timescale

L = length(dlist);
ds = mean(srange(2:end) - srange(1:end-1)) ; % width of s bins
fd = zeros(size(dlist));   % output, f(d)

s1 = srange(1:end-1);   % left-hand edges of s bins
s2 = srange(2:end);     % right-hand edges of s bins

for i=1:L
    d = dlist(i); 
    % calculate p(d|S) for each bin, averaging left and right bin edge:
    pd_Sj = 0.5*fr1*(ncxDistPdf(d, sig1squared, s1) + ncxDistPdf(d, sig1squared, s2)) + ...
        0.5*(1-fr1)*(ncxDistPdf(d, sig2squared, s1) + ncxDistPdf(d, sig2squared, s2)); 
    % calculate p(S) for each bin, averaging left and right bin edge:
    pSj = 0.5*(AnalyticalJumpModel(s1,b,t*mu) + AnalyticalJumpModel(s2,b,t*mu));
    % for each bin, calculate p(d|S)*p(S)*ds:
    f = ds .* pd_Sj .* pSj ;
    % sum up each bin and add in the discrete probability at s=0:
    fd(i) = sum(f) + (fr1*ncxDistPdf(d, sig1squared, 0) + ...
        (1-fr1)*ncxDistPdf(d, sig2squared, 0))*AnalyticalJumpModel(0,b,t*mu) ;
end
