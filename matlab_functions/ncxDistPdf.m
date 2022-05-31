function f = ncxDistPdf(d, sig2, s)
% computes the pdf of the non-centrally chi-distributed distances between 2
% gaussian distributed points. 
% This does not give a true non-central chi distribution, but is scaled by
% (sqrt(2)*sigma)^-1
% Parameters:
% f = f(d), pdf at some distance
% d = distance
% sig2 = sigma^2 for the gaussian distributed points
% s = true underlying distance (distance between the centers of the 2
% guassians)

twosigsq = 2 * sig2; 

inbessel = d.*s/twosigsq;
inbessel700 = inbessel>700; % matlab fails to solve the bessel for arguments > 700

f = (d/twosigsq) .* exp(-(0.5/twosigsq)*(d.*d+s.*s)) .* besseli(0, inbessel);

f2 = normpdf(s, d, sqrt(twosigsq)); % at high values, we can approximate by a gaussian 

f(inbessel700) = f2(inbessel700);
