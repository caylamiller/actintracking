function f = AnalyticalJumpModel(x,b,tmu)
f = zeros(size(x));
x0 = [x==0];
f(x0) = exp(-tmu) ; 
    
    f1 = exp(-tmu-b.*x(~x0)); 
    f2 = (tmu.*b./x(~x0)).^0.5;
    f3 = besseli(1, 2*(tmu.*b.*x(~x0)).^0.5); 
    f(~x0) = f1.*f2.*f3; 