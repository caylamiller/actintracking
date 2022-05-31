function [lsqparams, mleparams, mlecis] = FitWbl(results, ...
        noiseparams, fitts, paramseed, datamax,  ed, x1, xyres, spaceunits, ...
        frSep, timeunits, pop, map)
    
% initialize fit parameters:
Ntimes = max(fitts);
mleparams=zeros(Ntimes,2);
lsqparams=zeros(Ntimes,2);
mlecis=zeros(Ntimes,4);

figure;
Nmap = length(map);

% set range over which to integrate displacements in numerical 
% approximation of Bayes theorem, in spaceunits
srange  = linspace(0.0001,1000,100); 

% specify noise parameters:
f1 = noiseparams(1);
sig1 = noiseparams(2);
sig2 = noiseparams(3);

% initialize legend text:
legendtext  = {} ;
counter = 1; 

%%
for k=flip(fitts) % iterate over fit times, in reverse order
    % timescale, in timeunits:
    t=k*frSep;
    
    % initialize best fit parameters to zeros:
    phat=[0 0 ];
    mleci=[0 0 0 0];
    
    fprintf('fitting time %d', k)
    
    % find all displacements for the wanted timescale:
    steplist = [results.Nframes]==k ;
    
%% first use LSQ fitting to estimate parameters:
    myfunLSQ = @(p,xdata) WblWithNoise(xdata, sig1^2, ...
        sig2^2, f1, srange, p(1), p(2));
%    options = optimoptions(@lsqcurvefit,'FunctionTolerance',1e-9);
    
    % calculate the pdf for these displacements:
    pdf = histcounts([results(steplist).step]*xyres, ed, ...
        'Normalization','pdf');
    
    lsqparam = lsqcurvefit(myfunLSQ,paramseed,x1,pdf,[0 0],[]);
    
    lsqparams(k,:) = lsqparam;
    
%% Now do MLE fitting:
    % if there are more displacements than datamax, take a subset:
    if sum(steplist) > datamax % check if there are too many displacements
        % pull a random sample of displacements:
        steplist = datasample(find(steplist),datamax);
    end
    
    % function needs to be defined differently for the mle solver:
    myfunMLE = @(x,shape,scale) WblWithNoise(x, sig1^2, sig2^2,...
        f1, srange, shape, scale);

    try % try/catch loop ensures that loop continues onto next fit, even if one fit fails
        % do the mle solving:
        [phat,mleci] = mle([results(steplist).step]*xyres,...
            'pdf',myfunMLE,'start',lsqparams(k,:),'lowerbound',[0 0]);
    catch
    end
    
    mleparams(k,:) = phat(:);
    mlecis(k,:) = mleci(1:4);
    
    %%
    color = map(round((k-1)*(Nmap-1)/(Ntimes-1)+1),:);
    % plot measured pdf:
    plot(x1,pdf,'color',color,'linewidth',2);
    hold on;
    % and the corresponding best fit line
    bestfit= WblWithNoise(x1, sig1^2, sig2^2, f1, srange,...
        phat(1), phat(2));
    plot(x1,bestfit,'color',color,'linewidth',1);
    
    % update legend text:
    legendtext{counter} = [num2str(t) timeunits ' data'];
    legendtext{counter+1} = [num2str(t) timeunits ' fit'];
    counter = counter + 2; 
end

% Plot formatting:
title([pop ' displacement fits'])
xlabel(['displacement (' spaceunits ')'])
ylabel('pdf')
box off
xlim([ed(1) ed(end)])



