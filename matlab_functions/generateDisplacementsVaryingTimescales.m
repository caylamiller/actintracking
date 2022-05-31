function newresults = generateDisplacementsVaryingTimescales(allresults, Ntimes)

newresults={};
counter=1;

for i=1:length(allresults)
% iterate through tracks, identify population:

%    if allresults(i).cellMask>0     %FOR FIXED CELLS
%        pop=1;                      %FOR FIXED CELLS
   if allresults(i).paxactMask>0
       pop=1;
   elseif allresults(i).actMask>0  %ACT POP=1 USED FOR HUVECS (no pxn pop)
         pop=2;                    
   elseif allresults(i).paxMask>0
       pop=3;
   elseif allresults(i).cortMask>0
       pop=4;
    else
        continue
%        pop=2;                      %FOR HUVECS
    end
    L=allresults(i).length;
    imnum = (allresults(i).imNum);
    
    for k=1:min(L,Ntimes) % for each track, iterate over varying timescales
        % pull x and y coordinates from track:
        y=[allresults(i).y];
        x=[allresults(i).x];
        
        % calculate all displacements for this timescale:
        displacement = ((x((k+1):end) - x(1:end-k)).^2 + ...
            (y((k+1):end) - y(1:end-k)).^2 ).^0.5;
        
        % get rid of NaN displacements (due to missing pts/failed subpixel localizations:
        displacement = displacement(~isnan(displacement));
        
        if isempty(displacement)
            continue
        end
        
        for ii=1:length(displacement)
            newresults(counter).step = displacement(ii); % displacement, in px
            newresults(counter).Nframes = k; % number of frames over which displacement occured
            newresults(counter).pop = pop; % population (1-4)
            newresults(counter).imnum = imnum; 
            counter=counter+1;
        end
    end
end