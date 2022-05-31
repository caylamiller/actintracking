function [results] = calctracks(tracksfi, xyres, spaceunits, ...
    timeres, timeunits)
% This function pulls results from QFSM-created track file and stores them
% in an output structure for downstream analysis.
% INPUTS:
% tracksfi  = path to QFSM track file
% xyres     = image resolution in spaceunits/pixel
% spaceunits= units used to define xyres
% timeres   = time resolution in timeunits/frame
% timeunits = units used to define timeres
% OUTPUTS:
% results   = structure storing all track info:
% results(i).x = x position of all points in track i (pixel
% localization)
% results(i).y = y position of all points in track i (pixel
% localization)
% results(i).initfr = frame # in which track i appears
% results(i).length = number of frames for which track i is tracked
% results(i).v_unit = velocity units, defined as spaceunits/timeunits
% results(i).mean_v = mean frame-to-frame speed of track i

results     = struct([]); % initiate results structure
trackiter   = 1; % initiate track counter

load(tracksfi) ; % Load QFSM output file, giving tracks matrix MPM:
dims        = size(MPM);

fprintf('Pulling results and calculating frame-to-frame velocities.\n');

for row=1:dims(1)
    % initiate vectors to store x & y location and keep track of frame
    % number:
    x=[]; y=[]; fr=0;
    
    for col=1:2:(dims(2)-1) % iterate over every other column in MPM
        % (column format is ypos1, xpos1, ypos2, xpos2, ...
        
        if (MPM(row,col)+MPM(row,col+1))>0; % check for x & y > 0; x,y==0
            % is used in MPM to separate tracks
            y=[y,MPM(row,col)];  % add to y vector
            x=[x,MPM(row,col+1)];% add to x vector
            fr=(col+1)/2;        % keep track of frame number (# columns =
            % # frames * 2)
        else % if x,y==0, we've reached the end of a track and will begin a
            % new one on the next iteration
            tracklen=length(x); % define track length by the number of pts
            
            if tracklen > 1 % identify whether a true track (2+ pts) or an
                % isolated localization
                
                % calculate total frame-to-frame distances of all points:
                totdist     = (x(1:(end-1))-x(2:end)).^2 + ...
                    (y(1:(end-1))-y(2:end)).^2 ;
                totdist = sum(totdist.^0.5) ;
                % calculate average frame-to-frame velocity:
                avgv = totdist * xyres / ( timeres * (tracklen-1) ) ;
                if avgv < 0
                    pause % this should only catch if there is an error in
                    % input files, was used for troubleshooting but should
                    % be useless now...
                end
                % save all info in results structure:
                results(trackiter).x = x;
                results(trackiter).y = y;
                results(trackiter).initfr = fr+1-tracklen;
                results(trackiter).length = tracklen;
                results(trackiter).v_unit = [spaceunits, '/', timeunits];
                results(trackiter).mean_v = avgv; % mean velocity
                trackiter = trackiter + 1;
            end
            x=[]; y=[]; % re-initiate x & y
        end
    end
end
end
