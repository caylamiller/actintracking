function allresults = generateTrackFiles(imgnums, datadir, dirsuff, xyres, ...
    spaceunits, frSep, timeunits, n2vdir, driftTabledir, maskdir, ...
    outdir, outfile_suffix)
%% This function generates matlab files with tracks stored in structure 
% named "results" for each image, as well as outputing a structure 
% "allresults"  that both pools together all tracks from multiple images 
% and specifies track location by subcellular grouping (stress fibers vs 
% cortical vs adhesions) 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INPUTS:
% imgnums = list of image numbers that correspond to image names (vector of
%           integers)
% datadir = string, directory in which QFSM data is saved (individual output
%           folders with QFSM resuls should be inside)
% dirsuff = string, specifies format of QFSM output folders 
%           i.e. folders are in format "[imgnum][dirsuff]"
% n2vdir  = directory in which noise2voided actin stacks are saved
% maskdir = string, directory in which paxillin mask images are saved 
% driftTabledir = directory in which drift tables are saved
% outdir  = string, directory in which to save matlab output files
% outfile_suffix = string, to be used in output files of tracks
%
% xyres   = pixel size, in spaceunits
% spaceunits  = string defining units of pixel size 
% frSep   = time between frames, in timeunits
% timeunits   = string defining units of frSep
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% OUTPUT:
% allresults = structure which holds all tracks, with the following
% fields:
%          imNum: which image number the track came from
%         length: the length (in frames) of the track
%        paxMask, 
%        actMask, 
%     paxactMask, 
%       cortMask: 0 or 1, specify the subcellular grouping of the track
%              x: x-coordinates of each point in the track (NaN = failed
%              subpixel localiation)
%              y: y-coordinates of each point in the track (NaN = failed
%              subpixel localiation)
%             fr: frame numbers during which the track occurred

%% Initiate variables for counter and results:
allresults={};  % will become structure holding track info for all images
trackcount=1;     % keeps track of number of tracks in allresults

%% Iterate over images, pull QFSM result tracks file, and store in allresults
for k=1:length(imgnums)
    tic
    imgnum = num2str(imgnums(k));
    fprintf('Collecting tracks for image number %s\n' , imgnum)
% define directory and file name of QFSM track file:    
     trackdir= fullfile( datadir,[ imgnum dirsuff ],'QFSMPackage',...
         'speckleTracks' ) ;
     trackfi = [ imgnum '_filt_tracks.mat' ] ;
     
% pull tracks in to structure "results"
     results = calctracks(fullfile(trackdir,trackfi), ...
        xyres, spaceunits, frSep, timeunits) ;
    
% read in stack of actin images, to be used now for subpixel localization
% and later for masking:
    fprintf('reading in actin stack...\n');
    actinstack = tiffread2(fullfile(n2vdir,[ imgnum '.tif']));
    fprintf('subpixel localization...\n');
    results = rainstorm_subloc_func(results, actinstack );
% save results in the output directory    
    outfile = fullfile(outdir,[date '_' imgnum outfile_suffix '.mat']) ;
    fprintf('Saving...\n') ;
    save( outfile, '*dir*','outfile_suffix','frSep','xyres',...
        'trackfi','results','actinstack') ;
    fprintf( 'Results saved as:\n%s\n' , outfile );
    toc

    %% keep only multi-point tracks (no single localizations):
    results = results([results.length]>1) ;
    
    %% load drift table:
    T = readtable(fullfile(driftTabledir,[ imgnum '.csv']));
    xdriftvect = T.X_Drift_pixels_; 
    ydriftvect = T.Y_Drift_pixels_;
    
    clearvars T 
    %% PAXILLIN MASK
    adhesmaskfi = [imgnum '_pxnmask.tif' ];
    adhesmask  = imread( fullfile( maskdir , adhesmaskfi ) )>0 ;
    
    clearvars adhesmaskfi
    %% REFORMAT ACTIN STACK AS A 3D MATRIX %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    newstack= zeros(actinstack(1).height,actinstack(1).width,length(actinstack));
    for i=1:length(actinstack)
         newstack(:,:,i)=actinstack(i).data;
    end
    tprojactin=sum(newstack,3); % take the time projection of the actin stack
    
    clearvars actinstack newstack
    %% CELL MASK
    cellmaskfi = [imgnum '_cellmask.tif' ];
    incell  = imread( fullfile( maskdir , cellmaskfi ) )>0 ;
    
    clearvars cellmaskfi
    %% STRESS FIBER MASK
    actinvals = sort(tprojactin(:));
    cutoff=actinvals(round(0.98*length(actinvals))) ;
    
    actmask = (tprojactin>=cutoff) .* incell ;
    clearvars tprojactin cutoff actinvals
    %% CORTEX MASK
    % masks out values outside stress fibers and adhesions, but in cell:
    cortmask = (actmask==0) .* incell .* (adhesmask==0);
    
    %% COLLECT TRACKS
    dims = size(adhesmask);

    for i=1:length(results) % iterate over tracks
        
        ff  = results(i).initfr; % first frame track appears in
        L   = results(i).length; % number of frames track appears in
        
        % pull drift correction values for relevant frames:
        xdrift = xdriftvect(ff:(ff+L-1)); 
        ydrift = ydriftvect(ff:(ff+L-1));
        
        % apply drift correction to the x and y coordinates of the track:
        y   = [results(i).ysub] + ydrift';
        x   = [results(i).xsub] + xdrift';
        
        % calculate 
        xx = max(1, min(dims(2), round(x(1))));
        yy = max(1, min(dims(1), round(y(1))));
        
        paxcount = adhesmask(yy,xx) > 0 ; % is 1 if track originates over pxn adhesion
        actcount = actmask(yy,xx) > 0 ; % is 1 if track originates over stress fiber
        paxactcount = paxcount * actcount ; % is 1 if track originates over pxn adhesion AND over stress fiber
        cortcount =  cortmask(yy,xx) > 0 ; % is 1 if track originates in cortex
        cellcount = incell(yy,xx) > 0 ; % is 1 if track originates inside cell
        
        if cellcount % only keep tracks from inside cell mask
            allresults(trackcount).imNum = imgnum;
            allresults(trackcount).length = L;
            allresults(trackcount).paxMask = paxcount - paxactcount ;
            allresults(trackcount).actMask = actcount - paxactcount ;
            allresults(trackcount).paxactMask = paxactcount ;
            allresults(trackcount).cortMask = cortcount ;
            allresults(trackcount).x = x ;
            allresults(trackcount).y = y ;
            allresults(trackcount).fr = ff:(ff+L-1) ;
            trackcount = trackcount+1;
        end
    end
    
    clearvars x y adhesmask actmask cortmask incell xx yy ...
        paxcount actcount paxactcount cytcount edgcount cellcount ...
        ff l xdrift ydrift xdriftvect ydriftvect results
    
    toc
end

save(fullfile(outdir, [date '_all_tracks' outfile_suffix '.mat']))