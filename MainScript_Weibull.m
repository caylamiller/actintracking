%% SET UP VARIABLES: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clearvars; 
% Directories and file naming: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
outdir  = 'matlab_output' ;    % directory for matlab output files
datadir = 'sample_data' ;       % directory in which QFSM data is saved:
trackdirsuff = '_filt_output' ; % QFSM output directories in format 
                                % "#[trackdirsuff]"
maskdir = ['sample_data' filesep 'masks']; % directory in which mask images
                            % are saved; images in this dir are in the 
                            % format "#_cellmask.tif" and "#_pxnmask.tif"
n2vdir  = ['sample_data' filesep 'n2v'] ;   % directory in which actin 
                            % stacks (after noise2void) are saved; stacks
                            % in this dir are in the format "#.tif"
driftTabledir=  ['sample_data' filesep 'drift_tables']; % directory in 
                            % which drift tables are saved; drift tables in 
                            % this dir are in the format "#.csv"
trackOutfileSuffix = '_tracks_n2v' ; % suffix to be used in track file name

% image numbers, specified in image names and subdirectories:
imgnums  = [14 16 18] ;

% Imaging parameters: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
frSep       = 2 ; %time between frames, in s
timeunits   = 's' ;
xyres       = 64; % pixel size, in nm
spaceunits  = 'nm' ;

% Histogram parameters for displacement distributions: %%%%%%%%%%%%%%%%%%%%
edmin   = 0 ;                           % minimum step size (spaceunits)
edmax   = 900 ;                         % maximum step size (spaceunits)
w       = 25 ;                          % bin width (spaceunits)
ed      = edmin : w : edmax ;           % bin edges (spaceunits)
x1      = edmin+w/2 : w : edmax-w/2 ;   % bin centers (spaceunits)
Nbins   = length(x1);                   % number of bins

% Fitting parameters: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
noiseparams = [0.5717 27.0837 53.8431]; % noise parameters from fixed cell
                                        % fits, i.e. [f1, sig1, sig2]
paramseed = [1 50];   % starting parameters to seed analytical jump model 
                        % fit; [mean jump distance, mean time btw jumps],
                        % in [spaceunits, timeunits]
datamax = 5000;         % maximum number of points to fit (to save time, 
                        % set lower for preliminary fits)
fitts = [10 20 30 40];  % timescales for which the data will be fit, in 
                        % spaceunits; must be a multiple of frSep
fitfrs = fitts/frSep;  % Convert from timeunits to frames
Ntimes = max(fitfrs);   % total number of timescales to calculate
vrange=0:0.5:60;    % velocities over which to plot pdfs, in 
                        % spaceunits/timeunits
                        
% Set up colormap:  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load('batlow.mat'); % load colormap
map = flip(batlow); % flip direction
Nmap = size(map,1); % colormap length

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not need to change variables or functions below this line (for sample 
% data -- may need to modify for other data if files are stored/formatted 
% differently)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Take QFSM files, pull out tracks, do subpixel localization, and sort 
% tracks by subcellular location:
allresults = generateTrackFiles(imgnums, datadir, trackdirsuff, xyres, ...
    spaceunits, frSep, timeunits, n2vdir, driftTabledir, maskdir, ...
    outdir, trackOutfileSuffix);

%% Take tracks and calculate displacement on varying timescales:
stepresults = generateDisplacementsVaryingTimescales(allresults, Ntimes); 
% stepresults has 1 entry per measured displacement, rather than per track
% (it's generally ~an order of magnitude larger than allresults)
save(fullfile(outdir, [date '_all_displacements.mat']))

% population list, matching pop numbering introduced in the previous function:
pops={'SFs_adh','SFs','adhesions','cortical'}; 

%% Plot distributions of tracks
figure;
legendtext = {};
counter = 1;
for i = fitfrs
    steplist = [stepresults.Nframes]==i ;
    
    pdf = histcounts([stepresults([steplist]).step]*xyres, ed, 'Normalization','pdf');
    color = map(round((i-1)*(Nmap-1)/(Ntimes-1)+1),:);
    plot(x1,pdf,'color',color,'linewidth',2);
    hold on;
    legendtext{counter} = [num2str(i*frSep) timeunits];
    counter = counter + 1;
end

box off;
xlabel(['displacement (' spaceunits ')'])
ylabel('pdf')
legend(legendtext)
title('all populations')

saveas(gcf, fullfile(outdir,[date '_allpops_measured_displacements.fig']));

clearvars pdf color legendtext counter

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%% FIT TO ANALYTICAL JUMP MODEL %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for i=1:length(pops) % Iterate over the four population groups:
    
    % pull population name from list:
    pop = pops{i};
    fprintf(pop)
    
    % pull results that correspond to that population:
    results = stepresults([stepresults.pop]==i);
    
    % fit selected results to a model:
    [lsqparams, mleparams, mlecis] = FitWbl(results, ...
        noiseparams, fitfrs, paramseed, datamax,  ed, x1, xyres, spaceunits, ...
        frSep, timeunits, pop, map);
    
saveas(gcf, fullfile(outdir,[date '_' pop '_Wbl_fit_displacements.fig']));
    eval([pop 'paramsMLE = mleparams;']);
    eval([pop 'cisMLE = mlecis;']);
    eval([pop 'paramsLSQ = lsqparams;']);

%% Plot the output velocity distributions:
    figure;
    for k=fitfrs
        timesep = k*2;
        
        % calculate pdfs for vrange, dividing scale parameter by timescale
        % to plot velocity:
        yvals = wblpdf(vrange, mleparams(k,2)/timesep, ...
            mleparams(k,1));

        color = map(round((k-1)*(Nmap-1)/(Ntimes-1)+1),:);
        hold on;
        % plot pdf:
        plot(vrange, yvals, 'color', color, 'linewidth', 2)
        
    end
    xlabel(['velocity (' spaceunits '/' timeunits ')'])
    ylabel('pdf')
    title([pop ' best fit velocity distributions'])
    box off;
    
saveas(gcf, fullfile(outdir,[date '_' pop '_Wbl_velocity_distributions.fig']));

end

%% Plot the results of all fits:
figure;

for i=[1:4] % iterate over populations
    % pull parameters for this population:
    pop = pops{i};
    eval(['cis = ' pop 'cisMLE;']);
    eval(['params = ' pop 'paramsMLE;']);
    fprintf('plotting parameters for population: %s\n', pop)
    
    subplot(1,2,1);
    hold on;
    % plot the mean jump distance parameter:
    errorbar(fitts, params(fitfrs,1), [cis(fitfrs,2)-cis(fitfrs,1)]/2);
    xlabel(['timescale (' timeunits ')']);
    ylabel(['shape']);
    box off; 
    
    subplot(1,2,2);
    hold on;
    % plot the mean time between jumps parameter:
    errorbar(fitts, params(fitfrs,2)'./fitts, ...
        ((cis(fitfrs,4)-cis(fitfrs,3))'/2)./fitts);
    xlabel(['timescale (' timeunits ')']);
    ylabel(['scale (' spaceunits '/' timeunits ')']);
    box off;
end

legend(pops)
saveas(gcf, fullfile(outdir,[date '_Wbl_fit_params.fig']));

save(fullfile(outdir, [date '_Wbl_fit_params.mat']),...
    '*paramsMLE','*cisMLE','fitts')