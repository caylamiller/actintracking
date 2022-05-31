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
paramseed = [50, 10];   % starting parameters to seed analytical jump model 
                        % fit; [mean jump distance, mean time btw jumps],
                        % in [spaceunits, timeunits]
datamax = 5000;         % maximum number of points to fit (to save time, 
                        % set lower for preliminary fits)
fitts = [10 20 30 40];  % timescales for which the data will be fit, in 
                        % spaceunits; must be a multiple of frSep
fitfrs = fitts/frSep;  % Convert from timeunits to frames
Ntimes = max(fitfrs);   % total number of timescales to calculate
vrange=0.001:0.5:60;    % velocities over which to plot pdfs, in 
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

%% Check if output directory exists, create if not %%%%%%%%%%%%%%%%%%%%%%%%
if ~exist(outdir, 'dir')
    mkdir(outdir)
end
%% Take QFSM files, pull out tracks, do subpixel localization, and sort 
% tracks by subcellular location:
% Use this function to generate allresults structure from 
% allresults = generateTrackFiles(imgnums, datadir, trackdirsuff, xyres, ...
%     spaceunits, frSep, timeunits, n2vdir, driftTabledir, maskdir, ...
%     outdir, trackOutfileSuffix);

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
    [lsqparams, mleparams, mlecis] = Fit1DAnalyticalJumpModel(results, ...
        noiseparams, fitfrs, paramseed, datamax,  ed, x1, xyres, spaceunits, ...
        frSep, timeunits, pop, map);
    
    saveas(gcf, fullfile(outdir,[date '_' pop '_fit_displacements.fig']));
    eval([pop 'paramsMLE = mleparams;']);
    eval([pop 'cisMLE = mlecis;']);
    eval([pop 'paramsLSQ = lsqparams;']);
    
    %% Plot the output velocity distributions:
    figure;
    for k=fitfrs
        timesep = k*2;
        % convert vrange from velocities to displacments for the given
        % timescale (because the jump model is formulated in terms of
        % displacements):
        drange = vrange*timesep;
        
        % calculate probability of zero velocity:
        y0 = AnalyticalJumpModel(0,1/mleparams(k,1),timesep/mleparams(k,2)) ;
        
        % calculate pdfs for drange:
        yvals = AnalyticalJumpModel(drange, 1/mleparams(k,1), ...
            timesep/mleparams(k,2));
        % scale by the time to convert back to velocities:
        yvals = yvals*timesep;
        
        color = map(round((k-1)*(Nmap-1)/(Ntimes-1)+1),:);
        % plot discrete probability of zero displacement/velocity:
        plot(0,y0, '.', 'color', color);
        hold on;
        % plot pdf for remainder of vrange:
        plot(vrange, yvals, 'color', color, 'linewidth', 2)
        
    end
    xlabel(['velocity (' spaceunits '/' timeunits ')'])
    ylabel('pdf')
    title([pop ' best fit velocity distributions'])
    box off;
    
saveas(gcf, fullfile(outdir,[date '_' pop '_velocity_distributions.fig']));

end

% Plot the results of all fits:
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
    ylabel(['mean jump distance (' spaceunits ')']);
    box off; 
    
    subplot(1,2,2);
    hold on;
    % plot the mean time between jumps parameter:
    errorbar(fitts, params(fitfrs,2), [cis(fitfrs,4)-cis(fitfrs,3)]/2);
    xlabel(['timescale (' timeunits ')']);
    ylabel(['mean time btw jumps (' timeunits ')']);
    box off;
end

legend(pops)
saveas(gcf, fullfile(outdir,[date '_analytical_jump_model_fit_params.fig']));

save(fullfile(outdir, [date '_analytical_jump_model_fit_params.mat']),...
    '*paramsMLE','*cisMLE','fitts')