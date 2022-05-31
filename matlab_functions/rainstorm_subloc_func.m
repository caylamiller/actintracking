function results = rainstorm_subloc_func(results, actStack)
totFr = length(actStack); 
%% Initialize localisation and thresholding variables
initX0  = 0; 
rad     = 3;    % Radius for fitting. Square LN(-rad:rad)
tol     = 0.1;  % Tolerance for least squares fit 
maxIts  = 6;    % Number of fitting iterations 
lam_px  = 674/64 ; % emission wavelength 674 nm / 64 px
NA      = 1.49 ; 
initSig = 0.25*lam_px/NA;  % Initial guess of P.S.F. sigma
allowSig= [initSig/1.5 initSig*1.5]; % Reject fits with way-out sigma values (0.8--3)
allowX  = 2;  % Reject localisations if abs(x0) > over  (2)

clearvars lam_px NA
%% 
Ntracks = length(results);  % total # of tracks
lenlist = [results.length] ;%  the length of each track

maxL = max(lenlist) ;       % length of the longest track
allx = zeros(Ntracks, maxL) ;% to be filled. x-values for all tracks
ally = allx ;               % to be filled. y-values for all tracks

initfr  = [results.initfr] ;    % frame each track first appears in
finfr   = initfr + lenlist - 1; % frame each track appears last in
%% get x and y lists
for tracknum = 1:Ntracks % for each track
    currx = results(tracknum).x ; % x values of current track
    curry = results(tracknum).y ; % y values of current track
    L = lenlist(tracknum) ;       % length of current track
    allx( tracknum, 1:L) = currx ;  % fill matrix
    ally( tracknum, 1:L) = curry ;  % fill matrix
end
allxsize = size(allx); 
newx = NaN(allxsize); % to be re-written by subpix localization later
newy = NaN(allxsize); % to be re-written by subpix localization later
sigx = NaN(allxsize);
sigy = NaN(allxsize);
%%
for currfr = 1:totFr % loop through frames (not tracks!)
img = actStack(currfr).data;
dims = size(img); 
    
    subset =  [initfr' <= currfr] .* [finfr' >= currfr]; % find the tracks
                                         % which are present in this frame
    trN = find(subset); % indices of the subset tracks
    fr = currfr + 1 - initfr(trN)' ; % get the frame # (within a track) this 
                                 % frame corresponds to
    
    % pull x & y coords from tracks for current frame:
    idxlist1 = sub2ind(allxsize, trN, fr); 
    xlist = allx(idxlist1);
    ylist = ally(idxlist1);

%    not_edge = logical([xlist > rad] .* [ylist > rad] .* ...
%        [xlist < dims(2)-rad] .* [ylist < dims(1)-rad]);
%    kxy = find(not_edge); 
    
    idxlist2 = sub2ind(dims, ylist, xlist); 
    c3 = img(idxlist2) ;
    bkgdSig = std(double(img(img < mean(img(:))))); % Avoids signal
    
    %%
    [positions,params] = rainSTORM_fitLocGF(img,...
        [xlist ylist c3],...
        initX0,initSig,allowSig,rad,tol,allowX,bkgdSig,maxIts);
    
    fitted_indx = positions(:,1)~=-1; 
   
    %%
    newx(idxlist1((fitted_indx))) = positions(fitted_indx,1);
    newy(idxlist1((fitted_indx))) = positions(fitted_indx,2);
    sigx(idxlist1((fitted_indx))) = params(fitted_indx,4);
    sigy(idxlist1((fitted_indx))) = params(fitted_indx,5);
%     sigx(idxlist1((fitted_indx))) = params(fitted_indx,7);
%     sigy(idxlist1((fitted_indx))) = params(fitted_indx,8);
%% plot, for troubleshooting
%     figure;
%     imagesc(img)
%     axis image
%     colormap gray
%     hold on;
%     plot(positions(fitted_indx,1), positions(fitted_indx,2), 'ob')
%     plot(allx(idxlist1(~fitted_indx)),ally(idxlist1(~fitted_indx)), 'xr')
%     title(num2str(currfr))
%     
%     h=gcf;
%     h.Position = [1 41 1680 933];

end
%%

for tracknum = 1:Ntracks
    
    L = lenlist(tracknum) ;
    results(tracknum).xsub = newx( tracknum, 1:L) ;
    results(tracknum).ysub = newy( tracknum, 1:L) ;
    results(tracknum).sigx = sigx( tracknum, 1:L) ;
    results(tracknum).sigy = sigy( tracknum, 1:L) ;
end
%%

end
