% Gausian Fit. Local Background Subtraction. 
function [myFits,myParams] = rainSTORM_fitLocGF(myFrame,myPixels,initX0,initSig,allowSig,rad,tol,allowX,bkgdSig,maxIts); 

%% Work on rows then cols. Reject fits with far-out x0, sigX, or residual.
Npts    = size(myPixels,1);
myFits  = -ones(Npts,2);    % [row col] matrix for each candidate
myParams = -ones(Npts,7);   % Parameters for accepted fits
Nfails  = 0;        % Reset the number of rejected fits for this new frame
%%
for lpPx = 1:Npts % For local maxima in descending order
myCol = myPixels(lpPx,1);   % col = x
myRow = myPixels(lpPx,2);   % row = y
dims = size(myFrame);       % image dimensions

%% Added by CMM, 23 Jan 2018; accounts for ROIs near image edge. 
myROI = myFrame(max(1,myRow-rad):min(dims(1),myRow+rad), ...
    max(1,myCol-rad):min(dims(2),myCol+rad)) ;

mean_intensity = mean(myROI(:));
myROI = myROI - min(myROI(:));  % square region to fit. Subtract minimum.

offcent = (2*rad+1) - size(myROI); % check whether pt is off center 
% (because near an edge of image). offcent(1)=row=y, offcent(2)=col=x

if myRow <= rad % top edge
    offcent(1) = -offcent(1) ; 
end
if myCol <= rad % left
    offcent(2) = -offcent(2) ; 
end
%% 
flagRowFits = false;   % Begin by noting the centre-position is not fitted
flagColFits = false;

if offcent(1) < 0 % ROI abuts top edge
yy = ((-rad-offcent(1)):rad)'; % y-positions (rows) (pixel widths) as column vector
else
yy = (-rad:(rad-offcent(1)))'; % y-positions (rows) (pixel widths) as column vector
end
if offcent(2) < 0 % ROI abuts left edge
xx = ((-rad-offcent(2)):rad)'; % x-positions (cols) (pixel widths) as column vector
else
xx = (-rad:(rad-offcent(2)))'; % x-positions (cols) (pixel widths) as column vector   
end
yRows = sum(myROI,2);  % Sum together columns (dim=2) to get row intensities (y)
xCols = sum(myROI,1)'; % Sum together rows (dim=1) to get col intensities (x)

%% Fit a Gaussian to the observed intensities, binned by row: (x)
x0 = initX0;
sigX = initSig;
Cx  = xCols(rad+1+offcent(1)); % Guess height of f(x). Centre value is a good guess.
  fofX = Cx*exp(-(xx-x0).^2/(2*sigX^2)); % Initial guess of f(x)
  Beta = xCols - fofX; % Change needed in f(x)
  for lpLSF = 1:maxIts
  A = [fofX/Cx,fofX.*(xx-x0)/sigX^2,fofX.*(xx-x0).^2/sigX^3]; % Jacobian
  b = A'*Beta;
  a = A'*A;
  dL= a\b;
  Cx = Cx+dL(1);
  x0 = x0 + dL(2);
  sigX = abs(sigX + dL(3));
  fofX = Cx*exp(-(xx-x0).^2/(2*sigX^2)); 
  Beta = xCols - fofX;
    if(abs(x0)>allowX || (sigX < allowSig(1)) || (sigX > allowSig(2)) )
%fprintf('%d 1 %f %f\n',lpPx,x0, sigX)
      break; % Stop iterating if solution drifts too far
    end
  end
  %%
  % Judge the fit. Accept if residue is a small proportion of |y^2|, etc.
  residueCols = sum(Beta.^2)/sum(xCols.^2);
  if (residueCols<tol && abs(x0)<allowX && sigX>allowSig(1) && sigX<allowSig(2))
  fitColPos = double(myCol)+x0;%-0.5; % Note (-0.5) for image registration
  flagColFits = true; % Flag the row-direction fit as acceptable
%    else
%  fprintf('%d res=%f\ty=%f\tsig=%f\n',lpPx,residueCols,x0,sigX)
 %figure; imagesc(myROI);
 %hold on; plot(initX0+5,x0+5,'x')
 %pause()
  end
 %%
  % Fit a Gaussian to the observed intensities, binned by Col
  if(flagColFits) % Don't fit row-direction if the col-axis fit was rejected
  y0 = initX0;
  sigY = sigX;% Keep sigX from col-fit. It should match the row-fit.
  Cy = yRows(rad+1+offcent(2));
  fofX = Cy*exp(-(yy-y0).^2/(2*sigY^2)); % Initial guess of f(x)
  Beta = yRows - fofX; % Change needed in f(x)
  for lpLSF = 1:maxIts
  A = [fofX/Cy,fofX.*(yy-y0)/sigY^2,fofX.*(yy-y0).^2/sigY^3]; % Jacobian
  b = A'*Beta;
  a = (A'*A);
  dL= a\b;
  Cy = Cy+dL(1);
  y0 = y0 + dL(2);
  sigY = abs(sigY + dL(3));
  fofX = Cy*exp(-(yy-y0).^2/(2*sigY^2)); 
  Beta = yRows - fofX;
   if(abs(y0)>allowX || (sigY < allowSig(1)) || (sigY > allowSig(2)) )
%fprintf('%d 2 %f %f\n',lpPx,y0,sigY)
     break; % Stop iterating if solution drifts too far
   end
  end
  %%
  % Judge the column direction fit:
  residueRows = sum(Beta.^2)/sum(yRows.^2);
   if (residueRows<tol && abs(y0)<allowX && sigY>allowSig(1) && sigY<allowSig(2))
    fitRowPos = double(myRow)+y0;
    flagRowFits = true; % Flag the column-direction fit as acceptable
%    else
%  fprintf('%d res=%f\ty=%f\tsig=%f\n',lpPx,residueRows,y0,sigY)
%figure; imagesc(myROI);
%hold on; plot(y0+5,x0+5,'x')
   end

    end % End if, which only tries fitting Col-wise if Row-wise fitted OK
 
 if(flagRowFits && flagColFits )     % Accept iff Row and Col fits are OK
  myFits(lpPx,:)=[fitColPos,fitRowPos];
  myParams(lpPx,1)=myPixels(lpPx,3); % Magnitude of this signal's center
  myParams(lpPx,2)=(residueRows+residueCols)/2; % Mean critical tol for fit
  myParams(lpPx,3)=mean_intensity; % Mean of the signal for this fit
  myParams(lpPx,4)=sigX;  % X-width (sigma, rows, fitted) of this Gaussian
  myParams(lpPx,5)=sigY;  % Y-width (sigma, cols, fitted) of this Gaussian
  myParams(lpPx,6)=bkgdSig;  % Background for each ROI
%   myParams(lpPx,7)=Cx;  % Background for each ROI
%   myParams(lpPx,8)=Cy;  % Background for each ROI
 else 
  Nfails = Nfails+1;
 end
 % No halting condition for multiple fitting failures are applied in this

end  % Loop to the next local maximum

%myParams = myParams(myFits(:,1)~=-1,:); %Params of accepted fits only
%myFits = myFits(myFits(:,1)~=-1,:); % Return accepted fits only. (-1)s are rejected fits.

end