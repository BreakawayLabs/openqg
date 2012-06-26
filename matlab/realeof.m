function [V,D,Dperc,PCs]=realeof(data,neof)
% REALEOF Function to find the EOF of a data set.
%   [V,D,DPERC,PCS]= REALEOF(DATA,NEOF) finds the 
%   first NEOF EOFs of the dataset DATA. The dataset 
%   should have dimensions NT*NX, where the first dimension 
%   represents time, and the second dimension encasulates all 
%   the spatial dimensions. 
%   The function returns the EOFs (V) the variance explained (D)
%   the variance explained as a percentage (DPERC) and the
%   principal components (PCS).
%
%  v1.1 AH 6/9/2004

%   VERSION LOG
%   v1.0 - taken from HILBERTEOF by AH JDHA, DDC, which was 
%         based on code by RYH, PC and GDQ  -- AH 19/3/2004 
%   v1.1 - took nt & nx out of arg list - AH 6/9/04


  % Find the total variance of the dataset. We will use this below
  % in calculating the percentage of variance explained by each mode. 
  SD = sum(var(data));
      
  disp('Now finding covariance matrix ...');
  XX = cov(data);
        
  disp('Now finding eigenvalues ...');
  [V,Dtmp]=eigs(XX,neof);
  XX = [];clear XX
  
  % V is NX x NX array containing eigenfunctions
  
  % Dtmp is an NX x NX diagonal array containing the eigenvalues,
  % in ascending order of significance. We now write a column
  % vector containing the absolute values of the eigenvalues, and
  % call it D:
  disp('Now we calculate the relative variance of each mode ...');
  for ii=1:neof
    D(ii,1)=abs(Dtmp(ii,ii));
  end
  Dtmp=[];clear Dtmp
  Dperc=D(:)/SD*100;
    
  disp('Now we calculate the principal components of each mode ...');
  PCs = data*V;
  