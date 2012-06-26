function [V,D,Dperc,PCs]=hilberteof(data,neof)
% HILBERTEOF Function to find the Hilbert (or complex) EOF of a
%   data set.
%   [V,D,Dperc,PCs]= HILBERTEOF(DATA,NEOF) finds the 
%   first NEOF Hilbert EOFS of the dataset DATA. The dataset 
%   should have dimensions NT*NX, where the first dimension 
%   represents time, and the second dimension encasulates all 
%   the spatial dimensions. 
%
%   The function returns the EOFs (V) the variance explained (D)
%   the variance explained as a percentage (Dperc) and the
%   principal components (PCs).
%
%  v1.1 AH 6/9/2004

%   VERSION LOG
%   v1.0 - 29/10/2003 - AMH1, JDHA, DDC. Based on code by RYH, PC & GDQ
%   v1.1 - took nt & nx out of arg list - AH 6/9/04
  
  disp('Now finding Hilbert transform ...');
  ynew=hilbert(data);
  data=ynew;
  ynew=[];clear ynew
    
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
  