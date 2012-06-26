function [U,V,D,PCs1,PCs2]=jointsvd(data1,data2,neof,nt)
% JOINTSVD Function to find the joint SVD of two data sets.
%   [U,V,D,PCS1,PCS2]= JOINTSVD(DATA,NEOF,NT,NX) finds the 
%   first NEOF joint SVDs of the datasets DATA1 & DATA2. The 
%   datasets should have dimensions NT*NX, where the first dimension 
%   represents time, and the second dimension encasulates all 
%   the spatial dimensions. They should be normalised and de-meaned. 
%   The function returns the spatial modes (U,V) the variance 
%   explained (D) and the principal components (PCS1,PCS2).
%
%  v1.0 AH 3/9/2004

%  v1.0 -  created using recipe from Bretherton et al. (1992) 
%       -  AH 3/9/2004

  disp('Now finding cross-covariance matrix ...');
  Csz = data1'*data2;
        
  disp('Now finding SVD ...');
  [U,S,V]=svds(Csz,neof);
  Csz = [];clear Csz
  
  % U,V are NX x NEOF arrays containing spatial modes
  
  % We now write a column vector containing the square of the 
  % singular values which represents the relative variance of each mode. 
  % Call it D:
  disp('Now we calculate the relative variance of each mode ...');
  Dtmp=S.^2/nt;
  for ii=1:neof
    D(ii)=Dtmp(ii,ii);
  end
  Dtmp=[];clear Dtmp
    
  disp('Now we calculate the principal components of each mode ...');
  PCs1 = data1*U;
  PCs2 = data2*V;
return
