function qgfiltcpclags(file_dir,run,printflag,halfw)
% QGFILTcPCLAGS  Find and plot filtered CEOFS from OpenQG run
%   QGFILTPCLAGS(FILE_DIR,RUN,PRINTFLAG,HALFW) takes filtered 
%  CPCs from OpenQG (calculated by QGFILTCEOFS) held in the
%  FILE_DIR and does a lagged cross correlation between a 
%  number of quantities. Uses lagcorrceof.m.
%   RUN is the subdirectory for the data.
%   PRINTFLAG should be 1 if 
%  you want the plots printed to pdf files, or 0 otherwise.
%   HALFW is the half-width of the correlation length in yrs.
%
%  v1.1 AH 12/4/2007

%   VERSION LOG
%   v1.0 - created by AH, 2/9/04
%   v1.1 - updated for q-gcm v1.4 -  AH, 12/4/07

tic
disp('CALCULATING LAGGED CORR. COEFS OF CEOFS:')
disp('----------------------------------------')

% Load parameters from allvars.mat
matfile = [file_dir,run,'/','allvars.mat'];
load(matfile,'oceanonly','atmosonly','outflat','outfloc')
    
ceoffile=[file_dir,run,'/','filtceofs.mat'];
load(ceoffile)
filtfile = [file_dir,run,'/','filtdata.mat'];
load(filtfile,'to','ta')
nto = length(to);
nta = length(ta);

dto = to(2) - to(1);
dta = ta(2) - ta(1);
inc = ceil(halfw/dto);

if nta~=nto
  disp('Sorry, vectors should be equal length');
  return
end

if (atmosonly)
  ast1= squeeze(astpcs(:,1));
  ast2= squeeze(astpcs(:,2));
  ast3= squeeze(astpcs(:,3));
  ast4= squeeze(astpcs(:,4));
  pa1= squeeze(pa1pcs(:,1));
  pa2= squeeze(pa1pcs(:,2));
  pa3= squeeze(pa1pcs(:,3));
  pa4= squeeze(pa1pcs(:,4));  
  MM=[ast1 ast2 ast3 ast4 pa1 pa2 pa3 pa4];
  str={'ast1','ast2','ast3','ast4','pa1','pa2','pa3','pa4'};
  numfigs=2;
elseif (oceanonly)
  po1= squeeze(po1pcs(:,1));
  po2= squeeze(po1pcs(:,2));
  sst1= squeeze(sstpcs(:,1));
  sst2= squeeze(sstpcs(:,2)); 
  MM=[po1 po2 sst1 sst2];
  str={'po1','po2','sst1','sst2'};
  numfigs=1;
else
  po1= squeeze(po1pcs(:,1));
  po2= squeeze(po1pcs(:,2));
  sst1= squeeze(sstpcs(:,1));
  sst2= squeeze(sstpcs(:,2));
  ast1= squeeze(astpcs(:,1));
  ast2= squeeze(astpcs(:,2));
  ast3= squeeze(astpcs(:,3));
  ast4= squeeze(astpcs(:,4));
  pa1= squeeze(pa1pcs(:,1));
  pa2= squeeze(pa1pcs(:,2));
  pa3= squeeze(pa1pcs(:,3));
  pa4= squeeze(pa1pcs(:,4));
  MM=[po1 po2 sst1 sst2 ast1 ast2 ast3 ast4 pa1 pa2 pa3 pa4];
  str={'po1','po2','sst1','sst2','ast1','ast2','ast3','ast4','pa1','pa2','pa3','pa4'};
  numfigs=3;
end

str2=[run,': Hilbert PC lagged correlation coefficients'];
lagcorrcoef(MM,to,inc,str,str2)

if (printflag)
  for ii=1:numfigs
    print('-dpdf',[file_dir,run,'/','filtcpclags',num2str(ii),'.pdf'])
  end
end

t1 = toc;
disp(sprintf('Done (%5.1f sec)',t1));
disp(' ')
return
  
