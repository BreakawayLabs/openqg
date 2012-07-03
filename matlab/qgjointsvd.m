function qgjointsvd(base_dir,run,printflag,frecut,lvar,rvar)
% QGJOINTSVD  Find and plot joint Singular Value Decomposition 
% between LVAR and RVAR.
%   QGJOINTSVD(BASE_DIR,RUN,PRINTFLAG,FRECUT,LVAR,RVAR) takes filtered, 
%  normalised data from OpenQG (filtered by QGFFTFILT) held in the
%  BASE_DIR and finds EOFS. 
%   RUN is the subdirectory for the data.
%   PRINTFLAG should be 1 if 
%  you want the plots printed to pdf files, or 0 otherwise.
%   FRECUT is the filtering length in yrs^{-1}.  
%   LVAR is a string containing the left-hand variable (usually ocean SST)
%   RVAR is a string containing the right-hand variable 
%  (usually atmospheric height)
% This script will not work if you choose RVAR or LVAR which do not exist!!
%
%  v1.1 AH 13/4/2007

%   VERSION LOG
%   v1.0 - created using recipe from Bretherton et al. (1992)
%        - AH 3/9/04
%   v1.1 - updated for Q-GCM V1.4 - AH 13/4/07 

tic
disp(['CALCULATING JOINT ',lvar,'-',rvar,' SVDs:'])
disp('--------------------------------')
    
datfile = [base_dir,run,'/','filtdata.mat'];
outfile=[base_dir,run,'/','jointsvd_',lvar,'_',rvar,'.mat'];
matfile = [base_dir,run,'/','allvars.mat'];
  
  
% Load parameters from files
load(matfile,'oceanonly','atmosonly')
if ~(oceanonly)
  load(matfile,'nxta','nyta')
  load(datfile,'nsa','ta','xa','ya','xpa','ypa','xta','yta')
  nxsa = ceil(nxta/nsa); %% Size of subsampled coordinate vectors 
  nysa = ceil(nyta/nsa);  %% 
  t=ta;
end
if ~(atmosonly)
  load(matfile,'nxto','nyto')
  load(datfile,'nso','to','xo','yo','xpo','ypo','xto','yto')
  nxso = ceil(nxto/nso); %% Size of subsampled coordinate vectors 
  nyso = ceil(nyto/nso);  %%(Half of stored size) 
  t=to;
end  
nt = length(t);
dt=t(2)-t(1);        %yrs

%% Only save first 12 EOFS. Use this to initialise outfile
MM = [1:12];
save(outfile,'MM')

%% Load average of variables and store
lvar1 = [lvar,'bar'];
rvar1 = [rvar,'bar'];
load(datfile,lvar1,rvar1)
eval(['lvarbar = ',lvar,'bar;'])
eval(['rvarbar = ',rvar,'bar;'])

%% Load actual variables
lvar2 = [lvar,'new'];
rvar2 = [rvar,'new'];
load(datfile,lvar2,rvar2)
eval(['lvarnew = ',lvar,'new;'])
eval(['rvarnew = ',rvar,'new;'])

%% Figure out axes of left & right variables
if strcmp(lvar,'sst') == 1
  xl=xto;yl=yto;
  xsl=xo;ysl=yo;
  nxsl=nxso;nysl=nyso;
elseif  strcmp(lvar,'ho1') == 1
  xl=xpo;yl=ypo;
  xsl=xo;ysl=yo;
  nxsl=nxso;nysl=nyso;
elseif  strcmp(lvar,'po1') == 1
  xl=xpo;yl=ypo;
  xsl=xo;ysl=yo;
  nxsl=nxso;nysl=nyso;
elseif  strcmp(lvar,'ast') == 1
  xl=xta;yl=yta;
  xsl=xa;ysl=ya;
  nxsl=nxsa;nysl=nysa;
elseif  strcmp(lvar,'ha1') == 1
  xl=xpa;yl=ypa;
  xsl=xa;ysl=ya;
  nxsl=nxsa;nysl=nysa;
elseif  strcmp(lvar,'pa1') == 1
  xl=xpa;yl=ypa;
  xsl=xa;ysl=ya;
  nxsl=nxsa;nysl=nysa;
end
if strcmp(rvar,'sst') == 1
  xr=xto;yr=yto;
  xsr=xo;ysr=yo;
  nxsr=nxso;nysr=nyso;
elseif  strcmp(rvar,'ho1') == 1
  xr=xpo;yr=ypo;
  xsr=xo;ysr=yo;
  nxsr=nxso;nysr=nyso;
elseif  strcmp(rvar,'po1') == 1
  xr=xpo;yr=ypo;
  xsr=xo;ysr=yo;
  nxsr=nxso;nysr=nyso;
elseif  strcmp(rvar,'ast') == 1
  xr=xta;yr=yta;
  xsr=xa;ysr=ya;
  nxsr=nxsa;nysr=nysa;
elseif  strcmp(rvar,'ha1') == 1
  xr=xpa;yr=ypa;
  xsr=xa;ysr=ya;
  nxsr=nxsa;nysr=nysa;
elseif  strcmp(rvar,'pa1') == 1
  xr=xpa;yr=ypa;
  xsr=xa;ysr=ya;
  nxsr=nxsa;nysr=nysa;
end
nl = nxsl*nysl;
nr = nxsr*nysr;

data1 = reshape(lvarnew,nt,nl);
clear lvarnew

%% De-mean data1
DM = mean(data1);
for ii = 1:nt
  data1(ii,:) = data1(ii,:) - DM;
end    
clear DM

%% Normalise data point-by-point
vard = sqrt(sum(data1.*data1)/nt);
for ii=1:nt
  data1(ii,:) = data1(ii,:)./vard;
end
clear vard

data2 = reshape(rvarnew,nt,nr);
clear rvarnew

%% De-mean data2
DM = mean(data2);
for ii = 1:nt
  data2(ii,:) = data2(ii,:) - DM;
end    
clear DM

%% Normalise data point-by-point
vard = sqrt(sum(data2.*data2)/nt);
for ii=1:nt
  data2(ii,:) = data2(ii,:)./vard;
end
clear vard

%% Find Joint SVD
[U,V,D,pcs1,pcs2]=jointsvd(data1,data2,MM(end),nt);

%% Reshape back to original
uu = reshape(U,nysl,nxsl,MM(end));
vv = reshape(V,nysr,nxsr,MM(end));
clear U V

%% Save data
save(outfile,'uu','vv','D','pcs1','pcs2','-append')

%% Plot it all Left-hand data first
str1 = [run,': Joint ',lvar,'/',rvar,'SVDs -  ',lvar,' modes'];
eofplot(xsl,ysl,t,uu,D,zeros(1,12),pcs1,xl,yl,lvarbar,str1,frecut)
if printflag
  print('-dpdf',[base_dir,run,'/','svd_',lvar,'_',rvar,'_',lvar,'.pdf'])
end

%% Plot it all right-hand data next
str1 = [run,': Joint ',lvar,'/',rvar,'SVDs -  ',rvar,' modes'];
eofplot(xsr,ysr,t,vv,D,zeros(1,12),pcs2,xr,yr,rvarbar,str1,frecut)
if printflag
  print('-dpdf',[base_dir,run,'/','svd_',lvar,'_',rvar,'_',rvar,'.pdf'])
end

t1 = toc;
disp(sprintf('Done (%5.1f sec)',t1));
disp(' ')
return

