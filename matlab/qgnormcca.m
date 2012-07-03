function qgnormcca(base_dir,run,printflag,frecut,lvar,rvar)
% QGNORMCCA  Find and plot Canonical Correlation Analysis between
% LVAR and RVAR.
%   QGNORMCCA(BASE_DIR,RUN,PRINTFLAG,FRECUT) takes filtered, 
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
%  v1.1 AH 23/7/2004

%   VERSION LOG
%   v1.0 - created using recipe from Bretherton et al. (1992)
%%                                  -- AH 17/6/04
%   v1.1 - updated for Q-GCM V1.3 - AH 23/7/04 

close all

tic
disp(['CALCULATING JOINT ',lvar,'-',rvar,' CCAs:'])
disp('--------------------------------')
    
infile = [base_dir,run,'/','normeofs.mat'];
datfile = [base_dir,run,'/','filtdata.mat'];
outfile=[base_dir,run,'/','normcca_',lvar,'_',rvar,'.mat'];
matfile = [base_dir,run,'/','allvars.mat'];
  
  
% Load parameters from files
load(matfile,'oceanonly','atmosonly')
if ~(oceanonly)
  load(matfile,'nxta','nyta')
  load(datfile,'nsa','ta','xa','ya','xpa','ypa','xta','yta')
  nxsa = ceil(nxta/nsa); %% Size of subsampled coordinate vectors 
  nysa = ceil(nyta/nsa);  %% 
  nt = length(ta);
  dt=ta(2)-ta(1);        %yrs
end
if ~(atmosonly)
  load(matfile,'nxto','nyto')
  load(datfile,'nso','to','xo','yo','xpo','ypo','xto','yto')
  nxso = ceil(nxto/nso); %% Size of subsampled coordinate vectors 
  nyso = ceil(nyto/nso);  %%(Half of stored size) 
  nt = length(to);
  dt=to(2)-to(1);        %yrs
end  

%% Load all EOFs and PCs
load(infile)

%% Load average of variables and store
lvar1 = [lvar,'bar'];
rvar1 = [rvar,'bar'];
load(datfile,lvar1,rvar1)
eval(['lvarbar = ',lvar,'bar;'])
eval(['rvarbar = ',rvar,'bar;'])

%% Figure out aces of left & right variables
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

%% First normalise PCs
eval(['sqd = sqrt(',lvar,'Dext);'])
sqd=sqd';
for ii=1:nt
  eval(['pcsl(ii,:) =',lvar,'pcs(ii,:)./sqd;'])
end

eval(['sqd = sqrt(',rvar,'Dext);'])
sqd=sqd';
for ii=1:nt
  eval(['pcsr(ii,:) =',rvar,'pcs(ii,:)./sqd;'])
end

%% Now build matrices E & F
eval(['E = reshape(',lvar,'vv,nxsl*nysl,MM(end));'])
eval(['F = reshape(',rvar,'vv,nxsr*nysr,MM(end));'])

%% Now build diagonal matrices from Dext
eval(['K = sqrt(diag(',lvar,'Dext));'])
eval(['M = sqrt(diag(',rvar,'Dext));'])

%% Now find cross-correlation matrix
Cab = (pcsl'*pcsr/(nt-1));

%% Take SVD of this
[u,s,v]=svd(Cab);

%% Calculate expansion coeffs - not the normalised ones
eval(['leftpcs = ',lvar,'pcs*u;'])
eval(['rightpcs = ',rvar,'pcs*v;'])

%% Calculate patterns of variability
U = E*K*u;
uu = reshape(U,nysl,nxsl,MM(end));
clear U u
V = F*M*v;
vv = reshape(V,nysr,nxsr,MM(end));
clear V v

%% Calculate power spectrum of the PCs for plotting
fs = 1/dt;           %yrs^-1
for ii=1:8
  [stleft(:,ii),con,ft] = pmtm(leftpcs(:,ii),4,round(nt/2),fs,'adapt');
  [stright(:,ii),con,ft] = pmtm(rightpcs(:,ii),4,round(nt/2),fs,'adapt');
end

%% find axes for spectral plot
[y,ftmax]=min( abs(ft - frecut));

figure(1)
afig(3)
load politoliu

for ii=1:4
  
  %% Now do RH field
  pp=squeeze(vv(:,:,ii));
  
  subplot(5,2,ii*2),contourf(xsr,ysr,pp,20)
  caxis([min(pp(:)) max(pp(:))])
  shading flat
  hold on
  contour(xr,yr,rvarbar,10,'k')
  ylabel('Y (km)')
  set(gca,'dataaspectratio',[1 1 1])
  title(rvar)
    
  %% Now do LH field
  pp=squeeze(uu(:,:,ii));
  subplot(5,2,(ii-1)*2+1),contourf(xsl,ysl,pp,20)
  caxis([min(pp(:)) max(pp(:))])
  shading flat
  hold on
  contour(xl,yl,lvarbar,10,'k')
  ylabel('Y (km)')
  set(gca,'dataaspectratio',[1 1 1])
%  title(lvar)
      
  xte=xsl(end)+ 2000;
  yte=ysl(end) + 500;
%  t=text(xte,yte,sprintf('Mode %d:\n CC = %5.3f',ii,s(ii,ii)),'fontweight','bold','fontsize',12,'horizontalalignment','center');
  title(sprintf(strcat(lvar,',',' Mode %d: CC = %5.3f'),ii,s(ii,ii)));

  colormap(politoliu)
end

ss=suptitle([run,': Canonical correlation patterns, modes 1-4']);
set(ss,'interpreter','none')

subplot(529),loglog(ft,stleft(:,1:4))
grid on
ylabel('Power')
xlabel('Frequency (yrs^{-1})')
title(['Power Spectra: ',lvar])
axis([min(ft) ft(ftmax) min(min(stleft(1:ftmax,1:8))) max(max(stleft(1:ftmax,1:8)))])
legend(num2str([1:4]'),3)
  
subplot(5,2,10),loglog(ft,stright(:,1:4))
grid on
ylabel('Power')
xlabel('Frequency (yrs^{-1})')
title(['Power Spectra: ',rvar])
axis([min(ft) ft(ftmax) min(min(stright(1:ftmax,1:8))) max(max(stright(1:ftmax,1:8)))])
legend(num2str([1:4]'),3)
  
if printflag
  print('-dpdf',[base_dir,run,'/','cca_',lvar,'_',rvar,'_1.pdf'])
end

figure(2)
afig(3)
load politoliu

for ii=1:4
  
  %% Now do RH field
  pp=squeeze(vv(:,:,ii+4));
  
  subplot(5,2,ii*2),contourf(xsr,ysr,pp,20)
  caxis([min(pp(:)) max(pp(:))])
  shading flat
  hold on
  contour(xr,yr,rvarbar,10,'k')
  ylabel('Y (km)')
  set(gca,'dataaspectratio',[1 1 1])
  title(rvar)
    
  %% Now do LH field
  pp=squeeze(uu(:,:,ii+4));
  subplot(5,2,(ii-1)*2+1),contourf(xsl,ysl,pp,20)
  caxis([min(pp(:)) max(pp(:))])
  shading flat
  hold on
  contour(xl,yl,lvarbar,10,'k')
  ylabel('Y (km)')
  set(gca,'dataaspectratio',[1 1 1])
  title(lvar)
  
  xte=xsl(end)+ 2000;
  yte=ysl(end) + 500;
%  t=text(xte,yte,sprintf('Mode %d:\n CC = %5.3f',ii+4,s(ii+4,ii+4)),'fontweight','bold','fontsize',12,'horizontalalignment','center');
  title(sprintf(strcat(lvar,',',' Mode %d: CC = %5.3f'),ii,s(ii,ii)));
  
  colormap(politoliu)
end

ss=suptitle([run,': Canonical correlation patterns, modes 5-8']);
set(ss,'interpreter','none')

subplot(529),loglog(ft,stleft(:,5:8))
grid on
ylabel('Power')
xlabel('Frequency (yrs^{-1})')
title(['Power Spectra: ',lvar])
axis([min(ft) ft(ftmax) min(min(stleft(1:ftmax,1:8))) max(max(stleft(1:ftmax,1:8)))])
legend(num2str([5:8]'),3)
  
subplot(5,2,10),loglog(ft,stright(:,5:8))
grid on
ylabel('Power')
xlabel('Frequency (yrs^{-1})')
title(['Power Spectra: ',rvar])
axis([min(ft) ft(ftmax) min(min(stright(1:ftmax,1:8))) max(max(stright(1:ftmax,1:8)))])
legend(num2str([5:8]'),3)
  
if printflag
  print('-dpdf',[base_dir,run,'/','cca_',lvar,'_',rvar,'_2.pdf'])
end

save(outfile,'lvar','xsl','ysl','leftpcs','uu','xl','yl','lvarbar','rvar','xsr','ysr','rightpcs','vv','xr','yr','rvarbar','dt','nt','s')
t1 = toc;
disp(sprintf('Done (%5.1f sec)',t1));
disp(' ')
return

