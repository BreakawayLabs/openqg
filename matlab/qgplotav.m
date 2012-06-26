function qgplotav(base_dir,file_dir,run_name,fileflag,printflag)
% QGPLOTAV  Plot monitoring data from Q-GCM output
%   QGPLOTAV(BASE_DIR,FILE_DIR,RUN_NAME,FILEFLAG,PRINTFLAG) takes
%   data from either the 'avges.nc' files of a Q-GCM
%   run, or else from the 'allvars.mat' file, plots, 
%   and prints important fields. 
%   BASE_DIR is the directory where data is held. 
%   FILE_DIR is the directory where processed data is written. 
%   RUN_NAME is the subdirectory for the data.
%   FILEFLAG may be either 'nc' for using netCDF data 
%  file 'monit.nc', or else 'mat' to use the matlab 
%  data file 'allvars.mat'. 
%   PRINTFLAG should be 1 if you want the plots printed 
%  to pdf files, or 0 otherwise. 
%
%  v2.0 MW 20/7/2009
%       [Marshall Ward (marshall.ward@anu.edu.au)]


%   VERSION LOG
%   v1.0 - created from plot_avges.m by AH, 23/6/03
%   v1.1 - updated to use V1.2 Q_GCM data - AH 5/12/03
%   v1.2 - updated for including "run" -- AH 15/6/04
%   v1.3 - updated for Q-GCM V1.2. Now includes PV, & copes
%          with different ocean aspect ratios -- AH 26/8/04
%   v1.4 - separated data and file directories -- AH 1/11/05
%
%   v2.0 - The 'nc' flag is no longer supported - MW (20/7/09)
%          (The flag is now redundant, but is kept in for now)

tic
disp('PLOTTING MEAN DATA:')
disp('-------------------')
    
%% Load data
% if strcmp('nc',fileflag)
%     
%   file1=[base_dir,'/',run_name,'/','avges.nc'];
%   ncload(file1);
%      
%   run([base_dir,'/',run_name,'/','input_parameters.m']);
% %  ldir = pwd;
% %  eval(['cd ','/',run_name,'/',base_dir])
% %  input_parameters
% %  eval(['cd ','/',run_name,'/',ldir])
% elseif strcmp('mat',fileflag)
if strcmp('mat',fileflag)
  file1=[file_dir,'/',run_name,'/','allvars.mat'];
  load(file1)
else
  disp(['File type ',fileflag,' unknown'])
  return
end
	

%% If it's not an ocean only run, then plot atmosphere data
if ~oceanonly
  figure(1)
  clf
  if (xta(end)/yta(end) < 2.5) 
    afig(2)
  else
    afig(3)
  end
  s=suptitle(run_name);
  set(s,'interpreter','none')
  
  % Plot single layer variables first:
  subplot(521),[cs,h]=contourf(xta,yta,ast,20);
  shading flat
  title('Mixed layer temperature (K)')
  set(gca,'dataaspectratio',[1 1 1])
  colorbar('vert')
  ylabel('Y (km)')
  set(gca,'xticklabel','')
  
  subplot(522),contourf(xta,yta,fnetat,30);
  shading flat
  title('Heat flux (W/m^2)')
  set(gca,'dataaspectratio',[1 1 1])
  colorbar('vert')
  set(gca,'xticklabel','')
  
  subplot(523),pcolor(xpa,yta,uptpat);
  shading flat
  title('<u\prime T\prime> (K.m/s)')
  set(gca,'dataaspectratio',[1 1 1])
  colorbar('vert')
  ylabel('Y (km)')
  set(gca,'xticklabel','')
  
  subplot(524),pcolor(xta,ypa,vptpat);
  shading flat
  title('<v\prime T\prime> (K.m/s)')
  set(gca,'dataaspectratio',[1 1 1])
  colorbar('vert')
  set(gca,'xticklabel','')
  
  subplot(525),[cs,h]=contourf(xpa,ypa,tauxa,20);
  shading flat
  title('\tau^x ( m^2/s^2)')
  set(gca,'dataaspectratio',[1 1 1])
  colorbar('vert')
  ylabel('Y (km)')
  set(gca,'xticklabel','')
  
  subplot(526),[cs,h]=contourf(xpa,ypa,tauya,20);
  shading flat
  title('\tau^y (m^2/s^2)')
  set(gca,'dataaspectratio',[1 1 1])
  colorbar('vert')
  set(gca,'xticklabel','')
  
  for ii=1:2 % Plot bottom 2 layers in pressure
    subplot(5,2,6+ii),[cs,h]=contourf(xpa,ypa,squeeze(pa(ii,:,:)),20);
    shading flat
    if ii==1
      ylabel('Y (km)')
    end
    title(['Layer ',num2str(ii),' pressure (m^2/s^2)'])
    set(gca,'dataaspectratio',[1 1 1])
    colorbar('vert')
    set(gca,'xticklabel','')
  end
  
  for ii=1:2 % Plot bottom 2 layers in pv
    subplot(5,2,8+ii),[cs,h]=contourf(xpa,ypa,squeeze(qa(ii,:,:)),20);
    shading flat
    if ii==1
      ylabel('Y (km)')
    end
    xlabel('X (km)')
    title(['Layer ',num2str(ii),' PV (s^{-1})'])
    set(gca,'dataaspectratio',[1 1 1])
    colorbar('vert')
  end
  
  pause(0.0) 
  if printflag
    print('-dpdf',[file_dir,'/',run_name,'/','avgat.pdf'])
  end
end

%% If it's not an atmosphere only run then plot all ocean data
if ~atmosonly    
  %%Calculate upper & lower layer transport
  for ii=1:nlo
    pf(ii,:,:)=1e-6*hoc(ii)*squeeze(po(ii,:,:) - po(ii,1,1))/fnot;
  end
  pft=squeeze(sum(pf));

  disp(['Transport is ',num2str(-squeeze(pf(:,end,1))'),' Sv']) 
  
  %% Calculate wek
  [txx txy]=gradient(tauxo,dxo);
  [tyx tyy]=gradient(tauxo,dxo);
  wek = (tyx-txy)/fnot;
  
  % Figure out aspect ratio
  aspoc = xpo(end)/ypo(end);
  if (aspoc < 1)
    np1 = 3;
    np2 = 4;
  elseif (aspoc < 3)
    np1 = 4;
    np2 = 3;
  else
    np1 = 6;
    np2 = 2;
  end
  
  figure(2)
  afig(3)
  s=suptitle(run_name);
  set(s,'interpreter','none')
  
  kk=1;
  subplot(np1,np2,kk),[cs,h]=contourf(xto,yto,sst,20);
  shading flat
  title('SST (K)')
  set(gca,'dataaspectratio',[1 1 1])
  if mod(kk,np2)==1
    ylabel('Y (km)')
  end
  if kk<=((np1-1)*np2)
    set(gca,'xticklabel',[])
  else
    xlabel('X (km)')
  end
  pp = get(gca,'position');
  pp(1) = pp(1) + pp(3);
  pp(2) = pp(2) + pp(4)/5;
  pp(3) = 0.01;
  pp(4) = 0.5*pp(4);
  colorbar('position',pp);
  
  kk=2;
  subplot(np1,np2,kk),contourf(xto,yto,fnetoc,20);
  shading flat
  title('Heat flux (W/m^2)')
  set(gca,'dataaspectratio',[1 1 1])
  if mod(kk,np2)==1
    ylabel('Y (km)')
  end
  if kk<=((np1-1)*np2)
    set(gca,'xticklabel',[])
  else
    xlabel('X (km)')
  end
  pp = get(gca,'position');
  pp(1) = pp(1) + pp(3);
  pp(2) = pp(2) + pp(4)/5;
  pp(3) = 0.01;
  pp(4) = 0.5*pp(4);
  colorbar('position',pp);
  
  kk=3;
  subplot(np1,np2,kk),pcolor(xpo,yto,uptpoc);
  shading flat
  title('<u\prime T\prime> (K.m/s)')
  set(gca,'dataaspectratio',[1 1 1])
  if mod(kk,np2)==1
    ylabel('Y (km)')
  end
  if kk<=((np1-1)*np2)
    set(gca,'xticklabel',[])
  else
    xlabel('X (km)')
  end
%  colorbar
  
  kk=4;
  subplot(np1,np2,kk),pcolor(xto,ypo,vptpoc);
  shading flat
  title('<v\prime T\prime> (K.m/s)')
  set(gca,'dataaspectratio',[1 1 1])
  if mod(kk,np2)==1
    ylabel('Y (km)')
  end
  if kk<=((np1-1)*np2)
    set(gca,'xticklabel',[])
  else
    xlabel('X (km)')
  end
%  colorbar
  
  kk=5;
  subplot(np1,np2,kk),[cs,h]=contourf(xpo,ypo,tauxo,20);
  shading flat
  ylabel('Y (km)')
  title('\tau^x (m^2/s^2)')
  set(gca,'dataaspectratio',[1 1 1])
  if mod(kk,np2)==1
    ylabel('Y (km)')
  end
  if kk<=((np1-1)*np2)
    set(gca,'xticklabel',[])
  else
    xlabel('X (km)')
  end
  pp = get(gca,'position');
  pp(1) = pp(1) + pp(3);
  pp(2) = pp(2) + pp(4)/5;
  pp(3) = 0.01;
  pp(4) = 0.5*pp(4);
  colorbar('position',pp);
  
  kk=6;
  subplot(np1,np2,kk),[cs,h]=contourf(xpo,ypo,tauyo,20);
  shading flat
  title('\tau^y (m^2/s^2)')
  set(gca,'dataaspectratio',[1 1 1])
  if mod(kk,np2)==1
    ylabel('Y (km)')
  end
  if kk<=((np1-1)*np2)
    set(gca,'xticklabel',[])
  else
    xlabel('X (km)')
  end
  pp = get(gca,'position');
  pp(1) = pp(1) + pp(3);
  pp(2) = pp(2) + pp(4)/5;
  pp(3) = 0.01;
  pp(4) = 0.5*pp(4);
  colorbar('position',pp);
  
  kk=7;
  subplot(np1,np2,kk),[cs,h]=contourf(xpo,ypo,wek,20);
  shading flat
  title('w_{ek} (m/s)')
  set(gca,'dataaspectratio',[1 1 1])
  if mod(kk,np2)==1
    ylabel('Y (km)')
  end
  if kk<=((np1-1)*np2)
    set(gca,'xticklabel',[])
  else
    xlabel('X (km)')
  end
  pp = get(gca,'position');
  pp(1) = pp(1) + pp(3);
  pp(2) = pp(2) + pp(4)/5;
  pp(3) = 0.01;
  pp(4) = 0.5*pp(4);
  colorbar('position',pp);
  
  kk=8;
  subplot(np1,np2,kk),[cs,h]=contourf(xpo,ypo,pft,20);
  shading flat
  title('Total streamfunction (Sv)')
  set(gca,'dataaspectratio',[1 1 1])
  if mod(kk,np2)==1
    ylabel('Y (km)')
  end
  if kk<=((np1-1)*np2)
    set(gca,'xticklabel',[])
  else
    xlabel('X (km)')
  end
  pp = get(gca,'position');
  pp(1) = pp(1) + pp(3);
  pp(2) = pp(2) + pp(4)/5;
  pp(3) = 0.01;
  pp(4) = 0.5*pp(4);
  colorbar('position',pp);
    
  for ii=1:2
    kk=8+ii;
    subplot(np1,np2,kk),[cs,h]=contourf(xpo,ypo,squeeze(pf(ii,:,:)),20);
    shading flat
    if ii==1
      ylabel('Y (km)')
    end
    title(['Layer ',num2str(ii),' streamfunction (Sv)'])
    set(gca,'dataaspectratio',[1 1 1])
    if mod(kk,np2)==1
      ylabel('Y (km)')
    end
    if kk<=((np1-1)*np2)
      set(gca,'xticklabel',[])
    else
      xlabel('X (km)')
    end
    pp = get(gca,'position');
    pp(1) = pp(1) + pp(3);
    pp(2) = pp(2) + pp(4)/5;
    pp(3) = 0.01;
    pp(4) = 0.5*pp(4);
    colorbar('position',pp);
  end
  
  for ii=1:2
    kk=10+ii;
    subplot(np1,np2,kk),[cs,h]=contourf(xpo,ypo,squeeze(qo(ii,:,:)),20);
    shading flat
    if ii==1
      ylabel('Y (km)')
    end
    title(['Layer ',num2str(ii),' PV (s^{-1})'])
    set(gca,'dataaspectratio',[1 1 1])
    if mod(kk,np2)==1
      ylabel('Y (km)')
    end
    if kk<=((np1-1)*np2)
      set(gca,'xticklabel',[])
    else
      xlabel('X (km)')
    end
    pp = get(gca,'position');
    pp(1) = pp(1) + pp(3);
    pp(2) = pp(2) + pp(4)/5;
    pp(3) = 0.01;
    pp(4) = 0.5*pp(4);
    colorbar('position',pp);
  end
  
  
  pause(0.0)    
  if printflag
    print('-dpdf',[file_dir,'/',run_name,'/','avgoc.pdf'])
  end
end

t1 = toc;
disp(sprintf('Done (%5.1f sec)',t1));
disp(' ')
return
