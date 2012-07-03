function qgmonit(base_dir,file_dir,run,fileflag,printflag,frecut)
% QGMONIT  Plot monitoring data from OpenQG output
%   QGMONIT(BASE_DIR,FILE_DIR,RUN,FILEFLAG,PRINTFLAG,FRECUT) takes
%  data from either the 'monit.nc' files of a OpenQG
%  run, or else from the 'allvars.mat' file, plots, 
%  and prints all the data. 
%   BASE_DIR is the directory where data is held. 
%   FILE_DIR is the directory where processed data is written. 
%   RUN is the subdirectory for the data.
%   FILEFLAG may be either 'nc' for using netCDF data 
%  file 'monit.nc', or else 'mat' to use the matlab 
%  data file 'allvars.mat'. 
%   PRINTFLAG should be 1 if you want the plots printed 
%  to pdf files, or 0 otherwise. 
%   FRECUT is a freqeuncy cutoff for filtering the data, 
%  and is in units of years^{-1}.
%
%  v2.0 MW 20/7/2009
%       [Marshall Ward (marshall.ward@anu.edu.au)]

%   VERSION LOG
%   v1.0 - created from monitall.m by AH, 23/6/03
%   v1.1 - adapted for Q-GCM V1.2 --  AH, 5/12/03
%   v1.2 - updated for including "run" -- AH 15/6/04
%   v1.3 - updated from Q-GCM V1.3 -- AH 22/7/04
%   v1.4 - more updates for Q-GCM V1.3, and graphical changes
%                                   -- AH 24/8/04
%   v1.5 - Added spinup flag       -- AH 3/9/04
%   v1.6 - separated data and file directories -- AH 1/11/05

%   v2.0 - The 'nc' flag is no longer supported -- MW 17/7/09

  
tic
disp('PRINTING MONITORING DATA:')
disp('-------------------------')

%% Load data
%if strcmp('nc',fileflag)
%  file1=[base_dir,run,'/','monit.nc'];
%  ncload(file1)
%  ldir = pwd;
%  eval(['cd ',base_dir,run,'/'])
%  input_parameters
%  eval(['cd ',ldir])
%elseif strcmp('mat',fileflag)
if strcmp('mat',fileflag)
  file1=[file_dir,run,'/','allvars.mat'];
  load(file1)
elseif strcmp('spin',fileflag)
  file1=[file_dir,run,'/','spinup_vars.mat'];
  load(file1)
else
  disp(['File type ',fileflag,' unknown'])
  return
end

colorlist={'b-','r-','c-','g-','y-','k-'};

N=length(time);
filgap = ceil(frecut*N*dgnday/365);%% Assumes sampling interval of dgnday days!
dt=dgnday*3600*24;                      %% Time in seconds
fs=365/dgnday;                      %% Frequency in years^-1

%% If it's not an ocean only run, then plot atmosphere data
if ~oceanonly
  figure(1)
  afig(2)
  
  %% We don't plot ekman velocities at P points (wepmat, wapmat)
  
  %% Filter & plot ekman velocities at T points & layer 1 entrainment
  wetmat=onedee_filter(wetmat,fs,0,frecut,1);
  watmat=onedee_filter(watmat,fs,0,frecut,1);
  entmat=onedee_filter(entmat(:,1),fs,0,frecut,1);
  enamat=onedee_filter(enamat(:,1),fs,0,frecut,1);
    
  subplot(811),plot(time,wetmat,colorlist{1},time,watmat,colorlist{2},time,entmat,colorlist{3},time,enamat,colorlist{4})
  axis('tight')
  set(gca,'xticklabel','')
  l1=legend('<w ek>','|w ek|','<e1>','|e1|',1);
  title('Ekman velocities and layer 1 entrainment (m/s)')
  grid on
  
  % No longer plot hmlmat or tmlmat
  
  %% Filter and plot heat content of layers 
  hcmlat=onedee_filter(hcmlat,fs,0,frecut,1);
  for ii=1:nla-1
    etamat(:,ii)=onedee_filter(etamat(:,ii),fs,0,frecut,1);
  end
  heat_aqg = etamat*(tat(1)-tat(2))*rhoat*cpat;
  subplot(812), pr(1)=plot(time,hcmlat,colorlist{1});
  hold on
  leglist{1}='ML';
  for ii=1:nla-1
    pr(ii+1)=plot(time,heat_aqg(:,ii),colorlist{ii+1});
    leglist{ii+1}=['L',num2str(ii),',',num2str(ii+1)];
  end
  axis('tight')
  set(gca,'xticklabel','')
  legend(pr,leglist,1)
  title('Mean heat content (J/m^2)')
  grid on
  clear pr leglist
  
  %% no longer plot convection
  
  %% Now plot radiative fluxes
  arocav=onedee_filter(arocav,fs,0,frecut,1);
  arlaav=onedee_filter(arlaav,fs,0,frecut,1);
  olrtop=onedee_filter(olrtop,fs,0,frecut,1);
  subplot(813),plot(time,arocav,'b-',time,arlaav,'r-',time,olrtop,'c-')
  axis('tight')
  set(gca,'xticklabel','')
  title('Radiative flux (W m^{-2})')
  grid on
  legend('ML over oc','ML over land','TOA',1)
  
  %% Don't plot ermasa
  
  %% no longer plot fractional mass error:
  
  %% No longer plot interface height squared. 
  
  %% kinetic energy
  ketot = sum(ekeiat,2);
  dkedt = gradient(ketot,dt);        %% time derivative of KE (W/m^2)
  for ii=1:nla
    ekeiat(:,ii)=onedee_filter(ekeiat(:,ii),fs,0,frecut,1);
  end
  ketot=onedee_filter(ketot,fs,0,frecut,1);
  subplot(814),pr(1)=plot(time,ketot,colorlist{1});
  hold on
  leglist{1}='Tot.';
  for ii=1:nla
    pr(ii+1)=plot(time,ekeiat(:,ii),colorlist{ii+1});
    leglist{ii+1}=['L',num2str(ii)];
  end
  axis('tight')
  set(gca,'xticklabel','')
  legend(pr,leglist,1)
  title('Average KE (J)')
  grid on
  clear pr leglist
  
  %% Tendency of energetics
  dkedt=onedee_filter(dkedt,fs,0,frecut,1);
  ake1 = -pketat - pkenat;
  ake1=sum(ake1,2);
  ake1=onedee_filter(ake1,fs,0,frecut,1);
  dissat = - sum(ah4dat,2);  %% Dissipation
  dissat=onedee_filter(dissat,fs,0,frecut,1);
  utauat=onedee_filter(-utauat,fs,0,frecut,1);
  subplot(815),plot(time,dkedt,'b-',time,ake1,'r-',time,dissat,'c-',time,utauat','y')
  axis('tight')
  set(gca,'xticklabel','')
  legend('dkedt','source', 'diss','drag',1)
  title('Energy terms (W/m^2)')
  grid on
  
  % average PV
  for ii=1:nla
    qavgat(:,ii)=onedee_filter(qavgat(:,ii),fs,0,frecut,1);
  end
  subplot(816), hold on
  for ii=1:nla
    pr(ii)=plot(time,qavgat(:,ii),colorlist{ii});
    leglist{ii}=['L',num2str(ii)];
  end
  axis('tight')
  set(gca,'xticklabel','')
  legend(pr,leglist,1)
  title('Average PV (s^{-1})')
  grid on
  clear pr leglist
  
  %% Storm track velocity
  for ii=1:nla
    atstval(:,ii)=onedee_filter(atstval(:,ii),fs,0,frecut,1);
  end
  subplot(817), hold on
  for ii=1:nla
    pr(ii)=plot(time,atstval(:,ii),colorlist{ii});
    leglist{ii}=['L',num2str(ii)];
  end
  axis('tight')
  set(gca,'xticklabel','')
  legend(pr,leglist,1)
  title('Max storm track velocity (m/s)')
  grid on
  clear pr leglist
  
  %% Storm track position
  for ii=1:nla
    atstpos(:,ii)=onedee_filter(double(atstpos(:,ii)),fs,0,frecut,1);
  end
  subplot(818), hold on
  for ii=1:nla
    pr(ii)=plot(time,atstpos(:,ii),colorlist{ii});
    leglist{ii}=['L',num2str(ii)];
  end
  axis('tight')
  legend(pr,leglist,1)
  title('Storm track position (gridsquare)')
  grid on
  xlabel('time (years)')
  
  tit=['Atmosphere : ',run];
  s=suptitle(tit);
  set(s,'interpreter','none')
  
  legend(pr,leglist,1)
  clear pr leglist
  
  pause(0.0)    
  if printflag
    if strcmp('spin',fileflag)
      print('-dpdf',[file_dir,run,'/','spinat.pdf'])
    else
      print('-dpdf',[file_dir,run,'/','monat.pdf'])
    end
  end
end

%% If it's not an atmosphere only run then plot all ocean data
if ~atmosonly
  figure(2)
  afig(2)
  
  % don't plot ekman velocities at P points (wepmoc, wapmoc)
  
  %% Filter & plot ekman velocities at T points & entrainment
  wetmoc=onedee_filter(wetmoc,fs,0,frecut,1);
  watmoc=onedee_filter(watmoc,fs,0,frecut,1);
  entmoc=onedee_filter(entmoc,fs,0,frecut,1);
  enamoc=onedee_filter(enamoc,fs,0,frecut,1);
  subplot(811),plot(time,wetmoc,colorlist{1},time,watmoc,colorlist{2},time,entmoc,colorlist{3},time,enamoc,colorlist{4})
  axis('tight')
  legend('<w ek>','|w ek|','<e>','|e|',1)
  title('Ekman velocities and layer 1 entrainment (ms^{-1})')
  set(gca,'xticklabel','')
  grid on
  
  %% Filter and plot heat content of layers 
  for ii=1:nlo-1
    etamoc(:,ii)=onedee_filter(etamoc(:,ii),fs,0,frecut,1);
  end
  heat_oqg = etamoc*(tocc(2)-tocc(1))*rhooc*cpoc;
  tmlmoc=onedee_filter(tmlmoc,fs,0,frecut,1);
  heat_oml = tmlmoc*rhooc*cpoc*hmoc;
  subplot(812), pr(1)=plot(time,heat_oml,colorlist{1});
  hold on
  leglist{1}='ML';
  for ii=1:nlo-1
    pr(ii+1)=plot(time,heat_oqg(:,ii),colorlist{ii+1});
    leglist{ii+1}=['L',num2str(ii),',',num2str(ii+1)];
  end
  axis('tight')
  legend(pr,leglist,1)
  set(gca,'xticklabel','')
  title('Mean heat content (J/m^2)')
  grid on
  clear pr leglist
  
  %% no longer plot convection - unfiltered
  
  if ~oceanonly
    % Plot ocean surface heat flux & heat flux at equatorward boundary,
    slhfav=onedee_filter(slhfav,fs,0,frecut,1);
    oradav=onedee_filter(oradav,fs,0,frecut,1);
    subplot(813),pr(1) = plot(time,slhfav,colorlist{1});
    hold on
    leglist{1} = 'SL';
    pr(2) = plot(time,oradav,colorlist{2});
    leglist{2} = 'rad';
    if (hflxsb)
      ttmads=onedee_filter(ttmads,fs,0,frecut,1);
      ttmdfs=onedee_filter(ttmdfs,fs,0,frecut,1);
      pr(3) = plot(time,ttmads,colorlist{3})
      leglist{3} = 'adv';
      pr(4) = plot(time,ttmdfs,colorlist{4})
      leglist{4} = 'diff';
    elseif (hflxnb)
      ttmadn=onedee_filter(ttmadn,fs,0,frecut,1);
      ttmdfn=onedee_filter(ttmdfn,fs,0,frecut,1);
      pr(3) = plot(time,ttmadn,colorlist{3})
      leglist{3} = 'adv';
      pr(4) = plot(time,ttmdfn,colorlist{4})
      leglist{4} = 'diff';
    end
    axis('tight')
    set(gca,'xticklabel','')
    legend(pr,leglist,1)
    title('Ocean heat flux (W m^{-2})')
    grid on
  end
  %% no longer plot hfmloc
  %% no longer plot tmaooc    
  %% no longer plot vfsmad
  
  % Energy
  ketot = sum(ekeioc,2); 
  dkedt = gradient(ketot,dt);        %% time derivative of KE (W/m^2)
  for ii=1:nlo
    ekeioc(:,ii)=onedee_filter(ekeioc(:,ii),fs,0,frecut,1);
  end
  ketot=onedee_filter(ketot,fs,0,frecut,1);
  subplot(814), pr(1)=plot(time,ketot,colorlist{1});
  hold on
  leglist{1}='Tot';
  for ii=1:nlo
    pr(ii+1)=plot(time,ekeioc(:,ii),colorlist{ii+1});
    leglist{ii+1}=['L',num2str(ii)];
  end
  axis('tight')
  legend(pr,leglist,1)
  title('Average KE (J/m^2)')
  grid on
  set(gca,'xticklabel','')
  clear pr leglist
  disp(['KE is ',num2str(mean(ekeioc))])
  
  % Energy tendency
  dkedt=onedee_filter(dkedt,fs,0,frecut,1);
  pketoc=sum(pketoc,2);
  oke1 = -pketoc - pkenoc;
  oke1=onedee_filter(oke1,fs,0,frecut,1);
  diss2oc = -sum(ah2doc,2);
  diss4oc = -sum(ah4doc,2);
  diss2oc=onedee_filter(diss2oc,fs,0,frecut,1);
  diss4oc=onedee_filter(diss4oc,fs,0,frecut,1);
  utauoc=onedee_filter(utauoc,fs,0,frecut,1);
  btdgoc=onedee_filter(-btdgoc,fs,0,frecut,1);
  subplot(815),plot(time,dkedt,'b-',time,oke1,'r-',time,diss2oc,'c-',time,diss4oc,'y-',time,utauoc','g-',time,btdgoc,'k-')
  axis('tight')
  set(gca,'xticklabel','')
  legend('dkedt','source', 'diss2', 'diss4','input','drag',1)
  title('Energy terms (W/m^2)')
  grid on
  
  % average PV
  for ii=1:nlo
    qavgoc(:,ii)=onedee_filter(qavgoc(:,ii),fs,0,frecut,1);
  end
  subplot(816), hold on
  for ii=1:nlo
    pr(ii)=plot(time,qavgoc(:,ii),colorlist{ii});
    leglist{ii}=['L',num2str(ii)];
  end
  axis('tight')
  set(gca,'xticklabel','')
  legend(pr,leglist,1)
  title('Average PV (s^{-1})')
  grid on
  clear pr leglist
  
  %% ocean jet velocity 
  for ii=1:nlo
    ocjval(:,ii)=onedee_filter(ocjval(:,ii),fs,0,frecut,1);
  end
  subplot(817),
  hold on
  for ii=1:nlo
    pr(ii)=plot(time,ocjval(:,ii),colorlist{ii});
    leglist{ii}=['L',num2str(ii)];
  end
  axis('tight')
  legend(pr,leglist,1)
  title('Max ocean jet velocity (m/s)')
  grid on
  set(gca,'xticklabel','')
  clear pr leglist
  
  %% Don't plot ocjpos - it's flawed.
  
  % Now plotting streamfunction extrema. 
  % Note that streamfunction is relative to southern boundary - not absolute
  for ii=1:nlo
    osfmin(:,ii)=onedee_filter(osfmin(:,ii),fs,0,frecut,1);
    osfmax(:,ii)=onedee_filter(osfmax(:,ii),fs,0,frecut,1);
  end
  subplot(818),hold on
  for ii=1:nlo
    pr(ii)=plot(time,-osfmin(:,ii),colorlist{ii});
    leglist{ii}=['Min, L',num2str(ii)];
    pr(ii+3)=plot(time,osfmax(:,ii),colorlist{ii+3});
    leglist{ii+3}=['Max, L',num2str(ii)];
  end
  axis('tight')
  title('Ocean rel. streamfunction extrema (Sv)')
  legend(pr,leglist,1)
  grid on
  xlabel('time (years)')
  
  tit=['Ocean : ',run];
  s=suptitle(tit);
  set(s,'interpreter','none')
    
  legend(pr,leglist,1)
  clear pr leglist
  
  pause(0.0)  
  if printflag
    if strcmp('spin',fileflag)
      print('-dpdf',[file_dir,run,'/','spinoc.pdf'])
    else
      print('-dpdf',[file_dir,run,'/','monoc.pdf'])
    end
  end
end

t1 = toc;
disp(sprintf('Done (%5.1f sec)',t1));
disp(' ')
return
