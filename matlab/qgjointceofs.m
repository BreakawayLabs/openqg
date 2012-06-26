function qgjointceofs(base_dir,run,printflag,frecut)
% QGJOINTCEOFS  Find and plot filtered joint EOFS from Q-GCM run
%   QGJOINTCEOFS(BASE_DIR,RUN,PRINTFLAG,FRECUT) takes filtered data
%  from Q-GCM (filtered by QGFFTFILT) held in the
%  BASE_DIR and finds joint Hilbert EOFS -- uses pa1, ast, sst and po1. 
%   RUN is the subdirectory for the data.
%   PRINTFLAG should be 1 if 
%  you want the plots printed to pdf files, or 0 otherwise.
%   FRECUT is the filtering length in yrs^{-1}.  
%
%  v1.1 AH 13/4/2007

%   VERSION LOG
%   v1.0 - created from qgjointeofs_v1.0.m by AH, 1/9/04
%   v1.1 - updated for Q-GCM V1.4 - AH 13/4/07 

tic
disp('CALCULATING JOINT HILBERT EOFS OF DATA SET:')
disp('-------------------------------------------')
    
% Define incoming and outgoing filenames:
outfile = [base_dir,run,'/','jointceofs.mat'];
infile = [base_dir,run,'/','filtdata.mat'];
matfile = [base_dir,run,'/','allvars.mat'];
  
% Load parameters from files
load(matfile,'oceanonly','atmosonly','outflat','outfloc')
if ~(oceanonly)
  load(matfile,'nxta','nyta')
  load(infile,'nsa')
  nxsa = ceil(nxta/nsa); %% Size of subsampled coordinate vectors 
  nysa = ceil(nyta/nsa);  %%
end
if ~(atmosonly)
  load(matfile,'nxto','nyto')
  load(infile,'nso')
  nxso = ceil(nxto/(2*nso)); %% Size of subsampled coordinate vectors 
  nyso = ceil(nyto/(2*nso));  %%(Half of stored size)
end  

%% Only save first 12 EOFS. Use this to initialise outfile
MM = [1:12];
save(outfile,'MM')

data=[];

%% First do atmospheric stuff:
if ~(oceanonly)
  load(infile,'ta','xa','ya')
  nt=length(ta);
  
  if outflat(2)==1
    %% Load filtered data, reshape and store as data
    disp('  - Finding filtered PA1 data ... ')
    load(infile,'pa1new')
    data = [data reshape(pa1new,nt,nxsa*nysa)];
    clear pa1new
        
    %% Load average of variable and store as pa1bar
    load(infile,'pa1bar','xpa','ypa')
  end
   
  if outflat(1)==1
    %% Load filtered data, reshape and store as data
    disp('  - Finding filtered AST data ... ')
    load(infile,'astnew')
    data = [data reshape(astnew,nt,nxsa*nysa)];
    
    %% Load average of variable and store
    load(infile,'astbar','xta','yta')
  end
end

%% Now do oceanic stuff:
if ~(atmosonly)
  load(infile,'to','xo','yo')
  nt=length(to);
  
  if outfloc(2)==1
    %% Load filtered data, reshape and store as data
    disp('  - Finding filtered ocean pressure data ... ')
    load(infile,'po1new')
    data = [data reshape(po1new(:,1:2:end,1:2:end),nt,nxso*nyso)];
    clear po1new
    
    %% Load average of variable and store as po1bar
    load(infile,'po1bar','xpo','ypo')
  end
  
  if outfloc(1)==1
    %% Load filtered data, reshape and store as data
    disp('  - Finding filtered SST data ... ')
    load(infile,'sstnew')
    data = [data reshape(sstnew(:,1:2:end,1:2:end),nt,nxso*nyso)];
    clear sstnew
        
    %% Load average of variable and store
    load(infile,'sstbar','xto','yto')
  end
end

%% De-mean data
DM = mean(data);
for ii = 1:nt
  data(ii,:) = data(ii,:) - DM;
end
clear DM

%% Normalise data point-by-point
vard = sqrt(sum(data.*data)/nt);
for ii=1:nt
  data(ii,:) = data(ii,:)./vard;
end
clear vard

%% Now find first MM Hilbert EOFs and PCs
[V,Dext,Dperc,pcs] = hilberteof(data,MM(end));
clear data

%% Now regroup data back into original variables
nn=0;
if ~(oceanonly)  
  if outflat(2)==1
    pa1vv = reshape(V(nn+1:nn+nxsa*nysa,:),nysa,nxsa,MM(end));
    nn = nn + nxsa*nysa;
  end
  
  if outflat(1)==1
    astvv = reshape(V(nn+1:nn+nxsa*nysa,:),nysa,nxsa,MM(end));
    nn = nn + nxsa*nysa;
  end
end

if ~(atmosonly)  
  if outfloc(2)==1
    po1vv = reshape(V(nn+1:nn+nxso*nyso,:),nyso,nxso,MM(end));
    nn = nn + nxso*nyso;
  end
  
  if outfloc(1)==1
    sstvv = reshape(V(nn+1:nn+nxso*nyso,:),nyso,nxso,MM(end));
    nn = nn + nxso*nyso;
  end
end

cvir = max(abs(V));
Vvar = sum(V.*V);
clear V

%% Calculate power spectrum of the PCs for plotting
dt=to(2)-to(1);        %yrs
nt=length(to);
fs = 1/dt;           %yrs^-1
for ii=1:6
  [st(:,ii),con,ft] = pmtm(pcs(:,ii),4,round(nt/2),fs,'adapt');
end

%% find axes for spectral plot
[y,ftmax]=min( abs(ft - frecut));

%% Save some data
save(outfile,'Dext','Dperc','pcs','-append')

%% Now plot each variable
if ~(oceanonly)
  if outflat(2)==1
    save(outfile,'pa1vv','-append')
    str1 = [run,': Normalised joint CEOFS and PCs -  PA modes'];
    ceofplot(xa,ya,ta,pa1vv,Dext,Dperc,pcs,xpa,ypa,pa1bar,str1,frecut)
    if printflag
      print('-dpdf',[base_dir,run,'/','jointatpceofs.pdf'])
    end
  end
  
  if outflat(1)==1
    save(outfile,'astvv','-append')
    str1 = [run,': Normalised joint CEOFS and PCs - AST modes'];
    ceofplot(xa,ya,ta,astvv,Dext,Dperc,pcs,xta,yta,astbar,str1,frecut)
    if printflag
      print('-dpdf',[base_dir,run,'/','jointattceofs.pdf'])
    end
  end
end

if ~(atmosonly)
  if outfloc(2)==1
    save(outfile,'po1vv','-append')
    str1 = [run,': Normalised joint CEOFS and PCs - PO modes'];
    ceofplot(xo(1:2:end),yo(1:2:end),to,po1vv,Dext,Dperc,pcs,xpo,ypo,po1bar,str1,frecut)
    if printflag
      print('-dpdf',[base_dir,run,'/','jointocpceofs.pdf'])
    end
  end
  
  if outfloc(1)==1
    save(outfile,'sstvv','-append')
    str1 = [run,': Normalised joint CEOFS and PCs - SST modes'];
    ceofplot(xo(1:2:end),yo(1:2:end),to,sstvv,Dext,Dperc,pcs,xto,yto,sstbar,str1,frecut)
    if printflag
      print('-dpdf',[base_dir,run,'/','jointoctceofs.pdf'])
    end
  end
end

t1 = toc;
disp(sprintf('Done (%5.1f sec)',t1));
disp(' ')
return