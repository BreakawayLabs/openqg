function qgfftfilt(base_dir,file_dir,run,files,numfiles,frecut,nsa,nso)
% QGFFTFILT  Filter data from Q-GCM run
%   QGFFTFILT(BASE_DIR,FILE_DIR,RUN,FILES,NUMFILES,FRECUT,NSA,NSO) takes data
%  from a Q-GCM run, and
%     (a) combines segments of data into 1 file;
%     (b) applies an fft filter to the full data set;
%     (c) subsamples the data,subtracts the mean and 
%         saves as a matlab data file.
% 
%   BASE_DIR is the directory where data is held. 
%   FILE_DIR is the directory where processed data is written.
%   RUN is the subdirectory for the data.
%   FILES is a cell array holding the filenames for the 
%  individual run segments.
%   NUMFILES is the number of segments.
%   FRECUT is a freqeuncy cutoff for filtering the data, 
%  and is in units of years^{-1}.
%   NSA is atmospheric spatial subsampling interval, and must be
%  a factor of original ccords.
%   NSO is oceanic spatial subsampling interval.
%
%  v2.0 MW 20/7/2009
%       [Marshall Ward (marshall.ward@anu.edu.au)]

%   VERSION LOG
%   v1.0 - created from fftproc.m by AH, 23/6/03
%   v1.1 -  AH, 6/8/03
%   v1.2 - include full pressure field                 - AH, 20/8/03
%   v1.3 - alterations to update with Q-GCM v1.2       - AH, 19/3/04
%        - also uses function onedee_filter with zero-padding
%   v1.4 - added use of "run" 
%        - added ha1/ho1 processing
%        - also now save average fields in data file 
%        - no longer need to de-mean - now done in onedee_filter 
%        - using new decim3.m                           - AH 15/6/04
%   v1.5 - updated for Q-GCM V1.3 - AH 23/7/04 
%   v1.6 - separated data and file directories -- AH 1/11/05
%   v1.7 - updated for Q-GCM V1.4 - AH 12/4/07 
%
%   v2.0 - netCDF is now handled by snctools - MW 17/7/09
  
tic
disp('FILTERING DATA SET:')
disp('-------------------')
    
% Define incoming and outgoing filenames:
outfile = [file_dir,run,'/','filtdata.mat'];
matfile = [file_dir,run,'/','allvars.mat'];

% Load parameters from allvars.mat
load(matfile,'oceanonly','atmosonly','outfloc','outflat')
if ~(oceanonly)
  load(matfile,'adiday','nska','xta','yta','xpa','ypa','nxta','nyta','gpat') 
  xa = xpa(floor(nsa/2)+1:nsa:nxta+1);
  ya = ypa(floor(nsa/2)+1:nsa:nyta+1);
end
if ~(atmosonly)
  load(matfile,'odiday','nsko','xto','yto','xpo','ypo','nxto','nyto','gpoc')
  xo = xpo(floor(nso/2)+1:nso:nxto);
  yo = ypo(floor(nso/2)+1:nso:nyto);
end
  
% Initialise filtdata.mat with axes:
if (atmosonly)
  save(outfile,'xa','ya','nsa')
elseif (oceanonly)
  save(outfile,'xo','yo','nso')
else
  save(outfile,'xa','ya','nsa','xo','yo','nso')
end

%% Figure out subsampling
if (oceanonly)
  file1 = [base_dir,run,'/',files{1},'/ocpo.nc'];
  time = nc_varget(file1, 'time');
  dt=time(2)-time(1);
  subdeloc=floor(1/(4*frecut*dt));  % Make sure 4 frames per frecut are saved
elseif (atmosonly)
  file1 = [base_dir,run,'/',files{1},'/atpa.nc'];
  time = nc_varget(file1, 'time');
  dt=time(2)-time(1);
  subdelat=floor(1/(4*frecut*dt));  % Make sure 4 frames per frecut are saved
else
  file1 = [base_dir,run,'/',files{1},'/ocpo.nc'];
  time = nc_varget(file1, 'time');
  dt=time(2)-time(1);
  subdeloc=floor(1/(4*frecut*dt));  % Make sure 4 frames per frecut are saved
  file1 = [base_dir,run,'/',files{1},'/atpa.nc'];
  time = nc_varget(file1, 'time');
  dta=time(2)-time(1);
  subdelat=subdeloc*dt/dta;
end

who

if ~(oceanonly)
  ns=(nsa/nska);
  nxs = nxta/nsa; %% Size of subsampled coordinate vectors 
  nys = nyta/nsa;
  nx = floor(nxta/nska); %% size of data coordinate vectors
  ny = floor(nyta/nska);
  diday=adiday;
  
  ns
  nx
  ny
  
  %% First do the time vector:
  disp('Finding time coordinates:')
  timesave=[];
  for kk=1:numfiles,   %% Loop through each segment of the run
    n1 = min(kk,2);
    if outflat(2)==1
      file1 = [base_dir,run,'/',files{kk},'/atpa.nc'];
    elseif outflat(1)==1
      file1 = [base_dir,run,'/',files{kk},'/atast.nc'];
    end
    time = nc_varget(file1, 'time');
    timesave = [timesave; time(n1:end)];
  end
  clear time
  ta = timesave(1:subdelat:end);
  save(outfile,'ta','-append')
  clear ta
  clear timesave
  
  %% Now do pressure & height
  if outflat(2)==1
    tic    
    disp('Filter atmosphere layer 1 pressure & height:')
    
    %% Load averages and store
    disp('     - Find average of data')
    load(matfile,'pa')
    pa1bar = squeeze(pa(1,:,:));
    ha1bar = squeeze(pa(1,:,:) - pa(2,:,:))/gpat(1);
    clear pa
    save(outfile,'ha1bar','pa1bar','xpa','ypa','-append')
    clear ha1bar
    clear pa1bar
    clear xpa
    clear ypa
    
    disp('     - Now loop through each segment, and spatially subsample')
    for kk=1:numfiles,   %% Loop through each segment of the run
      file1 = [base_dir,run,'/',files{kk},'/atpa.nc'];
      disp(['         [ Opening ',file1,' ]'])
      n1 = min(kk,2);
      
      PP1 = nc_varget(file1, 'p', [0 0 0 0], [-1 1 -1 -1]);
      PP2 = nc_varget(file1, 'p', [0 1 0 0], [-1 1 -1 -1]);
      HH1 = (PP1 - PP2)/gpat(1);
      clear PP2
      
      pp1  = decim3(PP1(n1:end,:,:),ns,nx,ny);
      clear PP1
      hh1  = decim3(HH1(n1:end,:,:),ns,nx,ny);
      clear HH1;
      
      pwd
      
      if(~isdir('data'))
          mkdir('data');
      end
      
      save(['data/temp',num2str(kk),'.mat'],'pp1','hh1')
      clear hh1;
      clear pp1;
    end
    tt = toc;
    disp(sprintf('       <  that took %5.1f seconds  >',tt))
    
    tic
    % FILTER
    disp('     - Filter pressure dataset')
    pp1save=[];
    for kk=1:numfiles,   %% Loop through each segment of the run
      load(['data/temp',num2str(kk),'.mat'],'pp1')
      pp1save = [pp1save; pp1];
    end
    clear pp1
    pp1new=onedee_filter(pp1save,1/dt,0,frecut,1);
    clear pp1save
    pa1new = pp1new(1:subdelat:end,:,:);
    clear pp1new
    tt = toc;
    disp(sprintf('       <  that took %5.1f seconds  >',tt))
    save(outfile,'pa1new','-append')
    clear pa1new
    
    tic
    % FILTER
    disp('     - Filter height dataset')
    hh1save=[];
    for kk=1:numfiles,   %% Loop through each segment of the run
      load(['data/temp',num2str(kk),'.mat'],'hh1')
      hh1save = [hh1save; hh1];
    end
    clear hh1
    hh1new=onedee_filter(hh1save,1/dt,0,frecut,1);
    clear hh1save
    ha1new = hh1new(1:subdelat:end,:,:);
    clear hh1new
    tt = toc;
    disp(sprintf('       <  that took %5.1f seconds  >',tt))
    save(outfile,'ha1new','-append')
    clear ha1new
  else
    disp('PA & HA data not saved ...')
  end
  
  if outflat(1)==1
    tic
    disp('Filter AST:')    
    
    disp('     - Find average of data')
    load(matfile,'ast')
    astbar = ast;
    clear ast
    save(outfile,'astbar','xta','yta','-append')
    clear astbar
    clear xta
    clear yta
      
    disp('     - Now loop through each segment, and spatially subsample')
    for kk=1:numfiles,   %% Loop through each segment of the run
      file1 = [base_dir,run,'/',files{kk},'/atast.nc'];
      disp(['         [ Opening ',file1,' ]'])
      ast = nc_varget(file1, 'ast'); %ncload(file1,'ast')
      n1 = min(kk,2);      
      tt1 = decim3(ast(n1:end,:,:),ns,nx,ny);
      save(['data/temp',num2str(kk),'.mat'],'tt1')
    end
    clear ast
    tt = toc;
    disp(sprintf('       <  that took %5.1f seconds  >',tt))
      
    tic
    % FILTER
    disp('     - Filter ast dataset')
    tt1save=[];
    for kk=1:numfiles,   %% Loop through each segment of the run
      load(['data/temp',num2str(kk),'.mat'],'tt1')
      tt1save = [tt1save; tt1];
    end
    clear tt1
    ttnew=onedee_filter(tt1save,1/dt,0,frecut,1);
    clear tt1save
    astnew = ttnew(1:subdelat:end,:,:);
    clear ttnew
    tt = toc;
    disp(sprintf('       <  that took %5.1f seconds  >',tt))
    save(outfile,'astnew','-append')
    clear astnew
  else
    disp('AST data not saved ...')    
  end
end

clear xa ya
who

if ~(atmosonly)
  ns=(nso/nsko);
  nxs = nxto/nso; %% Size of subsampled coordinate vectors 
  nys = nyto/nso;
  nx = floor(nxto/nsko); %% size of data coordinate vectors
  ny = floor(nyto/nsko);
  diday=odiday;
  
  
  ns
  nx
  ny
  
  %% First do the time vector:
  disp('Finding time coordinates:')

  timesave=[];
  for kk=1:numfiles,   %% Loop through each segment of the run
    n1 = min(kk,2);
    if outfloc(2)==1
      file1 = [base_dir,run,'/',files{kk},'/ocpo.nc'];
    elseif outfloc(1)==1
      file1 = [base_dir,run,'/',files{kk},'/ocsst.nc'];
    end
    time = nc_varget(file1, 'time'); % ncload(file1,'time');
    timesave = [timesave; time(n1:end)];
  end
  clear time
  to = timesave(1:subdeloc:end);
  save(outfile,'to','-append')
  clear to
  clear timesave
  
  if outfloc(2)==1
    tic
    disp('Filter ocean layer 1 pressure & height:')
    
    %% Load averages and store
    disp('     - Find average of data')
    load(matfile,'po')
    po1bar = squeeze(po(1,:,:));
    ho1bar = squeeze(po(2,:,:) - po(1,:,:))/gpoc(1);
    clear po
    save(outfile,'ho1bar','po1bar','xpo','ypo','-append')
    clear ho1bar
    clear po1bar
    clear xpo
    clear ypo
      
    disp('     - Now loop through each segment, and spatially subsample')
    for kk=1:numfiles,   %% Loop through each segment of the run
      file1 = [base_dir,run,'/',files{kk},'/ocpo.nc'];
      disp(['         [ Opening ',file1,' ]'])
      n1 = min(kk,2);
      
      PP1 = nc_varget(file1, 'p', [0 0 0 0], [-1 1 -1 -1]);
      PP2 = nc_varget(file1, 'p', [0 1 0 0], [-1 1 -1 -1]);
      HH1 = (PP2 - PP1)/gpoc(1);
      clear PP2
      
      pp1  = decim3(PP1(n1:end,:,:),ns,nx,ny);
      clear PP1
      hh1  = decim3(HH1(n1:end,:,:),ns,nx,ny);
      clear HH1;
      
      save(['data/temp',num2str(kk),'.mat'],'pp1','hh1')
      clear hh1;
      clear pp1;	
    end
    tt = toc;
    disp(sprintf('       <  that took %5.1f seconds  >',tt))
    
    tic
    % FILTER
    disp('     - Filter pressure dataset')
    pp1save=[];
    for kk=1:numfiles,   %% Loop through each segment of the run
      load(['data/temp',num2str(kk),'.mat'],'pp1')
      pp1save = [pp1save; pp1];
    end
    clear pp1
    pp1new=onedee_filter(pp1save,1/dt,0,frecut,1);
    clear pp1save
    po1new = pp1new(1:subdeloc:end,:,:);
    clear pp1new
    tt = toc;
    disp(sprintf('       <  that took %5.1f seconds  >',tt))
    save(outfile,'po1new','-append')
    clear po1new
            
    tic
    % FILTER
    disp('     - Filter height dataset')
    hh1save=[];
    for kk=1:numfiles,   %% Loop through each segment of the run
      load(['data/temp',num2str(kk),'.mat'],'hh1')
      hh1save = [hh1save; hh1];
    end
    clear hh1
    hh1new=onedee_filter(hh1save,1/dt,0,frecut,1);
    clear hh1save
    ho1new = hh1new(1:subdeloc:end,:,:);
    clear hh1new
    tt = toc;
    disp(sprintf('       <  that took %5.1f seconds  >',tt))
    save(outfile,'ho1new','-append')    
    clear ho1new  
  else
    disp('PO & HO data not saved ...')
  end
  
  if outfloc(1)==1
    tic
    disp('Filter SST:')
    
    disp('     - Find average of data')
    load(matfile,'sst')
    sstbar = sst;
    clear sst
    save(outfile,'sstbar','xto','yto','-append')
    clear sstbar
    clear xto
    clear yto
    
    disp('     - Now loop through each segment, and spatially subsample')
    for kk=1:numfiles,   %% Loop through each segment of the run
      file1 = [base_dir,run,'/',files{kk},'/ocsst.nc'];
      disp(['         [ Opening ',file1,' ]'])
      sst = nc_varget(file1, 'sst'); % ncload(file1,'sst')
      n1 = min(kk,2);
      
      tt1 = decim3(sst(n1:end,:,:),ns,nx,ny);
      save(['data/temp',num2str(kk),'.mat'],'tt1')
    end
    clear sst
    tt = toc;
    disp(sprintf('       <  that took %5.1f seconds  >',tt))
    
    tic
    % FILTER
    disp('     - Filter sst dataset')
    tt1save=[];
    for kk=1:numfiles,   %% Loop through each segment of the run
      load(['data/temp',num2str(kk),'.mat'],'tt1')
      tt1save = [tt1save; tt1];
    end
    clear tt1
    ttnew=onedee_filter(tt1save,1/dt,0,frecut,1);
    clear tt1save
    sstnew = ttnew(1:subdeloc:end,:,:);
    clear ttnew
    tt = toc;
    disp(sprintf('       <  that took %5.1f seconds  >',tt))
    save(outfile,'sstnew','-append')
  else
    disp('SST data not saved ...')
  end
end

t1 = toc;
disp(sprintf('Done (%5.1f sec)',t1));
disp(' ')
return
  
