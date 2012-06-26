function qgnormceofs(base_dir,run,printflag,frecut)
% QGNORMCEOFS  Find and plot filtered COMPLEX EOFS from Q-GCM run
%   QGNORMCEOFS(BASE_DIR,PRINTFLAG,FRECUT) takes filtered data
%   from Q-GCM (filtered by QGFFTFILT) held in the
%   BASE_DIR and finds CEOFS. 
%   RUN is the subdirectory for the data.
%   PRINTFLAG should be 1 if 
%  you want the plots printed to pdf files, or 0 otherwise.
%   FRECUT is the filtering length in yrs^{-1}.  
%
%  v1.2 AH 13/4/2007

%   VERSION LOG
%   v1.0 - created from qgfiltceofs_v1.5.m by AH, 16/6/04
%   v1.1 - updated for Q-GCM V1.3 - AH 27/8/04
%        - now generalised for all aspect ratios
%   v1.2 - updated for Q-GCM V1.4 - AH 13/4/07 

tic
disp('CALCULATING HILBERT EOFS OF DATA NORMALISED SET:')
disp('------------------------------------------------')
    
% Define incoming and outgoing filenames:
outfile = [base_dir,run,'/','normceofs.mat'];
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
  nxso = ceil(nxto/nso); %% Size of subsampled coordinate vectors 
  nyso = ceil(nyto/nso);  
end  

%% Only save first 12 EOFS. Use this to initialise outfile
MM = [1:12];
save(outfile,'MM')

%% First do atmospheric stuff:
if ~(oceanonly)
  load(infile,'ta','xa','ya')
  nt=length(ta);
  
  if outflat(2)==1
    %% Load filtered data, reshape and store as data
    disp('   - Finding filtered atmosphere pressure data ... ')
    load(infile,'pa1new')
    data = reshape(pa1new,nt,nxsa*nysa);
    clear pa1new
    
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
    
    %% find first MM hilbert eofs and pcs
    [pa1V,pa1Dext,pa1Dperc,pa1pcs] = hilberteof(data,MM(end));
    clear data
    pa1vv=reshape(pa1V,nysa,nxsa,MM(end));
    clear pa1V
    save(outfile,'pa1vv','pa1Dext','pa1Dperc','pa1pcs','-append')
    
    %% Load pa1bar
    load(infile,'pa1bar','xpa','ypa')
    
    str1 = [run,': Normalised atmosphere pressure CEOFS and PCs'];
    ceofplot(xa,ya,ta,pa1vv,pa1Dext,pa1Dperc,pa1pcs,xpa,ypa,pa1bar,str1,frecut)
    if printflag
      print('-dpdf',[base_dir,run,'/','normatpceofs.pdf'])
    end
    clear pa1vv pa1Dext pa1Dperc pa1pcs pa1bar
    
    %% Load filtered data, reshape and store as data
    disp('   - Finding filtered atmosphere height data ... ')
    load(infile,'ha1new')
    data = reshape(ha1new,nt,nxsa*nysa);
    clear ha1new
    
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
    
    %% find first MM hilbert eofs and pcs
    [ha1V,ha1Dext,ha1Dperc,ha1pcs] = hilberteof(data,MM(end));
    clear data
    ha1vv=reshape(ha1V,nysa,nxsa,MM(end));
    clear ha1V
    save(outfile,'ha1vv','ha1Dext','ha1Dperc','ha1pcs','-append')
    
    %% Load pa1bar
    load(infile,'ha1bar')
    
    str1 = [run,': Normalised atmosphere height CEOFS and PCs'];
    ceofplot(xa,ya,ta,ha1vv,ha1Dext,ha1Dperc,ha1pcs,xpa,ypa,ha1bar,str1,frecut)
    if printflag
      print('-dpdf',[base_dir,run,'/','normathceofs.pdf'])
    end
    clear ha1vv ha1Dext ha1Dperc ha1pcs xpa ypa ha1bar
  end
  
  if outflat(1)==1
    %% Load filtered data, reshape and store as data
    disp('   - Finding filtered AST data ... ')
    load(infile,'astnew')
    data = reshape(astnew,nt,nxsa*nysa);
    clear astnew
    
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
    
    %% find first MM hilbert eofs and pcs
    [astV,astDext,astDperc,astpcs] = hilberteof(data,MM(end));
    clear data
    astvv=reshape(astV,nysa,nxsa,MM(end));
    clear astV
    save(outfile,'astvv','astDext','astDperc','astpcs','-append')
      
    %% Load average of variable and store
    load(infile,'astbar','xta','yta')
    
    str1 = 'Normalised AST CEOFS and PCs';
    ceofplot(xa,ya,ta,astvv,astDext,astDperc,astpcs,xta,yta,astbar,str1,frecut)
    if printflag
      print('-dpdf',[base_dir,run,'/','normattceofs.pdf'])
    end
    clear astvv astDext astDperc astpcs xta yta astbar
  end  
  clear xa ya ta
end
  
%% Now do oceanic stuff:
if ~(atmosonly)
  load(infile,'to','xo','yo')
  nt=length(to);
  
  if outfloc(2)==1
    %% Load filtered data, reshape and store as data
    disp('   - Finding filtered ocean pressure data ... ')
    load(infile,'po1new')
    data = reshape(po1new,nt,nxso*nyso);
    clear po1new
    
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
    
    %% find first MM hilbert eofs and pcs
    [po1V,po1Dext,po1Dperc,po1pcs] = hilberteof(data,MM(end));
    clear data
    po1vv=reshape(po1V,nyso,nxso,MM(end));
    clear po1V
    save(outfile,'po1vv','po1Dext','po1Dperc','po1pcs','-append')
    
    %% Load average of variable and store as po1bar
    load(infile,'po1bar','xpo','ypo')
    
    str1 = [run,': Normalised ocean pressure CEOFS and PCs'];
    ceofplot(xo,yo,to,po1vv,po1Dext,po1Dperc,po1pcs,xpo,ypo,po1bar,str1,frecut)
    if printflag
      print('-dpdf',[base_dir,run,'/','normocpceofs.pdf'])
    end
    clear po1vv po1Dext po1Dperc po1pcs po1bar
    
    %% Load filtered data, reshape and store as data
    disp('   - Finding filtered ocean height data ... ')
    load(infile,'ho1new')
    data = reshape(ho1new,nt,nxso*nyso);
    clear ho1new
    
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
    
    %% find first MM eofs and pcs
    [ho1V,ho1Dext,ho1Dperc,ho1pcs] = hilberteof(data,MM(end));
    clear data
    ho1vv=reshape(ho1V,nyso,nxso,MM(end));
    clear ho1V
    save(outfile,'ho1vv','ho1Dext','ho1Dperc','ho1pcs','-append')
    
    %% Load average of variable and store as ho1bar
    load(infile,'ho1bar')
    
    str1 = [run,': Normalised ocean height CEOFS and PCs'];
    ceofplot(xo,yo,to,ho1vv,ho1Dext,ho1Dperc,ho1pcs,xpo,ypo,ho1bar,str1,frecut)
    if printflag
      print('-dpdf',[base_dir,run,'/','normochceofs.pdf'])
    end
    clear ho1vv ho1Dext ho1Dperc ho1pcs xpo ypo ho1bar
  end
  
  if outfloc(1)==1
    %% Load filtered data, reshape and store as data
    disp('   - Finding filtered SST data ... ')
    load(infile,'sstnew')
    data = reshape(sstnew,nt,nxso*nyso);
    clear sstnew
    
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
    
    %% find first MM hilbert eofs and pcs
    [sstV,sstDext,sstDperc,sstpcs] = hilberteof(data,MM(end));
    clear data
    sstvv=reshape(sstV,nyso,nxso,MM(end));
    clear sstV
    save(outfile,'sstvv','sstDext','sstDperc','sstpcs','-append')
    
    %% Load average of variable and store
    load(infile,'sstbar','xto','yto')
    
    str1 = [run,': Normalised SST CEOFS and PCs'];
    ceofplot(xo,yo,to,sstvv,sstDext,sstDperc,sstpcs,xto,yto,sstbar,str1,frecut)
    if printflag
      print('-dpdf',[base_dir,run,'/','normoctceofs.pdf'])
    end
    clear sstvv sstDext sstDperc sstpcs xto yto sstbar
  end
  clear xo yo to
end

t1 = toc;
disp(sprintf('Done (%5.1f sec)',t1));
disp(' ')
return