function qgfilteofs(base_dir,run,printflag,frecut)
% QGFILTEOFS  Find and plot filtered EOFS from OpenQG run
%   QGFILTEOFS(BASE_DIR,RUN,PRINTFLAG,FRECUT) takes filtered data
%  from OpenQG (filtered by QGFFTFILT) held in the
%  BASE_DIR and finds EOFS. 
%   RUN is the subdirectory for the data.
%   PRINTFLAG should be 1 if 
%  you want the plots printed to pdf files, or 0 otherwise.
%   FRECUT is the filtering length in yrs^{-1}.  
%
%  v1.6 AH 12/4/2007

%   VERSION LOG
%   v1.0 - created from filteofs.m by AH, 23/6/03
%          Originally cribbed from code by RYH
%   v1.1 - used eigs rather than eig for 1st 12 EOFS only -AH 7/8/03
%   v1.2 - altered to cope with 3 layer input in ocean pressure - AH 26/8/03
%   v1.3 - eliminated 3 layer input in ocean pressure
%        - updated for version 1.2
%        - now using standardised eof functions  -- AH 19/3/04
%   v1.4 - added use of "run" 
%        - added ha1/ho1 processing
%        - also now have saved average fields in data file 
%        - Filtered data needs to be de-meaned          - AH 16/6/04
%   v1.5 - updated for Q-GCM V1.3 - AH 27/8/04 
%        - now using eofplot
%   v1.6 - updated for Q-GCM V1.4 - AH 12/4/07 
%   

tic
disp('CALCULATING EOFS OF DATA SET:')
disp('-----------------------------')
    
% Define incoming and outgoing filenames:
outfile=[base_dir,run,'/','filteofs.mat'];
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
    
    %% find first MM eofs and pcs
    [pa1V,pa1Dext,pa1Dperc,pa1pcs] = realeof(data,MM(end));
    clear data
    pa1vv=reshape(pa1V,nysa,nxsa,MM(end));
    clear pa1V
    save(outfile,'pa1vv','pa1Dext','pa1Dperc','pa1pcs','-append')
    
    %% Load pa1bar
    load(infile,'pa1bar','xpa','ypa')
    
    str1 = [run,': Filtered atmosphere pressure EOFS and PCs'];
    eofplot(xa,ya,ta,pa1vv,pa1Dext,pa1Dperc,pa1pcs,xpa,ypa,pa1bar,str1,frecut)
    if printflag
      print('-dpdf',[base_dir,run,'/','filtatpeofs.pdf'])
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
    
    %% find first MM eofs and pcs
    [ha1V,ha1Dext,ha1Dperc,ha1pcs] = realeof(data,MM(end));
    clear data
    ha1vv=reshape(ha1V,nysa,nxsa,MM(end));
    clear ha1V
    save(outfile,'ha1vv','ha1Dext','ha1Dperc','ha1pcs','-append')
    
    %% Load pa1bar
    load(infile,'ha1bar')
    
    str1 = [run,': Filtered atmosphere height EOFS and PCs'];
    eofplot(xa,ya,ta,ha1vv,ha1Dext,ha1Dperc,ha1pcs,xpa,ypa,ha1bar,str1,frecut)
    if printflag
      print('-dpdf',[base_dir,run,'/','filtatheofs.pdf'])
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
    
    %% find first MM eofs and pcs
    [astV,astDext,astDperc,astpcs] = realeof(data,MM(end));
    clear data
    astvv=reshape(astV,nysa,nxsa,MM(end));
    clear astV
    save(outfile,'astvv','astDext','astDperc','astpcs','-append')
    
    %% Load average of variable and store
    load(infile,'astbar','xta','yta')
    
    str1 = [run,': Filtered AST EOFS and PCs'];
    eofplot(xa,ya,ta,astvv,astDext,astDperc,astpcs,xta,yta,astbar,str1,frecut)
    if printflag
      print('-dpdf',[base_dir,run,'/','filtatteofs.pdf'])
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
    
    %% find first MM eofs and pcs
    [po1V,po1Dext,po1Dperc,po1pcs] = realeof(data,MM(end));
    clear data
    po1vv=reshape(po1V,nyso,nxso,MM(end));
    clear po1V
    save(outfile,'po1vv','po1Dext','po1Dperc','po1pcs','-append')
    
    %% Load average of variable and store as po1bar
    load(infile,'po1bar','xpo','ypo')
    
    str1 = [run,': Filtered ocean pressure EOFS and PCs'];
    eofplot(xo,yo,to,po1vv,po1Dext,po1Dperc,po1pcs,xpo,ypo,po1bar,str1,frecut)
    if printflag
      print('-dpdf',[base_dir,run,'/','filtocpeofs.pdf'])
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
    
    %% find first MM eofs and pcs
    [ho1V,ho1Dext,ho1Dperc,ho1pcs] = realeof(data,MM(end));
    clear data
    ho1vv=reshape(ho1V,nyso,nxso,MM(end));
    clear ho1V
    save(outfile,'ho1vv','ho1Dext','ho1Dperc','ho1pcs','-append')
    
    %% Load average of variable and store as ho1bar
    load(infile,'ho1bar')
    
    str1 = [run,': Filtered ocean height EOFS and PCs'];
    eofplot(xo,yo,to,ho1vv,ho1Dext,ho1Dperc,ho1pcs,xpo,ypo,ho1bar,str1,frecut)
    if printflag
      print('-dpdf',[base_dir,run,'/','filtocheofs.pdf'])
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
    
    %% find first MM eofs and pcs
    [sstV,sstDext,sstDperc,sstpcs] = realeof(data,MM(end));
    clear data
    sstvv=reshape(sstV,nyso,nxso,MM(end));
    clear sstV
    save(outfile,'sstvv','sstDext','sstDperc','sstpcs','-append')
    
    %% Load average of variable and store
    load(infile,'sstbar','xto','yto')
    
    str1 = [run,': Filtered SST EOFS and PCs'];
    eofplot(xo,yo,to,sstvv,sstDext,sstDperc,sstpcs,xto,yto,sstbar,str1,frecut)
    if printflag
      print('-dpdf',[base_dir,run,'/','filtocteofs.pdf'])
    end
    clear sstvv sstDext sstDperc sstpcs xto yto sstbar
  end
  clear xo yo to
end

t1 = toc;
disp(sprintf('Done (%5.1f sec)',t1));
disp(' ')
return
  
