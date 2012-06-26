function qgprocess(base_dir,file_dir,run_name,files,numfiles)
% QGPROCESS   Post-process data from Q-GCM.
%   QGPROCESS(BASE_DIR,FILE_DIR,RUN_NAME,FILES,NUMFILES) takes data from
%  a sequence of Q-GCM simulations.
%   BASE_DIR is the directory where data is held. 
%   FILE_DIR is the directory where processed data is written. 
%   RUN_NAME is the subdirectory for the data.
%   FILES is a cell array holding the filenames for the 
%  individual run segments.
%   NUMFILES is the number of segments.
%
%   Data is saved in the base directory in a file 
%  called 'allvars.mat'.
%
%  v2.0 MW 20/7/2009
%       [Marshall Ward (marshall.ward@anu.edu.au)]

%  NB: wekpo from avges.nc is not saved

%   VERSION LOG
%   v1.0 - created from process.m by AH, 23/6/03
%   v1.1 - adapted to read in data from Q-GCM V1.2 -- AH, 4/12/03
%   v1.2 - updated for including "run" -- AH 15/6/04
%   v1.3 - updated from Q-GCM V1.3 -- AH 22/7/04
%   v1.4 - finished updating for Q-GCM V1.3 -- AH 24/8/04
%   v1.5 - separated data and file directories -- AH 1/11/05
%
%   v2.0 - NetCDF is now handled by snctools
%        - The dataset from monit.nc is now dynamically generated
%        - v2.0 is *slower* than 1.x, for two reasons:
%        -     (*) excessive use of eval()
%        -     (*) var = [var; (new data)] calls
%        - These will be addressed in a later version
%        -   MW (20/7/09)

tic
disp('PROCESSING DATA:')
disp('----------------')
    
outfile = [file_dir,run_name,'/','allvars.mat'];

%% Check whether area-averaged diagnostics are present
aa = dir( [base_dir,run_name,'/',files{1},'/areas.nc'] );
areaflag = not(isempty(aa));

%% Let's go:
for ii=1:numfiles,   %% Loop through each segment of the run
  ii
  %% First load all the data
  run([base_dir,run_name,'/',files{ii},'/','input_parameters']);

  avges_file = [base_dir,run_name,'/',files{ii},'/avges.nc'];
  avges_info = nc_info(avges_file);
  n_avges = length(avges_info.Dataset);
  
  monit_file = [base_dir,run_name,'/',files{ii},'/monit.nc'];
  monit_info = nc_info(monit_file);
  n_monit = length(monit_info.Dataset);
  
  if (areaflag)
    areas_file = [base_dir,run_name,'/',files{ii},'/areas.nc'];
    areas_info = nc_info(areas_file);
  end
  
  %% load the avges grids
  for jj = 1:n_avges
    if(length(avges_info.Dataset(jj).Size(1)) == 1)
      varName = avges_info.Dataset(jj).Name;
      eval(strcat(avges_info.Dataset(jj).Name,' = nc_varget(avges_file, varName);'));
    end      
  end
  
  if ii==1,   %% If it's first time through
    
    %% save tini from input_parameters.m
    tinisave=tini;

    for jj = 1:n_avges
      if(length(avges_info.Dataset(jj).Size(1)) ~= 1)
        varName = avges_info.Dataset(jj).Name;
        eval(strcat(avges_info.Dataset(jj).Name,' = nc_varget(avges_file, varName);'));
      end
    end
    
    %% save time vector
    time = nc_varget(monit_file,'time');
    nt_file = length(time);            % (# of timesteps in each file)
    nt_total = 1 + numfiles*(nt_file-1);   % (total # of saved timesteps)
    
    %% save all the other timeseries from monit.nc
    %
    %  Using 'eval' is generally a bad idea, and is slowing things down, 
    %  but it helps maintain consistency with later scripts
    %
    %  'time' is loaded twice
    
    for jj = 1:n_monit
      if(monit_info.Dataset(jj).Size(1) == nt_file)
        varName = monit_info.Dataset(jj).Name;
        eval(strcat(monit_info.Dataset(jj).Name,' = nc_varget(monit_file, varName);'));
      end
    end
    
    %% save area averaged time-series from areas.nc
    if (areaflag)
        
      if (~oceanonly)
        atdata = nc_varget(areas_file,'atdata');
%        atdatasave = zeros(nt_total, 5);
%        atdatasave(1:(1+nt_file),:) = nc_varget(areas_file,'atdata');
      end
      if (~atmosonly)
        ocdata = nc_varget(areas_file,'ocdata');
%        ocdatasave = zeros(nt_total, 5);
%        ocdatasave(1:(1+nt_file),:) = nc_varget(areas_file,'ocdata');
      end
    end
    
  else %% if it's not the first time through

    %% add averaged data from avges.nc
    for jj = 1:n_avges
      if(length(avges_info.Dataset(jj).Size(1)) ~= 1)
        varName = avges_info.Dataset(jj).Name;
        eval(strcat(varName,' = ',varName,' + nc_varget(avges_file, varName);'));
      end
    end
    
    %% append area-averaged timeseries from areas.nc
    % Note: snctools starts at index i = 0
    if (areaflag)
      if (~oceanonly)
        atdata = [atdata; nc_varget(areas_file,'atdata',[1 0], [-1 -1]) ];
%        atdatasave((2+(ii-1)*nt_file):(1+ii*nt_file)) = nc_varget(areas_file,'atdata',[1 0], [-1 -1]);
      end
      if (~atmosonly)      
        ocdata = [ocdata; nc_varget(areas_file,'ocdata',[1 0], [-1 -1]) ];
%        ocdatasave((2+(ii-1)*nt_file):(1+ii*nt_file)) = nc_varget(areas_file,'ocdata',[1 0], [-1 -1]);
      end
    end
    
    %% append all the other timeseries from monit.nc
    for jj = 1:n_monit
      if(monit_info.Dataset(jj).Size(1) == nt_file)
        varName = monit_info.Dataset(jj).Name;
        varSize = size(monit_info.Dataset(jj).Size);
        varStart = zeros(varSize); varStart(1) = 1;
        varSpan = -1*ones(varSize);
        % Use this to check output
        % disp(strcat(varName,' = [',varName,';, nc_varget(monit_file, varName, [',num2str(varStart),'], [',num2str(varSpan),'] ) ];'));
        eval(strcat(varName,' = [',varName,';, nc_varget(monit_file, varName, [',num2str(varStart),'], [',num2str(varSpan),'] ) ];'));
      end
    end
  
  end
end

%% Done it all. 

%% copy saved tini value back to tini
tini=tinisave;

%% Clear unneeded data and save to matlab format file.
clear *save ans base_dir files ii jj ldir numfiles
clear avges_file avges_info monit_file monit_info areas_file areas_info
clear n_avges n_monit nt_file nt_total varName varSize varStart varSpan
save(outfile)

t1 = toc;
fprintf('Done (%5.1f sec)',t1);
disp(' ')
return