function qgspinproc(base_dir,file_dir,run_name,files,numfiles)
% QGPROCESS   Post-process spinup data from Q-GCM.
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
%  called 'spinup_vars.mat'.
%
%  v2.0 MW 20/7/2009
%       [Marshall Ward (marshall.ward@anu.edu.au)]

%   VERSION LOG
%   v1.0 - created from qgprocess+v1.4.m by AH, 3/9/04
%   v1.1 - separated data and file directories -- AH 1/11/05
%
%   v2.0 - NetCDF is now handled by snctools
%        - The dataset from monit.nc is now dynamically generated
%        - v2.0 is *slower* than 1.x, for two reasons:
%        -     (*) excessive use of eval()
%        -     (*) var = [var; (new data)] calls
%        - These will be addressed in a later version
%        -   MW (20/7/09)

tic
disp('PROCESSING SPINUP DATA:')
disp('-----------------------')
outfile = [file_dir,run_name,'/','spinup_vars.mat'];

%% Let's go:
for ii=1:numfiles,   %% Loop through each segment of the run
  
  %% First load all the data
  run([base_dir,run_name,'/',files{ii},'/','input_parameters']);
  
  monit_file = [base_dir,run_name,'/',files{ii},'/monit.nc'];
  monit_info = nc_info(monit_file);
  n_monit = length(monit_info.Dataset);
%  ldir = pwd;
%  eval(['cd ',base_dir,run_name,'/',files{ii}])
%  input_parameters
%  eval(['cd ',ldir])
%  ncload([base_dir,run_name,'/',files{ii},'/monit.nc'])
  
  if ii==1,   %% If it's first time through
    
    %% save tini
    tinisave=tini;
    
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
        
  else %% if it's not the first time through
    
    %% append all the other timeseries from monit.nc
    for jj = 1:n_monit
      if(monit_info.Dataset(jj).Size(1) == nt_file)
        varName = monit_info.Dataset(jj).Name;
        varSize = size(monit_info.Dataset(jj).Size);
        varStart = zeros(varSize); varStart(1) = 1;
        varSpan = -1*ones(varSize);
%        disp(strcat(varName,' = [',varName,';, nc_varget(monit_file, varName, [',num2str(varStart),'], [',num2str(varSpan),'] ) ];'));
        eval(strcat(varName,' = [',varName,';, nc_varget(monit_file, varName, [',num2str(varStart),'], [',num2str(varSpan),'] ) ];'));
      end
    end
    
  end
end

%% Done it all. 

%% copy saved tini value back to tini
tini=tinisave;

%% Clear unneeded data and save to matlab format file.
clear *save ans base_dir files ii ldir numfiles
clear monit_file monit_info n_monit nt_file nt_total
clear varName varSize varStart varSpan
save(outfile)

t1 = toc;
disp(sprintf('Done (%5.1f sec)',t1));
disp(' ')
return