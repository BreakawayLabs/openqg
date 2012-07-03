%% Matlab script which calls routines to process OpenQG model data
%% Updated to cope with all latest Q-GCM_V1.3 scripts - AH 3/9/04
%% Updated from NetCDF Toolbox to SNCTools - MW 17/7/09
%%

clear all
close all

%% The run base directory:
%
% base_dir is the home location of all OpenQG data files
% file_dir is the output location of all OpenQG plots and other output
% run_name is the name given for the folder of a particular run within
%           these locations
%
% Example:
% base_dir = '/Volumes/GFD/andy/'
% file_dir = '/Users/andy/Data/'
% run_name = 'SO-winds/SO_winds_udiff_lovis_40km/'
%
% Raw OpenQG data is stored in base_dir/run_name/
% Processed data is placed in file_dir/run_name/

base_dir = '../examples/double_gyre_vayu/'
file_dir = '../examples/double_gyre_vayu/'
run_name = 'double_gyre_test' 

aa = dir([file_dir,run_name]);
if (isempty(aa))
  mkdir([file_dir,run_name]);
end

%% Parameters to set::
monfrecut = 5 %% Monitoring Frequency cutoff in yrs^{-1}
frecut = 2 %% Data Frequency cutoff in yrs^{-1}
printflag = 1  %% 1 to print, 0 to not print
nsa = 16   %% spatial subsampling for atmos filtering
nso = 16   %% spatial subsampling for ocean filtering

%% SPINUP DATA:
%% Now go to the base directory and figure out which sub-directories
%% the data is being held in.
% Currently seems to assume that there is always a spinup set of data

aa = dir([base_dir,run_name,'/','spinup_yrs*']);
if(~isempty(aa))
    ldir=pwd;
    eval(['cd ',base_dir,run_name,'/']);
    tmp=ls('-d1','spinup_yrs*');
    eval(['cd ',ldir]);
    dirs = reshape(tmp,18,[])';
    files = cellstr(dirs)
    numfiles = length(files);
    
    qgspinproc(base_dir,file_dir,run_name,files,numfiles)
    qgmonit(base_dir,file_dir,run_name,'spin',printflag,monfrecut)
end

%% FULL DATA
% Now go to the base directory and figure out which sub-directories
% the data is being held in.

ldir=pwd;
eval(['cd ',base_dir,run_name,'/']);
tmp=ls('-d1','yrs*');
eval(['cd ',ldir]);
dirs = reshape(tmp,11,[])';
files = cellstr(dirs)
numfiles = length(files);

%% These ones to show monitoring and averages:
qgprocess(base_dir,file_dir,run_name,files,numfiles)
qgmonit(base_dir,file_dir,run_name,'mat',printflag,monfrecut)
qgplotav(base_dir,file_dir,run_name,'mat',printflag)

%% These ones to look at spatial variability [Optional]
qgfftfilt(base_dir,file_dir,run_name,files,numfiles,frecut,nsa,nso)
qgfilteofs(file_dir,run_name,printflag,frecut)
qgfiltpclags(file_dir,run_name,printflag,1)
qgfiltceofs(file_dir,run_name,printflag,frecut)
qgfiltcpclags(file_dir,run_name,printflag,1)

%% Extra ones for variability and correlations.[Very Optional]
%qgnormeofs(file_dir,run_name,printflag,frecut)
%qgnormceofs(file_dir,run_name,printflag,frecut)
%qgjointeofs(file_dir,run_name,printflag,frecut)
%qgjointceofs(file_dir,run_name,printflag,frecut)
%qgjointsvd(file_dir,run_name,printflag,frecut,'sst','ast')
%qgjointsvd(file_dir,run_name,printflag,frecut,'sst','ha1')
%qgjointsvd(file_dir,run_name,printflag,frecut,'ho1','ast')
%qgjointsvd(file_dir,run_name,printflag,frecut,'ho1','ha1')
%qgnormcca(file_dir,run_name,printflag,frecut,'sst','ast')
%qgnormcca(file_dir,run_name,printflag,frecut,'sst','ha1')
%qgnormcca(file_dir,run_name,printflag,frecut,'ho1','ast')
%qgnormcca(file_dir,run_name,printflag,frecut,'ho1','ha1')

!rm data/temp*.mat

%exit
