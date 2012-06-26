function IMF=onedee_filter(IM,fs,fc1,fc2,zeropad)
% ONEDEE_FILTER generic 1-D filter for N-D arrays.
%
%   IMF = ONEDEE_FILTER(IM,FS,FC1,FC2,ZEROPAD) filters the first
%   dimension of the matrix IM, sampled at frequencies FS (=
%   1/resolution), passing frequencies between FC1 and FC2
%   (boundaries included) in frequency. If ZEROPAD=1 then
%   zero-padding will be used during the FFT to eliminate the
%   periodicity of the signal. 
%
%   ONEDEE_FILTER allows the implementation of different types of 1-D
%   filters. FC1 can be set to 0, thus resulting in a low-pass
%   filter. FC2 can be set to Inf, thus resulting in a high-pass. For
%   instance,ONEDEE_FILTER(IM,FS,0,FC2,0) is a low-pass filter and
%   ONEDEE_FILTER(IM,FS,FC1,Inf,0) is a high pass filter.
%   
%   ONEDEE_FILTER works by Fourier-transforming, masking in 
%   frequency space and antitransforming.
%
%   v1.2 AH 15/6/2004

%   VERSION History
%   v1.0 - Stolen from Cipo's twodee_filter V1.0.
%   v1.1 - Amended to include filtering of first dimension of ND arrays
%   v1.2 - Amended zero-padding to de-mean and re-mean -- AH 15/6/04

if nargin~=5, disp('Wrong number of arguments'); IMF=[]; return, end

n=size(IM,1);
nx=size(IM);
nd=length(nx);
if nd>3
  disp('Matrix IM can be no more than 3 dimensions'); 
  IMF=[]; 
  return
end
if ( nx(2)==1 )
  nd=1;
end

if zeropad==0
  N = n;
elseif zeropad==1
  N = floor(1.3*n);
  %% Also de-mean signal
  IM_mean = mean(IM);
  if nd==1
    IM = IM - IM_mean;
  elseif nd==2
    for ij = 1:n
      IM(ij,:) = IM(ij,:) - IM_mean;
    end
  elseif nd==3
    for ij = 1:n
      IM(ij,:,:) = IM(ij,:,:) - IM_mean;
    end
  end  
else
  disp('Zeropad should be either 1 or 0'); 
  IMF=[]; 
  return
end

fmax=fs/2;

% fix cutoff frequencies in case they are outside permitted values
fc1=max(0,fc1);
fc2=min(fmax,fc2);

% generate frequency axes
if nd==1
  fx=fmax*freqspace(N,'whole')';
elseif nd==2
  [fx,y]=ndgrid(fmax*freqspace(N,'whole'),1:nx(2));
elseif nd==3
  [fx,y,z]=ndgrid(fmax*freqspace(N,'whole'),1:nx(2),1:nx(3));
end

% transform image to frequency space
FT=fft(IM,N);

% build filter mask
fmask=( ((fx>=fc1 & fx<=fc2) | (fx>=(fs-fc2) & fx<=(fs-fc1))) );

% filter
FT=FT.*fmask;

% transform back to original domain
IMF=real(ifft(FT));
if nd==1
  IMF=IMF(1:n);
elseif nd==2
  IMF=IMF(1:n,:);
elseif nd==3
  IMF=IMF(1:n,:,:);
end

if zeropad==1
  %% re-mean signal
  if nd==1
    IMF = IMF + IM_mean;
  elseif nd==2
    for ij = 1:n
      IMF(ij,:) = IMF(ij,:) + IM_mean;
    end
  elseif nd==3
    for ij = 1:n
      IMF(ij,:,:) = IMF(ij,:,:) + IM_mean;
    end
  end  
end

return
