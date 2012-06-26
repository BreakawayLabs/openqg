function IMF = decim3(IM,ns,nx,ny)
% DECIM3 Two-dimensional decimation function for 3-D matrix
%   DECIM3(IM,ns,nx,ny) takes a 3-D matrix IM in which the first
%   dimension is time. It subsamples in
%   both spatial directions (dimensions 2 & 3) at an interval of NS
%   gridpoints. NX and NY are the matrix sizes
%
%   V1.1  - 29/10/2004

%   VERSION LOG
%   v1.0 - created by AH from decim2.m - 15/6/04
%   v1.1 - added facility for nx,ny not multiples of ns - AH 29/10/04

nt=size(IM,1);
tmp=zeros(nt,ceil(ny/ns),ceil(nx/ns));
nn=zeros(nt,ceil(ny/ns),ceil(nx/ns));
for ii=1:ns
  for jj=1:ns
    wrk = IM(:,jj:ns:ny,ii:ns:nx);
    ws = size(wrk);
    wo = ones(ws);
    nn(1:ws(1),1:ws(2),1:ws(3))=nn(1:ws(1),1:ws(2),1:ws(3))+wo;
    tmp(1:ws(1),1:ws(2),1:ws(3)) = tmp(1:ws(1),1:ws(2),1:ws(3)) + wrk;
  end
end
IMF=tmp./nn;