function lagcorrcoef(vars,t,inc,str,str2)
% LAGCORRCOEF Lagged correlation between vectors
%   LAGCORRCOEF(VARS,T,INC,STR,STR2) takes a matrix (called
%  VARS) and does a lagged cross-correlation between the
%  columns.
%   T is the time vector
%   INC is the number of points to lag by
%   STR includes variable names for each column.
%   STR2 is the title of each page
%
%  v1.0 AH 2/9/2004

%   VERSION LOG
%   v1.0 - created by AH, 2/9/04

close
figure(1)

dt = t(2)-t(1);
x = [-inc:inc]*dt;
numvar = size(vars,2);

numplots=sum(1:numvar-1);
if (numplots > 20)
  afig(2)
  iplot=4;
  jplot=6;
elseif (numplots > 15)
  afig(2)
  iplot=4;
  jplot=5;
elseif (numplots >12)
  afig(2)
  iplot=3;
  jplot=5;
elseif (numplots > 8)
  afig(2)
  iplot=3;
  jplot=4;
elseif (numplots > 6)
  afig(2)
  iplot=2;
  jplot=4;
else
  afig(1)
  iplot=3;
  jplot=2;
end

npleft=numplots;
jk=0; %% jk is the plot number.
for jj=1:numvar
  for kk=jj+1:numvar
    jk=jk+1;
    if (jk==iplot*jplot+1)
      s=suptitle(str2);
      set(s,'interpreter','none')
      figure
      afig(2)
      jk=1;
      npleft=npleft-iplot*jplot;
    end
    for ii=-inc:inc
      M = squeeze([vars(inc+1:end-inc,jj) vars(inc+1+ii:end-inc+ii,kk)]);
      R=corrcoef(M);
      lcc(inc+1+ii) = real(R(1,2));
    end
    subplot(jplot,iplot,jk),plot(x,lcc)
    if (jk >= (min(iplot*jplot,npleft)-iplot+1))
      xlabel('lag (yrs)')
    end
    axis tight
    grid on
    title([char(str(jj)),' vs ',char(str(kk))])
  end
end
s=suptitle(str2);
set(s,'interpreter','none')
return
