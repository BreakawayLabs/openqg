function eofplot(x,y,t,vv,Dext,Dperc,pcs,xb,yb,vb,str1,frecut)
% EOFPLOT Plots EOF results on one page
%   EOFPLOT(X,Y,T,VV,DEXT,DPERC,PCS,XB,YB,VB,STR1,FRECUT) takes reshaped 
%   eigenfunctions VV(T,Y,X) and plots them with the politoliu colour 
%   scale. The mean fields VB(YB,XB) are also plotted, and the spectrum 
%   of the first 6 PCs. The title for the page is given by STR1;
%   information about the relative strength of the modes is also
%   included. FRECUT is the filtering length in yrs^{-1}.
%
%  v1.0 AH 27/8/2004
%   

%   VERSION LOG
%   v1.0 - created from eofplot1.m by AH, 27/8/04
%        - designed to cope with all aspect ratios.

  close all
  figure(1)
  load politoliu
  
  asp = x(end)/y(end);
  if (asp < 2)
    afig(2)
    np1 = 4;
    np2 = 2;
    npnum = [1 2 3 4 5 7];
  elseif (asp < 4)
    afig(3)
    np1 = 4;
    np2 = 2;
    npnum = [1 2 3 4 5 7];
  else
    afig(2)
    np1 = 9;
    np2 = 1;
    npnum = [1 2 3 4 5 6];
  end
    
  %% Plot some eigenfunctions in real space:
  for ii=1:6
    pp=squeeze(vv(:,:,ii));
    cvir =max(abs(pp(:)));
  
    subplot(np1,np2,npnum(ii))
    contourf(x,y,pp,20)
    shading flat
    caxis([-cvir cvir])
    hold on
    contour(xb,yb,vb,10,'k')
    str = sprintf('Mode %d -  D=%8.4g (%4.2f%%).',ii,Dext(ii),Dperc(ii));
    title(str)    
    set(gca,'dataaspectratio',[1 1 1])
    if (asp < 4)
      if ((ii==2)|(ii==4))
	set(gca,'yticklabel',[])
      else
	ylabel('Y (km)')
      end
      if ((ii==4)|(ii==6))
	xlabel('X (km)')
      else
	set(gca,'xticklabel',[])
      end
    else
      ylabel('Y (km)')
      if (ii==6)
	xlabel('X (km)')
      else
	set(gca,'xticklabel',[])
      end
    end
    
  end
  colormap(politoliu)
  suptmp = sprintf(': These make up %8.4g%%.',sum(Dperc(1:6)));
  s=suptitle([str1,suptmp]);
  set(s,'interpreter','none')

  %% Calculate power specturm of the PCs for plotting
  dt=t(2)-t(1);        %yrs
  nt=length(t);
  fs = 1/dt;           %yrs^-1
  for ii=1:6
    [st(:,ii),con,ft] = pmtm(pcs(:,ii),4,round(nt/2),fs,'adapt');
  end
  
  %% find axes for spectral plot
  [y,ftmax]=min( abs(ft - frecut));
  
  if (asp < 4)
    subplot(224),    
  else
    subplot(313),
  end
  loglog(ft,st(:,1:6))
  grid on
  legend(num2str([1:6]'),3)
  ylabel('Power')
  xlabel('Frequency (yrs^{-1})')
  title('Power Spectra: First 6 PCs')
  axis([min(ft) ft(ftmax) min(min(st(1:ftmax,1:6))) max(max(st(1:ftmax,1:6)))])
  return