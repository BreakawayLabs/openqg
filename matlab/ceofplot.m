function ceofplot(x,y,t,vv,Dext,Dperc,pcs,xb,yb,vb,str1,frecut)
% CEOFPLOT Plots CEOF results on one page
%   CEOFPLOT(X,Y,T,VV,DEXT,DPERC,PCS,XB,YB,VB,STR1,FRECUT) takes reshaped 
%   eigenfunctions VV(T,Y,X) and plots them with the politoliu colour 
%   scale. The mean fields VB(YB,XB) are also plotted, and the spectrum 
%   of the first 4 PCs. The title for the page is given by STR1;
%   information about the relative strength of the modes is also
%   included. FRECUT is the filtering length in yrs^{-1}. 
%
%  v1.0 AH 27/8/2004
%   

%   VERSION LOG
%   v1.0 - created from ceofplot1.m by AH, 27/8/04
%        - designed to cope with all aspect ratios.

  close all
  figure(1)
  load politoliu
   
  asp = x(end)/y(end);
  if (asp < 2.5)
    afig(2)
    np1 = 3;
    np2 = 2;
    npnum = [1 2 3 4];
  else
    afig(3)
    np1 = 3;
    np2 = 2;
    npnum = [1 2 3 5];
  end
    
  x = x/1000;
  y = y/1000;
  xb = xb/1000;
  yb = yb/1000;
  %% Plot some eigenfunctions in real space:
  for ii=1:4
    ppr=real(squeeze(vv(:,:,ii)));
    ppi=imag(squeeze(vv(:,:,ii)));
    if (asp<1.5)
      pp = [ppr ppi];
      y1 = y;
      x1 = [x; x+x(end)];
      vb1 = [vb vb];
      xb1 = [xb; xb+xb(end)];
      yb1 = yb;
    else
      pp= [ppr; ppi];
      x1 = x;
      y1 = [y; y+y(end)];
      vb1 = [vb; vb];
      yb1 = [yb; yb+yb(end)];
      xb1=xb;
    end
    cvir = max(abs(pp(:))); 
      
    subplot(np1,np2,npnum(ii))
    contourf(x1,y1,pp,20)
    shading flat
    caxis([-cvir cvir])
    hold on
    contour(xb1,yb1,vb1,10,'k')
    if (asp<1.5)
      l=line([x(end) x(end)],[y(1) y(end)]);
      set(l,'color','k','linewidth',2)
      if (ii>=3)
	g=get(gca,'xtick');
	for kk=1:length(g)
	  if (g(kk)> (x(end)))
	    g(kk) = g(kk) - squeeze(x(end));
	  end
	end
	set(gca,'xticklabel',g)
	xlabel('X (000 km)')
      else
	set(gca,'xticklabel',[])
      end
      if ((ii==1)|(ii==3))
	ylabel('Y (000 km)')
      else
	set(gca,'yticklabel',[])
      end
    elseif (asp<2.5)
      l=line([x(1) x(end)],[y(end) y(end)]);
      set(l,'color','k','linewidth',2)
      if ((ii==1)|(ii==3))
	g=get(gca,'ytick');
	for kk=1:length(g)
	  if (g(kk)> (y(end)))
	    g(kk) = g(kk) - squeeze(y(end));
	  end
	end
	set(gca,'yticklabel',g)
	ylabel('Y (000 km)')
      else
	set(gca,'yticklabel',[])
      end
      if (ii>=3)
	xlabel('X (000 km)')
      else
	set(gca,'xticklabel',[])
      end
    else
      l=line([x(1) x(end)],[y(end) y(end)]);
      set(l,'color','k','linewidth',2)
      if ((ii==1)|(ii==3)|(ii==4))
	g=get(gca,'ytick');
	for kk=1:length(g)
	  if (g(kk)> (y(end)))
	    g(kk) = g(kk) - squeeze(y(end));
	  end
	end
	set(gca,'yticklabel',g)
	ylabel('Y (000 km)')
      else
	set(gca,'yticklabel',[])
      end
      if ((ii==2)|(ii==4))
	xlabel('X (000 km)')
      else
	set(gca,'xticklabel',[])
      end
    end
    str = sprintf('Mode %d - Real; Imag -  D=%8.4g (%4.2f%%).',ii,Dext(ii),Dperc(ii));
    title(str)    
    set(gca,'dataaspectratio',[1 1 1])
    
  end
  colormap(politoliu)
  suptmp = sprintf(': These make up %8.4g%%.',sum(Dperc(1:4)));
  s=suptitle([str1,suptmp]);
  set(s,'interpreter','none')

  %% Calculate power specturm of the PCs for plotting
  dt=t(2)-t(1);        %yrs
  nt=length(t);
  fs = 1/dt;           %yrs^-1
  for ii=1:4
    [st(:,ii),con,ft] = pmtm(pcs(:,ii),4,round(nt/2),fs,'adapt');
  end
  
  %% find axes for spectral plot
  [y,ftmax]=min( abs(ft - frecut));
  
  if (asp > 2.5)
    subplot(224),    
  else
    subplot(313),
  end
  loglog(ft,st(:,1:4))
  grid on
  legend(num2str([1:4]'),3)
  ylabel('Power')
  xlabel('Frequency (yrs^{-1})')
  title('Power Spectra: First 4 Complex PCs')
  axis([min(ft) ft(ftmax) min(min(st(1:ftmax,1:4))) max(max(st(1:ftmax,1:4)))])
  
  
  return