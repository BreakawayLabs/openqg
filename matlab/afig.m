function P = afig(X);
% X = 1 - a5 landscape
% X = 2 - a4 portrait
% X = 3 - a4 landscape
% X = 4 - a5 portrait
%% We want to work in centimetres:
  set(gcf,'paperunits','centimeters','paperpositionmode','auto');
  %% First set the paper type and orientation so that we get
  %% the correct pdf output.  While we're there, set W and H of 
  %% the figure on the screen.
  if X==1,
    W=21;
    H=14.8;
    set(gcf,'papertype','A5');
    set(gcf,'paperorientation', 'landscape');
  elseif X==2,
    W=21;
    H=29.7;
    set(gcf,'papertype','A4');
    set(gcf,'paperorientation', 'portrait');
  elseif X==3,
    W=29.7;
    H=21;
    set(gcf,'papertype','A4');
    set(gcf,'paperorientation', 'landscape');
  elseif X==4,
    H=21;
    W=14.8;
    set(gcf,'papertype','A5');
    set(gcf,'paperorientation', 'portrait');
  end 
  %% Specify the size of the image on the screen to get correct
  %% font size etc.
  set(gcf,'units','centimeters')
  set(gcf,'position',[3 0 W H])
  return;
  
