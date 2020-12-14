function [h, c] = morOnePhaseOutputDisplay(ml,which,varargin)

    % main fields
    f(1).str={'U'};           f(1).func=@plotHvel; 
    f(end+1).str={'W'};       f(end).func=@plotVvel;
    f(end+1).str={'P'};       f(end).func=@plotPressure; 
    f(end+1).str={'T'};       f(end).func=@plotPotentialTemperature; 
    f(end+1).str={'eta'};     f(end).func=@plotViscosity;
    f(end+1).str={'e_2'};     f(end).func=@plotStrainRateSecInv;
    f(end+1).str={'A'};       f(end).func=@plotGrainSize;
    f(end+1).str={'Difn'};    f(end).func=@plotDiffusion;
    f(end+1).str={'Disl'};    f(end).func=@plotDislocation;
    f(end+1).str={'Plastic'}; f(end).func=@plotPlasticity;
    % overlays 
    f(end+1).str={'stream'};  f(end).func=@plotStreamlines;
    % residuals
    f(end+1).str={'Ur'};      f(end).func=@plotHvelR; 
    f(end+1).str={'Wr'};      f(end).func=@plotVvelR; 
    f(end+1).str={'Pr'};      f(end).func=@plotPressureR; 
    f(end+1).str={'pr'};      f(end).func=@plotPotentialTemperatureR;
    f(end+1).str={'Ar'};      f(end).func=@plotGrainSizeR;
   
    if nargin<2; helpfunc(f,@morOutputDisplay); h=[]; return; end;
    if ischar(ml); ml=loadMorOutput(ml); end
    if ~isstruct(ml); error('Invalid data structure.'); end;
    
    for i=1:length(f)
        if find(strmatch(which,f(i).str,'exact'))
            func = f(i).func;
            h = func(ml,varargin{:});
            return;
        end
    end
    helpfunc(f,@morOutputDisplay); h=[];
    
function helpfunc(f,fname)
  
   display(sprintf('\n* %s(<struct or filename>, <plotname>, [optional arguments]) *\n',func2str(fname)));
   display('Select among the following plotnames:');
   for i=1:length(f);
      display(sprintf('plotname: %15s \t function: %30s ',f(i).str{1},func2str(f(i).func)));
   end
   
    
%%%%%%%%%%%%%%%%%%% MAIN FIELDS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
function h = plotHvel(A,varargin);
    iin = [2:A.par.ni-3]; jin = [3:A.par.nj-3];
    h.im = imagesc(A.coord.xe(iin),A.coord.zc(jin),A.grid.U(jin,iin)*A.par.U0);
    h.pl = plotLabel(A,'Horizontal solid velocity, cm/yr'); axis image; colorbar;
   
function h = plotVvel(A,varargin);
    iin = [3:A.par.ni-3]; jin = [2:A.par.nj-3];
    h.im = imagesc(A.coord.xc(iin),A.coord.ze(jin),A.grid.W(jin,iin)*A.par.U0);
    h.pl = plotLabel(A,'Vertical solid velocity, cm/yr'); axis image; colorbar;
   
function h = plotPressure(A,varargin);
    iin = [3:A.par.ni-3]; jin = [3:A.par.nj-3];
    h.im = imagesc(A.coord.xc(iin),A.coord.zc(jin),A.grid.P(jin,iin));
    h.pl = plotLabel(A,'Dynamic pressure, dimensionless'); axis image; colorbar;
  
function h = plotPotentialTemperature(A,varargin);
    iin = [3:A.par.ni-3]; jin = [3:A.par.nj-3];
    h.im = imagesc(A.coord.xc(iin),A.coord.zc(jin),A.grid.T(jin,iin)*A.par.Tp);
    h.pl = plotLabel(A,'Potential temperature, deg C'); axis image; colorbar;

function h = plotViscosity(A,varargin);
    iin = [3:A.par.ni-3]; jin = [3:A.par.nj-3];
    h.im = imagesc(A.coord.xc(iin),A.coord.zc(jin),log10(A.grid.eta(jin,iin)*A.par.eta0));
    h.pl = plotLabel(A,'log_{10}(Viscosity), Pa-sec'); axis image; colorbar;
    
function h = plotStrainRateSecInv(A,varargin);
    iin = [3:A.par.ni-3]; jin = [3:A.par.nj-3];
    sec_per_year = 60*60*24*365.25;
    h.im = imagesc(A.coord.xc(iin),A.coord.zc(jin),log10(A.grid.e_2(jin,iin)*sec_per_year));
    h.pl = plotLabel(A,'log_{10}(Strain-rate_{II}), year^{-1}'); axis image; colorbar;
    
function h = plotGrainSize(A,varargin);
    iin = [3:A.par.ni-3]; jin = [3:A.par.nj-3];
    h.im = imagesc(A.coord.xc(iin),A.coord.zc(jin),log10(A.par.grain_size.*exp(A.grid.A(jin,iin))));
    h.pl = plotLabel(A,'log_{10}(Grain size), meters'); axis image; colorbar;
    
function h = plotDislocation(A,varargin);
    iin = [3:A.par.ni-3]; jin = [3:A.par.nj-3];
    h.im = imagesc(A.coord.xc(iin),A.coord.zc(jin),A.grid.dislocation(jin,iin));
    h.pl = plotLabel(A,'Fractional dominence of dislocation viscosity'); axis image; colorbar;
    
function h = plotDiffusion(A,varargin);
    iin = [3:A.par.ni-3]; jin = [3:A.par.nj-3];
    h.im = imagesc(A.coord.xc(iin),A.coord.zc(jin),A.grid.diffusion(jin,iin));
    h.pl = plotLabel(A,'Fractional dominence of diffusion viscosity'); axis image; colorbar;
    
function h = plotPlasticity(A,varargin);
    iin = [3:A.par.ni-3]; jin = [3:A.par.nj-3];
    h.im = imagesc(A.coord.xc(iin),A.coord.zc(jin),A.grid.plastic(jin,iin));
    h.pl = plotLabel(A,'Fractional dominence of plasticity'); axis image; colorbar;
       
 %%%%%%%%%%%%%%%%%%% RESIDUALS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   
function h = plotHvelR(A,varargin);
    iin = [2:A.par.ni-3]; jin = [3:A.par.nj-3];
    h.im = imagesc(A.coord.xe(iin),A.coord.zc(jin),A.residual.U(jin,iin));
    h.pl = plotLabel(A,'Horizontal solid velocity, RESIDUAL'); axis image; colorbar;
   
function h = plotVvelR(A,varargin);
    iin = [3:A.par.ni-3]; jin = [2:A.par.nj-3];
    h.im = imagesc(A.coord.xc(iin),A.coord.ze(jin),A.residual.W(jin,iin));
    h.pl = plotLabel(A,'Vertical solid velocity, RESIDUAL'); axis image; colorbar;
   
function h = plotPressureR(A,varargin);
    iin = [3:A.par.ni-3]; jin = [3:A.par.nj-3];
    h.im = imagesc(A.coord.xc(iin),A.coord.zc(jin),A.residual.P(jin,iin));
    h.pl = plotLabel(A,'Dynamic pressure, RESIDUAL'); axis image; colorbar;
  
function h = plotPotentialTemperatureR(A,varargin);
    iin = [3:A.par.ni-3]; jin = [3:A.par.nj-3];
    h.im = imagesc(A.coord.xc(iin),A.coord.zc(jin),A.residual.T(jin,iin));
    h.pl = plotLabel(A,'Potential temperature, RESIDUAL'); axis image; colorbar;
    
function h = plotGrainSizeR(A,varargin);
    iin = [3:A.par.ni-3]; jin = [3:A.par.nj-3];
    h.im = imagesc(A.coord.xc(iin),A.coord.zc(jin),A.residual.A(jin,iin));
    h.pl = plotLabel(A,'Grain size, RESIDUAL'); axis image; colorbar;
    
    
%%%%%%%%%%%%%%%%%%% OVERLAY %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function h = plotStreamlines(A,varargin);
    hold on;
%    xstart  = linspace(-1,1,35)*0.98*A.par.width/2;
    xstart  = linspace(0,1,15)*0.98*A.par.width;
    zstart  = ones(size(xstart))*A.par.height*0.99;
%     zstart  = ones(size(xstart))*A.par.height*0.45;
    U = interp2(A.coord.xe,A.coord.zc',A.grid.U,A.coord.xc,A.coord.zc');
    W = interp2(A.coord.xc,A.coord.ze',A.grid.W,A.coord.xc,A.coord.zc');
    h = streamline(A.coord.xc,A.coord.zc,U,W,xstart,zstart);
    hold off;
    
%%%%%%%%%%%%%%%%%%% UTILS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function pl = plotLabel(A, tstring)
   pl(1) = title(tstring);
   pl(2) = ylabel('Depth, km'); 
   pl(3) = xlabel('Distance from ridge axis, km');

function hs = fnamestring(A);
   hs = text(A.par.width-A.par.xridge,0,A.filename,'FontSize',12, ...
	     'VerticalAlignment','top','HorizontalAlignment','right',...
	     'color','w','interpreter','none');
