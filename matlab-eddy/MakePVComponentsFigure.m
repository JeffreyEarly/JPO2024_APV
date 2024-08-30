%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%% Specify stratification and dimensions
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

shouldSavePlots = 0;


% AMS figure widths, given in picas, converted to points (1 pica=12 points)
scaleFactor = 1;
LoadFigureDefaults
mycolormap = flip(cmocean('balance'));


Lz = 4000;
N0 = 3*2*pi/3600; % buoyancy frequency at the surface, radians/seconds
L_gm = 1300; % thermocline exponential scale, meters
N2 = @(z) N0*N0*exp(2*z/L_gm);

Lx = 750e3;
Ly = 750e3;

Nx = 64;
Ny = 64;
Nz = 40;

wvt = WVTransformHydrostatic([Lx, Ly, Lz], [Nx, Ny, Nz], N2=N2,latitude=35);

wvt.addOperation(EtaTrueOperation());
wvt.addOperation(APVOperation());

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%% Add a cyclonic and anticyclonic eddy
%
% Positive sea-surface anomaly, negative density anomaly (warmer water),
% which means that eta is also negative. That is, if I am a parcel of fluid
% at z=-1000 meters, but my no-motion depth is z=-500 meters, then I need
% to add 500 meters, or really subtract -500 meters per the definition.
%
% Also check out "Vertical structure of mesoscale eddies in the eastern
% South Pacific Ocean: A composite analysis from altimetry and Argo
% profiling floats"
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% "Shallow eddy"
% Density anomaly sits close to the surface
U = -0.30; % m/s
psibar = @(z) (pi*Le*Le/(wvt.Lx*wvt.Ly))*U*(Le/sqrt(2))*exp(1/2)*exp(-(z/He).^2 );
psi = @(x,y,z) U*(Le/sqrt(2))*exp(1/2)*exp(-((x-x0)/Le).^2 -((y-y0)/Le).^2 -(z/He/sqrt(2)).^2 ) - psibar(z);
wvt.setGeostrophicStreamfunction(psi);
rho = shiftdim(wvt.rhobar,-2)+wvt.rho_prime;
if ( any(rho(:) < min(wvt.rhobar)) || any(rho(:) > max(wvt.rhobar)) )
    warning('Invalid initial condition!')
end

zeta_x = wvt.diffY(wvt.w) - wvt.diffZF(wvt.v); % w_y - v_z
zeta_y = wvt.diffZF(wvt.u) - wvt.diffX(wvt.w);  % u_z - w_x
zeta_z = wvt.diffX(wvt.v) - wvt.diffY(wvt.u);  % v_x - u_y

eta = wvt.eta_true;
APVrv = zeta_z;
APVvs = - wvt.f * wvt.diffZG(eta);
APVnl = zeta_x .* wvt.diffX(eta) + zeta_y .* wvt.diffY(eta) + zeta_z .* wvt.diffZG(eta);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%% Figure
%
pvlim = max(abs(APVvs(:)/wvt.f))*[-1 1];

xlimits = [175 575];
ylimits = [-2000 0];
sliceIndex = floor(wvt.Ny/2);

FigureSize = [50 50 1000 500];
fig1 = figure(position=[50 50 figure_width_2col+8 300*scaleFactor],name='Potential vorticity components',PaperUnits="points");
set(gcf, 'Color', 'w');

tl = tiledlayout(1,3,TileSpacing="compact",Padding="tight");

val = squeeze(APVrv(:,sliceIndex,:)/wvt.f).';

ax = nexttile;
set( gca, 'FontSize', figure_axis_tick_size);
pcolor(wvt.x/1000,wvt.z,val); clim(0.1*pvlim), shading interp
colormap(ax,(mycolormap));
cb = colorbar;
cb.Label.String = 'f_0';
cb.Label.FontName = figure_font;
cb.Label.FontSize = figure_axis_tick_size;
cb.Location = "south";
title('\zeta^z',FontSize=figure_title_size,FontName=figure_font)
xlabel('distance (km)','FontSize', figure_axis_label_size, 'FontName', figure_font)
ylabel('depth (m)','FontSize', figure_axis_label_size, 'FontName', figure_font)
xlim(xlimits)
ylim(ylimits)

ax = nexttile;
set( gca, 'FontSize', figure_axis_tick_size);

val = squeeze(APVvs(:,sliceIndex,:)/wvt.f).';

set( gca, 'FontSize', figure_axis_tick_size);
pcolor(wvt.x/1000,wvt.z,val); clim(pvlim), shading interp
colormap(ax,(mycolormap));
cb = colorbar;
cb.Label.String = 'f_0';
cb.Label.FontName = figure_font;
cb.Label.FontSize = figure_axis_tick_size;
cb.Location = "south";
title('-f_0 \partial_z \eta',FontSize=figure_title_size,FontName=figure_font)
xlabel('distance (km)','FontSize', figure_axis_label_size, 'FontName', figure_font)
ax.YTickLabel = [];
xlim(xlimits)
ylim(ylimits)


val = squeeze(APVnl(:,sliceIndex,:)/wvt.f).';

ax = nexttile;
set( gca, 'FontSize', figure_axis_tick_size);
pcolor(wvt.x/1000,wvt.z,val); clim(0.1*pvlim), shading interp
colormap(ax,(mycolormap));
cb = colorbar;
cb.Label.String = 'f_0';
cb.Label.FontName = figure_font;
cb.Label.FontSize = figure_axis_tick_size;
cb.Location = "south";
title('non-linear',FontSize=figure_title_size,FontName=figure_font)
xlabel('distance (km)','FontSize', figure_axis_label_size, 'FontName', figure_font)
ax.YTickLabel = [];
xlim(xlimits)
ylim(ylimits)


% tl.TileSpacing = 'tight';
% tl.Padding = 'tight';

% title(tl,'Anticyclonic eddy');

if shouldSavePlots == 1
    print('pv-components.png','-dpng','-r300');
end