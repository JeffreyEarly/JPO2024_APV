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
Le = 120e3;
He = 350;
x0 = (1/2)*max(wvt.x); y0=max(wvt.y)/2; z0=-wvt.Lz/2; Lz = wvt.Lz;

U = 0.30; % m/s
psibar = @(z) (pi*Le*Le/(wvt.Lx*wvt.Ly))*U*(Le/sqrt(2))*exp(1/2)*exp(-(z/He).^2 );
psi = @(x,y,z) U*(Le/sqrt(2))*exp(1/2)*exp(-((x-x0)/Le).^2 -((y-y0)/Le).^2 -(z/He/sqrt(2)).^2 ) - psibar(z);
wvt.setGeostrophicStreamfunction(psi);
rho = shiftdim(wvt.rhobar,-2)+wvt.rho_prime;
if ( any(rho(:) < min(wvt.rhobar)) || any(rho(:) > max(wvt.rhobar)) )
    warning('Invalid initial condition!')
end

rho_total1 = wvt.rho_total;
rho_prime1 = wvt.rho_prime;

U = -0.30; % m/s
psibar = @(z) (pi*Le*Le/(wvt.Lx*wvt.Ly))*U*(Le/sqrt(2))*exp(1/2)*exp(-(z/He).^2 );
psi = @(x,y,z) U*(Le/sqrt(2))*exp(1/2)*exp(-((x-x0)/Le).^2 -((y-y0)/Le).^2 -(z/He/sqrt(2)).^2 ) - psibar(z);
wvt.setGeostrophicStreamfunction(psi);
rho = shiftdim(wvt.rhobar,-2)+wvt.rho_prime;
if ( any(rho(:) < min(wvt.rhobar)) || any(rho(:) > max(wvt.rhobar)) )
    warning('Invalid initial condition!')
end

rho_total2 = wvt.rho_total;
rho_prime2 = wvt.rho_prime;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%% Figure
%

xlimits = [175 575];
ylimits = [-2000 0];
sliceIndex = floor(wvt.Ny/2);

val1 = squeeze(rho_prime1(:,sliceIndex,:)).';
val2 = squeeze(rho_prime2(:,sliceIndex,:)).';

FigureSize = [50 50 1000 500];
fig1 = figure(position=[50 50 figure_width_2col+8 300*scaleFactor],name='Cyclonic & Anticyclonic Eddy',PaperUnits="points");
set(gcf, 'Color', 'w');

tl = tiledlayout(1,3);

ax = nexttile;
set( gca, 'FontSize', figure_axis_tick_size);

zq = linspace(-wvt.Lz,0,500).';

plot(sqrt(N2(zq))*3600/(2*pi),zq,LineWidth=2,Color=0*[1 1 1])
xlabel('frequency (cph)','FontSize', figure_axis_label_size, 'FontName', figure_font)
ylabel('depth (m)','FontSize', figure_axis_label_size, 'FontName', figure_font)
title('N(z)')
xlim([0 3])
ylim(ylimits)

ax = nexttile;
set( gca, 'FontSize', figure_axis_tick_size);

pcolor(wvt.x/1000,wvt.z,val1), shading interp, hold on
clim([-1 1]*max(abs(val1(:))))
cb = colorbar;
cb.Label.String = 'kg/m^{3}';
cb.Label.FontName = figure_font;
cb.Label.FontSize = figure_axis_tick_size;
cb.Location = "south";
contour(wvt.x/1000,wvt.z,squeeze(rho_total1(:,sliceIndex,:)).',linspace(min(rho_total1(:)),max(rho_total1(:)),10),'w','LineWidth',1.5);
contour(wvt.x/1000,wvt.z,squeeze(rho_total1(:,sliceIndex,:)).',linspace(min(rho_total1(:)),max(rho_total1(:)),10),'k','LineWidth',1.0);
colormap(ax,mycolormap);
title('\rho_e',FontSize=figure_title_size,FontName=figure_font)
xlabel('distance (km)','FontSize', figure_axis_label_size, 'FontName', figure_font)
% ylabel('depth (m)')
ax.YTickLabel = [];
xlim(xlimits)
ylim(ylimits)



ax = nexttile;
set( gca, 'FontSize', figure_axis_tick_size);

pcolor(wvt.x/1000,wvt.z,val2), shading interp, hold on
clim([-1 1]*max(abs(val2(:))))
cb = colorbar;
cb.Label.String = 'kg/m^{3}';
cb.Label.FontName = figure_font;
cb.Label.FontSize = figure_axis_tick_size;
cb.Location = "south";
contour(wvt.x/1000,wvt.z,squeeze(rho_total2(:,sliceIndex,:)).',linspace(min(rho_total2(:)),max(rho_total2(:)),10),'w','LineWidth',1.5);
contour(wvt.x/1000,wvt.z,squeeze(rho_total2(:,sliceIndex,:)).',linspace(min(rho_total2(:)),max(rho_total2(:)),10),'k','LineWidth',1.0);
colormap(ax,mycolormap);
title('\rho_e',FontSize=figure_title_size,FontName=figure_font)
xlabel('distance (km)','FontSize', figure_axis_label_size, 'FontName', figure_font)
ax.YTickLabel = [];
xlim(xlimits)
ylim(ylimits)


tl.TileSpacing = 'tight';
tl.Padding = 'tight';

% title(tl,'Anticyclonic eddy');

if shouldSavePlots == 1
    print('cyclonic-anticyclonic.png','-dpng','-r300');
end