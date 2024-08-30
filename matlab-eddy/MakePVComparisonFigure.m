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

zeta_x = wvt.diffY(wvt.w) - wvt.diffZF(wvt.v); % w_y - v_z
zeta_y = wvt.diffZF(wvt.u) - wvt.diffX(wvt.w);  % u_z - w_x
zeta_z = wvt.diffX(wvt.v) - wvt.diffY(wvt.u);  % v_x - u_y

eta = wvt.eta_true;
APVrv1 = zeta_z;
APVvs1 = - wvt.f * wvt.diffZG(eta);
APVnl1 = zeta_x .* wvt.diffX(eta) + zeta_y .* wvt.diffY(eta) + zeta_z .* wvt.diffZG(eta);
APV1 = wvt.apv;
QGPV1 = wvt.qgpv;

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
APVrv2 = zeta_z;
APVvs2 = - wvt.f * wvt.diffZG(eta);
APVnl2 = zeta_x .* wvt.diffX(eta) + zeta_y .* wvt.diffY(eta) + zeta_z .* wvt.diffZG(eta);
APV2 = wvt.apv;
QGPV2 = wvt.qgpv;


xIndex = floor(wvt.Nx/2);
yIndex = floor(wvt.Ny/2);

C = orderedcolors("gem");

ylimits = [-2000 0];

FigureSize = [50 50 1000 500];
fig1 = figure(position=[50 50 figure_width_2col+8 300*scaleFactor],name='Potential vorticity components',PaperUnits="points");
set(gcf, 'Color', 'w');

tl = tiledlayout(1,2,TileSpacing="tight",Padding="tight");

nexttile
plot([0 0],ylimits,Color='black',LineWidth=1.0); hold on
p1 = plot(squeeze(QGPV1(xIndex,yIndex,:))/wvt.f,wvt.z,Color=C(2,:),LineStyle='--',LineWidth=2);
p2 = plot(squeeze(APV1(xIndex,yIndex,:))/wvt.f,wvt.z,Color=C(2,:),LineStyle='-',LineWidth=2);
p3 = plot(squeeze(QGPV2(xIndex,yIndex,:))/wvt.f,wvt.z,Color=C(1,:),LineStyle='--',LineWidth=2);
p4 = plot(squeeze(APV2(xIndex,yIndex,:))/wvt.f,wvt.z,Color=C(1,:),LineStyle='-',LineWidth=2);
xlim([-2 2])
xlabel('potential vorticity (f_0)','FontSize', figure_axis_label_size, 'FontName', figure_font)
ylabel('depth (m)','FontSize', figure_axis_label_size, 'FontName', figure_font)
ylim(ylimits)
lg = legend([p1 p2 p3 p4],'qgpv, anticylone','apv, anticyclone', 'qgpv, cyclone', 'apv, anticyclone',Location='southeast');
lg.FontName = figure_font;
lg.FontSize = figure_axis_tick_size;


ax = nexttile;
plot([0 0],ylimits,Color='black',LineWidth=1.0); hold on
ax.ColorOrderIndex = 4;
p1 = plot(squeeze(APVrv2(xIndex,yIndex,:))/wvt.f,wvt.z,LineWidth=2);
p2 = plot(squeeze(APVvs2(xIndex,yIndex,:))/wvt.f,wvt.z,LineWidth=2);
p3 = plot(squeeze(APVnl2(xIndex,yIndex,:))/wvt.f,wvt.z,LineWidth=2);
xlim([-2 2])
xlabel('potential vorticity (f_0)','FontSize', figure_axis_label_size, 'FontName', figure_font)
ax.YTickLabel = [];
ylim(ylimits)
lg = legend([p1 p2 p3],'relative vorticity','vortex stretching', 'nonlinear',Location='southeast');
lg.FontName = figure_font;
lg.FontSize = figure_axis_tick_size;

if shouldSavePlots == 1
    print('pv-comparison.eps','-depsc2');
end