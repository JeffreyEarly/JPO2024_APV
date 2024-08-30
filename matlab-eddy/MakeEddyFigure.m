

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%% Specify stratification and dimensions
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

shouldSavePlots = 0;
stratificationName = 'realistic';
stratificationName = 'exponential';
mycolormap = flip(cmocean('balance'));

if strcmp(stratificationName,'realistic')
    lat0 = 30; lon0 = -179;
    lat0 = 30; lon0 = -120;
    lat0 = -50.0; lon0 = -25.0;
    % lat0 = -39.0; lon0 =  10.5;
    % lat0 = 45.0; lon0 = -150;
    % lat0 = 31.0; lon0 = -59;
    % atlas = VerticalModeAtlas('/Users/jearly/Data/VerticalModeAtlas/VerticalModeAtlas-01.nc');
    atlas = VerticalModeAtlas('/Volumes/MoreStorage/Data/VerticalModeAtlas/VerticalModeAtlas-01.nc');
    [N2_atlas,z_atlas] = atlas.N2(lat0,lon0);
    Lz = max(z_atlas)-min(z_atlas);
    N2 = @(zq) interp1(z_atlas,N2_atlas,zq);
    shouldApplyFilter = 1;
elseif strcmp(stratificationName,'exponential')
    Lz = 4000;
    N0 = 3*2*pi/3600; % buoyancy frequency at the surface, radians/seconds
    L_gm = 1300; % thermocline exponential scale, meters
    N2 = @(z) N0*N0*exp(2*z/L_gm);
    shouldApplyFilter = 0;
else
    error('no');
end

Lx = 750e3;
Ly = 750e3;

Nx = 64;
Ny = 64;
Nz = 40;

wvt = WVTransformHydrostatic([Lx, Ly, Lz], [Nx, Ny, Nz], N2=N2,latitude=35);

wvt.addOperation(EtaTrueOperation());
wvt.addOperation(APVOperation());

%% Plot stratification

zq = linspace(-wvt.Lz,0,500).';
figure(Name="Stratification")
plot(sqrt(N2(zq))*3600/(2*pi),zq,LineWidth=2,Color=0*[1 1 1])
xlabel('frequency (cph)')
ylabel('depth (m)')
title('stratification')
if shouldSavePlots == 1
    print(sprintf('%s-stratification.png',stratificationName),'-dpng','-r300');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%% Plot stratification bounds
%
Le = 120e3;
He = wvt.Lz/5;

Ae = wvt.f*Le*sqrt(2)*exp(1/2)*(1-pi*Le*Le/(wvt.Lx*wvt.Ly))/He/He;
f = @(z) (1/(Ae))*exp((z/He).^2 ) .* ( (L_gm*N0*N0/2)*(1-exp(2*z/L_gm)) - wvt.g*((wvt.rhobar(1) - wvt.rhobar(end))/wvt.rho0) )./z;
figure(name="Stratification bounds")
plot(f(wvt.z),wvt.z/wvt.Lz)

%, hold on, plot(zeros(size(wvt.z)),wvt.z/wvt.Lz)
xlim([-200 200])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%% Add an anticyclonic eddy
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
psi = @(x,y,z) U*(Le/sqrt(2))*exp(1/2)*exp(-((x-x0)/Le).^2 -((y-y0)/Le).^2 -(z/He/sqrt(2)).^2 ) - psibar(z); % .* exp(-(((z-z0)/(Lz/4)).^4)/2);
wvt.setGeostrophicStreamfunction(psi);
if shouldApplyFilter == 1
    [Qk,Ql,Qj] = ExponentialFilter(wvt,nDampedModes=20,shouldAntiAlias=1);
    wvt.A0 = Qj.*wvt.A0;
end

rv = wvt.diffX(wvt.v) - wvt.diffY(wvt.u);
fprintf('min-rv: %.2f f, max-rv: %.2f f\n',min(rv(:))/wvt.f,max(rv(:))/wvt.f)
ssh = wvt.seaSurfaceHeight;
fprintf('min-ssh: %.2f cm, max-ssh: %.2f cm\n',min(ssh(:))*100,max(ssh(:))*100)

rho = shiftdim(wvt.rhobar,-2)+wvt.rho_prime;
if ( any(rho(:) < min(wvt.rhobar)) || any(rho(:) > max(wvt.rhobar)) )
    warning('Invalid initial condition!')
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Make plots
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Notes: May 26th. I like the stratification from lat0 = -50.0; lon0 =
% -25.0; and I also like showing both the positive and negative anomaly,
% which have different errors due to the asymmetry, just like shallow
% water.

xlimits = [175 575];
ylimits = [-4500 0];
sliceIndex = floor(wvt.Ny/2);

% Shouldn't this go to zero?
% mean(mean(trapz(wvt.z,eta,3)/wvt.Lz,2),1)

FigureSize = [50 50 1000 500];
fig1 = figure('Units', 'points', 'Position', FigureSize,'Name','Anticyclonic Eddy');
set(gcf, 'Color', 'w');

tl = tiledlayout(1,4,"TileSpacing","compact");

eta = wvt.eta_true;
rho_total = wvt.rho_total;
val = squeeze(eta(:,sliceIndex,:)).';

ax = nexttile;
plot(sqrt(N2(zq))*3600/(2*pi),zq,LineWidth=2,Color=0*[1 1 1])
xlabel('frequency (cph)')
ylabel('depth (m)')
title('N(z)')
xlim([0 3])
ylim(ylimits)

ax = nexttile;
pcolor(wvt.x/1000,wvt.z,val), shading interp, hold on
clim([-1 1]*max(abs(val(:))))
cb = colorbar;
cb.Label.String = 'm';
cb.Location = "south";
contour(wvt.x/1000,wvt.z,squeeze(rho_total(:,sliceIndex,:)).',linspace(min(rho_total(:)),max(rho_total(:)),10),'w','LineWidth',1.5);
contour(wvt.x/1000,wvt.z,squeeze(rho_total(:,sliceIndex,:)).',linspace(min(rho_total(:)),max(rho_total(:)),10),'k','LineWidth',1.0);
colormap(ax,mycolormap);
title('\eta')
xlabel('distance (km)')
% ylabel('depth (m)')
ax.YTickLabel = [];
xlim(xlimits)
ylim(ylimits)

eta_e = wvt.eta;
val2 = squeeze(eta_e(:,sliceIndex,:)).';

ax = nexttile;
pcolor(wvt.x/1000,wvt.z,val2), shading interp, hold on
clim([-1 1]*max(abs(val(:))))
cb = colorbar;
cb.Label.String = 'm';
cb.Location = "south";
contour(wvt.x/1000,wvt.z,squeeze(rho_total(:,sliceIndex,:)).',linspace(min(rho_total(:)),max(rho_total(:)),10),'w','LineWidth',1.5);
contour(wvt.x/1000,wvt.z,squeeze(rho_total(:,sliceIndex,:)).',linspace(min(rho_total(:)),max(rho_total(:)),10),'k','LineWidth',1.0);
colormap(ax,mycolormap);
title('\eta_e')
xlabel('distance (km)')
ax.YTickLabel = [];
xlim(xlimits)
ylim(ylimits)

eta_error_pct = eta_e./eta;
eta_error_pct(abs(eta)<1) = 1;
eta_error_pct = 100*(eta_error_pct-1);

val = squeeze(eta_error_pct(:,sliceIndex,:)).';

ax = nexttile;
pcolor(wvt.x/1000,wvt.z,val), shading interp
clim([-1 1]*max(abs(val(:))))
cb = colorbar;
cb.Label.String = '%';
cb.Location = "south";
colormap(ax,(mycolormap));
title('\eta_e error')
ax.YTickLabel = [];
xlabel('distance (km)')
xlim(xlimits)
ylim(ylimits)

tl.TileSpacing = 'compact';
tl.Padding = 'compact';

title(tl,'Anticyclonic eddy');

if shouldSavePlots == 1
    print(sprintf('%s-anticyclonic.png',stratificationName),'-dpng','-r300');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%% Add a cyclonic eddy
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% "Shallow eddy"
% Density anomaly sits close to the surface
U = -0.30; % m/s
psibar = @(z) (pi*Le*Le/(wvt.Lx*wvt.Ly))*U*(Le/sqrt(2))*exp(1/2)*exp(-(z/He).^2 );
psi = @(x,y,z) U*(Le/sqrt(2))*exp(1/2)*exp(-((x-x0)/Le).^2 -((y-y0)/Le).^2 -(z/He/sqrt(2)).^2 ) - psibar(z); %.* exp(-(((z-z0)/(Lz/4)).^4)/2);
wvt.setGeostrophicStreamfunction(psi);
if shouldApplyFilter == 1
    [Qk,Ql,Qj] = ExponentialFilter(wvt,nDampedModes=20,shouldAntiAlias=1);
    wvt.A0 = Qj.*wvt.A0;
end

rv = wvt.diffX(wvt.v) - wvt.diffY(wvt.u);
fprintf('min-rv: %.2f f, max-rv: %.2f f\n',min(rv(:))/wvt.f,max(rv(:))/wvt.f)
ssh = wvt.seaSurfaceHeight;
fprintf('min-ssh: %.2f cm, max-ssh: %.2f cm\n',min(ssh(:))*100,max(ssh(:))*100)
if ( any(rho(:) < min(wvt.rhobar)) || any(rho(:) > max(wvt.rhobar)) )
    warning('Invalid initial condition!')
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Make plots
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Notes: May 26th. I like the stratification from lat0 = -50.0; lon0 =
% -25.0; and I also like showing both the positive and negative anomaly,
% which have different errors due to the asymmetry, just like shallow
% water.


eta = wvt.eta_true;
eta_e = wvt.eta;
eta_diff = eta-eta_e; 

fprintf('Total volume-integrated density perturbation:\n')
fprintf('\tint N^2 eta_e dV: %g\n',mean(mean(trapz(wvt.z,eta_e.*shiftdim(wvt.N2,-2),3)/wvt.Lz,2),1)/sqrt(mean(mean(trapz(wvt.z,(eta_e.*shiftdim(wvt.N2,-2)).^2,3)/wvt.Lz,2),1)))
fprintf('\tint rho_e dV: %g\n',mean(mean(trapz(wvt.z,wvt.rho_prime,3)/wvt.Lz,2),1))
fprintf('\tint eta_true dV: %g\n',mean(mean(trapz(wvt.z,eta,3)/wvt.Lz,2),1))

mean(mean(trapz(wvt.z,eta_e,3)/wvt.Lz,2),1)

%%
eta_bar = squeeze(mean(mean(eta,2),1));
eta_e_bar = squeeze(mean(mean(eta_e,2),1));
N2eta_e_bar = squeeze(mean(mean(eta_e.*shiftdim(wvt.N2,-2)/wvt.N2(end),2),1));
figure(Name="Horizontally averaged density"), plot(eta_e_bar,wvt.z,LineWidth=2), hold on
plot(eta_bar,wvt.z,LineWidth=2)
plot(N2eta_e_bar,wvt.z,LineWidth=2)
legend('eta_{e}','eta_{true}','N2*eta_{e}')
title('horizontally averaged')
%%
% eta_e = eta_e+0.5*shiftdim(wvt.dLnN2,-2).*eta_e.^2; % second order correction

mean(mean(trapz(wvt.z,eta_e.*shiftdim(wvt.N2,-2),3)/wvt.Lz,2),1)
% Shouldn't this go to zero?
mean(mean(trapz(wvt.z,eta,3)/wvt.Lz,2),1)

rho_total = wvt.rho_total;
sliceIndex = floor(wvt.Ny/2);

figure(name="Cyclonic Eddy")
tl = tiledlayout(1,3,"TileSpacing","compact");

val = squeeze(eta(:,sliceIndex,:)).';

ax = nexttile;
% pcolor(wvt.x/1000,wvt.z,squeeze(eta(:,sliceIndex,:)).'); clim([min(eta(:)),max(eta(:))]), colorbar; shading interp, hold on
pcolor(wvt.x/1000,wvt.z,val), shading interp, hold on
clim([-1 1]*max(abs(val(:))))
cb = colorbar;
cb.Label.String = 'm';
cb.Location = "south";
contour(wvt.x/1000,wvt.z,squeeze(rho_total(:,sliceIndex,:)).',linspace(min(rho_total(:)),max(rho_total(:)),10),'w','LineWidth',1.5);
contour(wvt.x/1000,wvt.z,squeeze(rho_total(:,sliceIndex,:)).',linspace(min(rho_total(:)),max(rho_total(:)),10),'k','LineWidth',1.0);
colormap(ax,mycolormap);
title('\eta')
xlabel('distance (km)')
ylabel('depth (m)')
xlim(xlimits)


val2 = squeeze(eta_e(:,sliceIndex,:)).';

ax = nexttile;
% pcolor(wvt.x/1000,wvt.z,squeeze(eta(:,sliceIndex,:)).'); clim([min(eta(:)),max(eta(:))]), colorbar; shading interp, hold on
pcolor(wvt.x/1000,wvt.z,val2), shading interp, hold on
clim([-1 1]*max(abs(val(:))))
cb = colorbar;
cb.Label.String = 'm';
cb.Location = "south";
contour(wvt.x/1000,wvt.z,squeeze(rho_total(:,sliceIndex,:)).',linspace(min(rho_total(:)),max(rho_total(:)),10),'w','LineWidth',1.5);
contour(wvt.x/1000,wvt.z,squeeze(rho_total(:,sliceIndex,:)).',linspace(min(rho_total(:)),max(rho_total(:)),10),'k','LineWidth',1.0);
colormap(ax,mycolormap);
title('\eta_e')
xlabel('distance (km)')
ylabel('depth (m)')
xlim(xlimits)


eta_error_pct = eta_e./eta;
eta_error_pct(abs(eta)<1) = 1;
eta_error_pct = 100*(eta_error_pct-1);
val = squeeze(eta_error_pct(:,sliceIndex,:)).';

ax = nexttile;
pcolor(wvt.x/1000,wvt.z,val), shading interp
clim([-1 1]*max(abs(val(:))))
cb = colorbar;
cb.Label.String = '%';
cb.Location = "south";
colormap(ax,(mycolormap));
title('\eta_e error')
ax.YTickLabel = [];
xlabel('distance (km)')
xlim(xlimits)

title(tl,'Cyclonic eddy');

if shouldSavePlots == 1
    print(sprintf('%s-cyclonic.png',stratificationName),'-dpng','-r300');
end



%%
% Ekj = wvt.transformToRadialWavenumber(wvt.A0_TE_factor .* abs(wvt.A0).^2);
% figure, pcolor(wvt.kRadial,wvt.j,log10(Ekj).')
% figure, plot(wvt.j,sum(Ekj,1)), ylog

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Constituent parts
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

zeta_x = wvt.diffY(wvt.w) - wvt.diffZF(wvt.v); % w_y - v_z
zeta_y = wvt.diffZF(wvt.u) - wvt.diffX(wvt.w);  % u_z - w_x
zeta_z = wvt.diffX(wvt.v) - wvt.diffY(wvt.u);  % v_x - u_y

APVrv = zeta_z;
APVvs = - wvt.f * wvt.diffZG(eta);
APVnl = zeta_x .* wvt.diffX(eta) + zeta_y .* wvt.diffY(eta) + zeta_z .* wvt.diffZG(eta);
apv = wvt.apv;

pvlim = max(abs(APVvs(:)/wvt.f))*[-1 1];

figure(name='Potential vorticity components')
tl = tiledlayout(1,3);
title(tl,'Potential vorticity components');

val = squeeze(APVrv(:,sliceIndex,:)/wvt.f).';

ax = nexttile;
pcolor(wvt.x/1000,wvt.z,val); clim(0.1*pvlim), shading interp
colormap(ax,(mycolormap));
cb = colorbar;
cb.Label.String = 'f_0';
cb.Location = "south";
title('\zeta_z')
xlabel('distance (km)')
ylabel('depth (m)')
xlim(xlimits)

val = squeeze(APVvs(:,sliceIndex,:)/wvt.f).';

ax = nexttile;
pcolor(wvt.x/1000,wvt.z,val); clim(pvlim), shading interp
colormap(ax,(mycolormap));
cb = colorbar;
cb.Label.String = 'f_0';
cb.Location = "south";
title('-f_0 \partial_z \eta')
xlabel('distance (km)')
ylabel('depth (m)')
xlim(xlimits)

val = squeeze(APVnl(:,sliceIndex,:)/wvt.f).';

ax = nexttile;
pcolor(wvt.x/1000,wvt.z,val); clim(0.1*pvlim), shading interp
colormap(ax,(mycolormap));
cb = colorbar;
cb.Label.String = 'f_0';
cb.Location = "south";
title('non-linear')
xlabel('distance (km)')
ylabel('depth (m)')
xlim(xlimits)

if shouldSavePlots == 1
    print(sprintf('%s-pv-parts.png',stratificationName),'-dpng','-r300');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Compute PV
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

U = -0.30; % m/s
psi = @(x,y,z) U*(Le/sqrt(2))*exp(1/2)*exp(-((x-x0)/Le).^2 -((y-y0)/Le).^2 -(z/He/sqrt(2)).^2 );
wvt.setGeostrophicStreamfunction(psi);
if shouldApplyFilter == 1
    [Qk,Ql,Qj] = ExponentialFilter(wvt,nDampedModes=20,shouldAntiAlias=1);
    wvt.A0 = Qj.*wvt.A0;
end

wvt.removeEnergyFromAliasedModes;

qgpv = wvt.qgpv;
apv = wvt.apv;
pvdiff = apv-qgpv;

m = cmocean('balance');
figure(name="Potential vorticity")

tl = tiledlayout(1,3);
ax = nexttile;
pcolor(wvt.x/1000,wvt.z,squeeze(apv(:,sliceIndex,:)/wvt.f).'); clim(max(abs(apv(:)/wvt.f))*[-1 1]), shading interp
colormap(ax,(cmocean('balance')));
cb = colorbar;
cb.Label.String = 'f_0';
cb.Location = "south";
title('apv')
xlabel('distance (km)')
ylabel('depth (m)')
xlim(xlimits)

ax = nexttile;
pcolor(wvt.x/1000,wvt.z,squeeze(qgpv(:,sliceIndex,:)/wvt.f).'); clim(max(abs(apv(:)/wvt.f))*[-1 1]), shading interp
colormap(ax,(cmocean('balance')));
cb = colorbar;
cb.Label.String = 'f_0';
cb.Location = "south";
title('qgpv')
ax.YTickLabel = [];
xlabel('distance (km)')
ylabel('depth (m)')
xlim(xlimits)


eta_error_pct = qgpv./apv;
eta_error_pct(abs(apv/wvt.f)<0.1) = 1;
eta_error_pct = 100*(eta_error_pct-1);
val = squeeze(eta_error_pct(:,sliceIndex,:)).';

ax = nexttile;
pcolor(wvt.x/1000,wvt.z,val), shading interp
colormap(ax,(cmocean('balance')));
% clim([-1 1]*max(abs(val(:))))
clim([-1 1]*40)
cb = colorbar;
cb.Label.String = '%';
cb.Location = "south";
colormap(ax,(cmocean('balance')));
title('\eta_e error')
ax.YTickLabel = [];
xlabel('distance (km)')
title('apv-qgpv')
xlim(xlimits)

title(tl,'Potential vorticity');

if shouldSavePlots == 1
    print(sprintf('%s-pv.png',stratificationName),'-dpng','-r300');
end

return




