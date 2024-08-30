% AMS figure widths, given in picas, converted to points (1 pica=12 points)
scaleFactor = 1;
LoadFigureDefaults

% Parameters for all figures
latitude = 35;
Le = 120e3;
He = 350;

f = 2 * 7.2921E-5 * sin( latitude*pi/180 );
Hz = @(z) -(z/He/He) .* exp( -z.^2/(2*He*He));

N0 = 3*2*pi/3600; % buoyancy frequency at the surface, radians/seconds
b = 1300; % thermocline exponential scale, meters
N2exp= @(z) N0*N0*exp(2*z/b);


Us = (He/Le)*(N0*N0*b)/(sqrt(2)*f);
ssh = f*Us*Le/(sqrt(2)*9.81);
fprintf('Us= %.1f cm/s, ssh=%.1f cm\n',100*Us,100*ssh);


D = 4000;
z = linspace(-D,0,100).';
rhoe = @(z) -Hz(z)/Hz(-He);

lb = @(z) (exp(2*z/b) - 1)*Us;
rb = @(z) (exp(2*z/b) - exp(-2*D/b))*Us;

% anti-cyclone/cyclone amplitude
% For 500 m
% Aac = .93*Us*exp(-1/2)*He*2/b;
% Ac = 0.85*Us*(exp(-2*He/b) - exp(-2*D/b));

% For 350 m
Aac = 0.98*Us*exp(-1/2)*He*2/b;
Ac = 0.92*Us*(exp(-2*He/b) - exp(-2*D/b));

% For 800 m
% Aac = 0.84*Us*exp(-1/2)*He*2/b;
% Ac = 0.62*Us*(exp(-2*He/b) - exp(-2*D/b));

figure(position=[50 50 figure_width_1col+8 300*scaleFactor],name=sprintf('Eddy bounds: L_e=%d km, H_e=%d m',round(Le/1e3),round(He)),PaperUnits="points")
set(gcf, 'Color', 'w');
plot(-Ac*rhoe(z)*100,z,LineWidth=2), hold on
plot(Aac*rhoe(z)*100,z,LineWidth=2)
plot(lb(z)*100,z,Color=0*[1 1 1],LineWidth=2)
plot(rb(z)*100,z,Color=0*[1 1 1],LineWidth=2)
plot([0 0],[-2000 0],LineWidth=1,Color=[0 0 0])
plot([-Ac*rhoe(-He) Aac*rhoe(-He)]*100, [-He -He],LineWidth=1,Color=[0 0 0])

ylim([-2000 0])
ax.YTickLabel = [];
xlabel('U (cm/s)','FontSize', figure_axis_label_size, 'FontName', figure_font)
set( gca, 'FontSize', figure_axis_tick_size);

print('eddy-bounds.eps','-depsc2');

% title(sprintf('L_e=%d km, H_e=%d m',round(Le/1e3),round(He)))