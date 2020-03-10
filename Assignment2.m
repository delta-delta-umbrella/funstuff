% Script for Assignment #2, MECH 5304, CFD
% Drafted by: Dan Duong
% Carleton student number: 100960337, uO student number: 7533396  
% Draft date: Feb. 26th, 2020
% Last date modified: Mar. 10th, 2020
% Assignment due: Mar. 18th, 2020

% Progress to date:
% -Relatively complete.

% Assignment Description:
% Stokes flow in 1-D, using FTCS method, in explicit time,
% solve and plot results for different times

% Initial given parameters, all variables in standard SI units as indicated
h = 0.01; %[m]
Re = [0.1, 1, 10];
nu = 1*10^-6; % [m^2/s] kinematic viscosity
omega = Re.*nu./h^2; % [s^-1] the frequency from the Re number

%%

% Analysis for Re(1)
ny = 70; % number of spatial steps 
dy = h/ny; % [m] incremental step size
double_T = 2*2*pi/omega(1); % [s] total time, calculated for 2 periods 
T = double_T/2; % Period of the oscillation
nts = 1250000; % number of time steps
dt = double_T/nts; % [s] incremental time size

limit = 2*nu*dt/dy^2; % CFL criterion set out in Assignment

for i = 1:ny+1
    y(i) = (i-1)*dy; % Position
    u(i,1) = 0; % Inital spatial conditions at t=0
end

y = flip(y); % Flip for top plate stationary

for k=1:nts+1
    t(k) = (k-1)*dt; % Time
    u(1,k) = 0; % Initial boundary time conditions at x = top plate
    u(ny+1,k) = sin(omega(1)*t(k)); % Initial boundary time condition for moving plate
end

% Utilize FTCS scheme to solve for u(i,k+1)
for k=1:nts
    for i=2:ny
        u(i,k+1) = u(i,k) + 0.5*limit*(u(i-1,k)+u(i+1,k)-2.*u(i,k));
    end
end

% Plotting

figure(1)
hold on
box on
grid on
axis([-1.1,1.1,0,1])
plot(u(:,1), y/h, 'k')
plot(u(:,52063) ,y/h,'k--')
plot(u(:,156205), y/h,'k-.')
plot(u(:,364479), y/h,'k:')
plot(u(:,468613), y/h,'k.')

xlabel('$$\it{u}/{\it{U_{max}}}$$','FontSize',14,'FontName','Times New Roman','Interpreter','Latex');
ylabel('$$\it{y}/{\it{h}}$$','FontSize',14,'FontName','Times New Roman','Interpreter','Latex');
legend('t = 0T [s]', 't = 0.0833T [s]', 't = 0.25T [s]', 't = 0.5832T [s]', 't = 0.75T [s]')

clearvars u y t
%%

% Analysis for Re(2)
ny = 100; % number of spatial steps 
dy = h./ny; % [m] incremental step size
double_T = 2*2*pi./omega(2); % [s] total time, calculated for 2 cycles
T = double_T./2; % Period of the oscillation
nts = 300000; % number of time steps
dt = double_T./nts; % [s] incremental time size

limit = 2*nu*dt./dy^2; % CFL criterion set out in Assignment

for i = 1:ny+1
    y(i) = (i-1)*dy; % Position
    u(i,1) = 0; % Inital spatial conditions at t=0
end

y = flip(y); % Flip for top plate stationary

for k=1:nts+1
    t(k) = (k-1)*dt; % Time
    u(1,k) = 0; % Initial boundary time conditions at x = top plate
    u(ny+1,k) = sin(omega(2)*t(k)); % Initial boundary time condition for moving plate
end

% Utilize FTCS scheme to solve for u(i,k+1)
for k=1:nts
    for i=2:ny
        u(i,k+1) = u(i,k) + 0.5*limit*(u(i-1,k)+u(i+1,k)-2.*u(i,k));
    end
end

% Plotting

figure(2)
hold on
box on
grid on
axis([-1.1,1.1,0,1])
plot(u(:,1), y/h, 'k')
plot(u(:,12495) ,y/h,'k--')
plot(u(:,37489), y/h,'k-.')
plot(u(:,87478), y/h,'k:')
plot(u(:,112467), y/h,'k.')

xlabel('$$\it{u}/{\it{U_{max}}}$$','FontSize',14,'FontName','Times New Roman','Interpreter','Latex');
ylabel('$$\it{y}/{\it{h}}$$','FontSize',14,'FontName','Times New Roman','Interpreter','Latex');
legend('t = 0T [s]', 't = 0.0833T [s]', 't = 0.25T [s]', 't = 0.5832T [s]', 't = 0.75T [s]')

% Spatial Grid analysis plots

figure(9)
hold on
box on
grid on
axis([-1.1,1.1,0,1])

plot(u(:,37489), y/h,'k--')
plot(u(:,87478), y/h,'k--')

xlabel('$$\it{u}/{\it{U_{max}}}$$','FontSize',14,'FontName','Times New Roman','Interpreter','Latex');
ylabel('$$\it{y}/{\it{h}}$$','FontSize',14,'FontName','Times New Roman','Interpreter','Latex');

% Legend to be added later

% Temporal Grid analysis plots

figure(11)
hold on
box on
grid on
axis([-1.1,1.1,0,1])

plot(u(:,37489), y/h,'k')
plot(u(:,87478), y/h,'k')

xlabel('$$\it{u}/{\it{U_{max}}}$$','FontSize',14,'FontName','Times New Roman','Interpreter','Latex');
ylabel('$$\it{y}/{\it{h}}$$','FontSize',14,'FontName','Times New Roman','Interpreter','Latex');

%Legend to be added after all data plotted

clearvars u y t
%%

% Analysis for Re(3)
ny = 100; % number of spatial steps 
dy = h/ny; % [m] incremental step size
double_T = 2*2*pi/omega(3); % [s] total time, calculated for 2 cycles 
T = double_T/2; % Period of the oscillation
nts = 50000; % number of time steps
dt = double_T/nts; % [s] incremental time size

limit = 2*nu*dt/dy^2; % CFL criterion set out in Assignment

for i = 1:ny+1
    y(i) = (i-1)*dy; % Position
    u(i,1) = 0; % Inital spatial conditions at t=0
end

y = flip(y); % Flip for top plate stationary

for k=1:nts+1
    t(k) = (k-1)*dt; % Time
    u(1,k) = 0; % Initial boundary time conditions at x = top plate
    u(ny+1,k) = sin(omega(3)*t(k)); % Initial boundary time condition for moving plate
end

% Utilize FTCS scheme to solve for u(i,k+1)
for k=1:nts
    for i=2:ny
        u(i,k+1) = u(i,k) + 0.5*limit*(u(i-1,k)+u(i+1,k)-2.*u(i,k));
    end
end

% Plotting

figure(3)
hold on
box on
grid on
axis([-1.1,1.1,0,1])
plot(u(:,1), y/h, 'k')
plot(u(:,2083) ,y/h,'k--')
plot(u(:,6248), y/h,'k-.')
plot(u(:,14579), y/h, 'k:')
plot(u(:,18745), y/h, 'k.')
plot(u(:,27083), y/h, 'r--')

xlabel('$$\it{u}/{\it{U_{max}}}$$','FontSize',14,'FontName','Times New Roman','Interpreter','Latex');
ylabel('$$\it{y}/{\it{h}}$$','FontSize',14,'FontName','Times New Roman','Interpreter','Latex');
legend('t = 0T [s]', 't = 0.0833T [s]', 't = 0.25T [s]', 't = 0.5832T [s]', 't = 0.75T [s]', 't = 1.0833T [s]')

% Surface plot, for interest

figure(4)
surf(t(1:1000:50001)./double_T,y(1:2:101)./h,u((1:2:101),(1:1000:50001)))
box on
grid on
axis([0,1,0,1,-1.1,1.1])
xlabel('$$\it{t}/{\it{2T}}$$','FontSize',14,'FontName','Times New Roman','Interpreter','Latex');
ylabel('$$\it{y}/{\it{h}}$$','FontSize',14,'FontName','Times New Roman','Interpreter','Latex');
zlabel('$$\it{u}/{\it{U_{max}}}$$','FontSize',14,'FontName','Times New Roman','Interpreter','Latex');

% Comparative Analysis to Analytical Solution plots

figure(12)
hold on
box on
grid on
axis([-1.1,1.1,0,1])
plot(u(:,27083), y/h, 'k-.')
plot(u(:,39580), y/h, 'k--')

xlabel('$$\it{u}/{\it{U_{max}}}$$','FontSize',14,'FontName','Times New Roman','Interpreter','Latex');
ylabel('$$\it{y}/{\it{h}}$$','FontSize',14,'FontName','Times New Roman','Interpreter','Latex');

clearvars u y t
%%

% Spatial grid analysis 1 for Re(2)
ny = 95; % number of spatial steps 
dy = h/ny; % [m] incremental step size
nu = 1*10^-6; % [m^2/s] kinematic viscosity
double_T = 2*2*pi/omega(2); % [s] total time, calculated for 2 cycles 
T = double_T/2; % Period of the oscillation
nts = 300000; % number of time steps
dt = double_T/nts; % [s] incremental time size

limit = 2*nu*dt/dy^2; % CFL criterion set out in Assignment

for i = 1:ny+1
    y(i) = (i-1)*dy; % Position
    u(i,1) = 0; % Inital spatial conditions at t=0
end

y = flip(y); % Flip for top plate stationary

for k=1:nts+1
    t(k) = (k-1)*dt; % Time
    u(1,k) = 0; % Initial boundary time conditions at x = top plate
    u(ny+1,k) = sin(omega(2)*t(k)); % Initial boundary time condition for moving plate
end

% Utilize FTCS scheme to solve for u(i,k+1)
for k=1:nts
    for i=2:ny
        u(i,k+1) = u(i,k) + 0.5*limit*(u(i-1,k)+u(i+1,k)-2.*u(i,k));
    end
end

figure(9)
hold on
box on
grid on
axis([-1.1,1.1,0,1])

plot(u(:,37489), y/h,'b--')
plot(u(:,87478), y/h,'b--')

xlabel('$$\it{u}/{\it{U_{max}}}$$','FontSize',14,'FontName','Times New Roman','Interpreter','Latex');
ylabel('$$\it{y}/{\it{h}}$$','FontSize',14,'FontName','Times New Roman','Interpreter','Latex');

% Legend to be added after all data is plotted

clearvars u y t
%%

% Spatial grid analysis 2 for Re(2)
ny = 105; % number of spatial steps 
dy = h/ny; % [m] incremental step size
nu = 1*10^-6; % [m^2/s] kinematic viscosity
double_T = 2*2*pi/omega(2); % [s] total time, calculated for 2 cycles 
T = double_T/2; % Period of the oscillation
nts = 300000; % number of time steps
dt = double_T/nts; % [s] incremental time size

limit = 2*nu*dt/dy^2; % CFL criterion set out in Assignment

for i = 1:ny+1
    y(i) = (i-1)*dy; % Position
    u(i,1) = 0; % Inital spatial conditions at t=0
end

y = flip(y); % Flip for top plate stationary

for k=1:nts+1
    t(k) = (k-1)*dt; % Time
    u(1,k) = 0; % Initial boundary time conditions at x = top plate
    u(ny+1,k) = sin(omega(2)*t(k)); % Initial boundary time condition for moving plate
end

% Utilize FTCS scheme to solve for u(i,k+1)
for k=1:nts
    for i=2:ny
        u(i,k+1) = u(i,k) + 0.5*limit*(u(i-1,k)+u(i+1,k)-2.*u(i,k));
    end
end

% Plotting

figure(9)
hold on
box on
grid on
axis([-1.1,1.1,0,1])

plot(u(:,37489), y/h,'r--')
plot(u(:,87478), y/h,'r--')

xlabel('$$\it{u}/{\it{U_{max}}}$$','FontSize',14,'FontName','Times New Roman','Interpreter','Latex');
ylabel('$$\it{y}/{\it{h}}$$','FontSize',14,'FontName','Times New Roman','Interpreter','Latex');


% All spatial grid analysis plotted, hence legend added
legend('dy = 1.0x10^{-4} m', '', 'dy = 1.053x10^{-4} m', '', 'dy = 0.952x10^{-4}')
p = findobj(gca,'Type','line');
set(get(get(p(5),'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
set(get(get(p(1),'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
set(get(get(p(3),'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
title('Re = 1.0','FontSize',14,'FontName','Times New Roman','Interpreter','Latex')

clearvars u y t
%%

% Temporal grid analysis 1 for Re(2)
ny = 100; % number of spatial steps 
dy = h/ny; % [m] incremental step size
nu = 1*10^-6; % [m^2/s] kinematic viscosity
double_T = 2*2*pi/omega(2); % [s] total time, calculated for 2 cycles 
T = double_T/2; % Period of the oscillation
nts = 350000; % number of time steps
dt = double_T/nts; % [s] incremental time size

limit = 2*nu*dt/dy^2; % CFL criterion set out in Assignment

for i = 1:ny+1
    y(i) = (i-1)*dy; % Position
    u(i,1) = 0; % Inital spatial conditions at t=0
end

y = flip(y); % Flip for top plate stationary

for k=1:nts+1
    t(k) = (k-1)*dt; % Time
    u(1,k) = 0; % Initial boundary time conditions at x = top plate
    u(ny+1,k) = sin(omega(2)*t(k)); % Initial boundary time condition for moving plate
end

% Utilize FTCS scheme to solve for u(i,k+1)
for k=1:nts
    for i=2:ny
        u(i,k+1) = u(i,k) + 0.5*limit*(u(i-1,k)+u(i+1,k)-2.*u(i,k));
    end
end

% Plotting

figure(11)
hold on
box on
grid on
axis([-1.1,1.1,0,1])

plot(u(:,43737), y/h,'b')
plot(u(:,102054), y/h,'b')

xlabel('$$\it{u}/{\it{U_{max}}}$$','FontSize',14,'FontName','Times New Roman','Interpreter','Latex');
ylabel('$$\it{y}/{\it{h}}$$','FontSize',14,'FontName','Times New Roman','Interpreter','Latex');

% Legend to be added after all data is plotted

clearvars u y t
%%

% Temporal grid analysis 2 for Re(2)
ny = 100; % number of spatial steps 
dy = h/ny; % [m] incremental step size
nu = 1*10^-6; % [m^2/s] kinematic viscosity
double_T = 2*2*pi/omega(2); % [s] total time, calculated for 2 cycles 
T = double_T/2; % Period of the oscillation
nts = 400000; % number of time steps
dt = double_T/nts; % [s] incremental time size

limit = 2*nu*dt/dy^2; % CFL criterion set out in Assignment

for i = 1:ny+1
    y(i) = (i-1)*dy; % Position
    u(i,1) = 0; % Inital spatial conditions at t=0
end

y = flip(y); % Flip for top plate stationary

for k=1:nts+1
    t(k) = (k-1)*dt; % Time
    u(1,k) = 0; % Initial boundary time conditions at x = top plate
    u(ny+1,k) = sin(omega(2)*t(k)); % Initial boundary time condition for moving plate
end

% Utilize FTCS scheme to solve for u(i,k+1)
for k=1:nts
    for i=2:ny
        u(i,k+1) = u(i,k) + 0.5*limit*(u(i-1,k)+u(i+1,k)-2.*u(i,k));
    end
end

% Plotting

figure(11)
hold on
box on
grid on
axis([-1.1,1.1,0,1])

plot(u(:,49986), y/h,'r')
plot(u(:,116633), y/h,'r')

xlabel('$$\it{u}/{\it{U_{max}}}$$','FontSize',14,'FontName','Times New Roman','Interpreter','Latex');
ylabel('$$\it{y}/{\it{h}}$$','FontSize',14,'FontName','Times New Roman','Interpreter','Latex');

% All temporal grid analysis plotted, hence legend added
legend('dt = 0.0042 s', '', 'dt = 0.0036 s', '', 'dt = 0.0031s')
p = findobj(gca,'Type','line');
set(get(get(p(5),'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
set(get(get(p(1),'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
set(get(get(p(3),'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
title('Re = 1.0','FontSize',14,'FontName','Times New Roman','Interpreter','Latex')

clearvars u y t

%%
% Comparative Analysis to Stokes second problem solution

% Solve Stokes second problem solution for field
% Load up the data for Re(3) to obtain y and t

for i = 1:length(t)
    for j =1:length(y)
        SSPS(j,i) = exp(-sqrt(0.5*omega(3)/nu)*y(j))*sin(omega(3)*t(i)-(sqrt(0.5*omega(3)/nu)*y(j)));
    end
end

figure(12)
hold on
box on
grid on
axis([-1.1,1.1,0,1])

plot(SSPS(:,27083), y/h, 'b')
plot(SSPS(:,39580), y/h, 'b')

legend('Numerical Solution', '', 'Stokes Second Problem Solution', '')
p = findobj(gca,'Type','line');
set(get(get(p(2),'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
set(get(get(p(4),'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
title('Re = 10.0','FontSize',14,'FontName','Times New Roman','Interpreter','Latex')
