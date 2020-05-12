% Script for Assignment #3, MECH 5304, CFD
% Drafted by: Dan Duong
% Draft date: Mar. 31st, 2020
% Last date modified: Apr. 3rd, 2020
% Assignment due: Apr. 8th, 2020

% Progress to date:
% Relatively complete
% Really needs to be cleaned up; got lazy and started hand-coding (¬_¬)

% ???   (?°?°??? ???   ???

% Assignment Description:
% Transient 1D Heat conduction, using FVM methods and TDMA to solve
% Method 1: Crank Nicolson
% Method 2: Fully implicit

% For the solving, the parameters that are required will be asked as
% inputs for the solving, as well as the method

% Note: TDMA solver is a separate .m function file

addpath('D:\Dan\Classes\MCG5334\A3')

T_init = 1;
T_bc = 0;
t_fin = 1;

mode = input('Mode? 1 = Calculation, 2 = Plotting, 3 = Optimization ');

%%
% Calculations

if mode == 1
    
cv = input('Number of control volumes? ');
    if isnumeric(cv) == 1
    elseif isnumeric(cv) == 0
        error('Error: Number of CVs entered not numeric.');
    end

timestep = input('Time step? ');
    if isnumeric(timestep) == 1
	elseif isnumeric(timestep) == 0
        error('Error: Timestep size entered not numeric.');
    end

method = input('Solving method? ''1'' = Crank-Nicholson, ''2'' = Fully Implicit ');

% Indexing parameter is to save the data into separate cell structures to
% be called upon later: For the problem described in the assignment:
% Below are indices for the Crank-Nicolson temperature cell structure
% a) CV = 21, timestep = 0.0005, index = 1
% a) CV = 21, timestep = 0.0010, index = 2
% a) CV = 21, timestep = 0.0020, index = 3
% b) timestep = 0.0010, CV = 11, index = 4
% b) timestep = 0.0010, CV = 21, index = 2
% b) timestep = 0.0010, CV = 101, index = 5
% c) timestep = 0.0010, CV = 21, index = 2
if method == 1
index = input('CN - Index number? ');
end
%%%%%%%%%%%%%%%%%%%%%%%%%% Initialize the field %%%%%%%%%%%%%%%%%%%%%%%%%%

% Initial spatial conditions
for i = 1:cv
    x(i) = (i-1)/(cv-1); % Position
    T(i,1) = T_init; % Initial thermal conditions at t=0
end

maxn = t_fin/timestep+1;
for n = 1:maxn
    t(n) = (n-1)*timestep; % Time
    T(1,n) = 0; % Initial boundary time conditions at x = left side of wall
    T(cv,n) = 0; % Initial boundary time conditions at x = right side of wall
end

%%%%%%%%%%%%%%%%%%%%%% Implement Solving Schemes%% %%%%%%%%%%%%%%%%%%%%%%%

%Crank Nicolson method
if method == 1

% R parameter, or the stability parameter:
R = timestep/(2*(1/cv)^2);    
    
    for n = 1:maxn
        for i = 2:cv-1
        f(i,n) = (1-2*R).*T(i,n) + R.*T(i-1,n) + R.*T(i+1,n);
        end
        f(i+1,n) = 0;
        
        a = (1+2*R)*ones(cv, 1);
        b = -R*ones(cv,1);
        c = -R*ones(cv,1);
        T(:,n+1) = TDMA(a, b, c, f(:,n));
    end   

T_CN{index} = T;
x_CN{index} = x;
t_CN{index} = t;
end

% Fully Implicit method
if method == 2

% R parameter, or the stablity parameter:
R = (1/cv)^2/(timestep*2); 

    for n = 1:maxn
        for i = 2:cv-1
        f(i,n) = (R).*T(i,n);
        end
        f(i+1,n) = 0;
        
        a = (1+R)*ones(cv, 1);
        b = -0.5*ones(cv,1);
        c = -0.5*ones(cv,1);
        T(:,n+1) = TDMA(a, b, c, f(:,n));
    end   

T_FI = T;    
x_FI = x;
t_FI = t;    
end

clearvars T x t f a b c

% Analytical solution field for set x, t

x_an = [0:1/200:1];
t_an = [0.008, 0.2, 0.4, 1.0];

for q = 1:length(t_an)
    for p = 1:length(x_an)
        for n = 1:10
            T_an_xi(n) = (4./(3.1415962).*(1./(2.*n-1))).*sin((2.*n-1)*(3.1415962.*x_an(p)))*exp(-(2.*n-1).^2*(3.1415962).^2.*t_an(q));
        end
        T_an(p,q) = sum(T_an_xi);
    end
end

end

% Save data command
% save('D:\Dan\Classes\MCG5334\A3\A3Data.mat', '-v7.3')

%%
% Plotting

% The following is the indices for plotting:
% For timestep = 0.0005 s
% 17 = 0.008 s
% 401 = 0.2 s
% 801 = 0.4 s
% 2001 = 1.0 s

% For timestep = 0.001 s
% 9 = 0.008 s
% 201 = 0.2 s
% 401 = 0.4 s
% 1001 = 1.0 s

% For timestep = 0.002 s
% 5 = 0.008 s
% 101 = 0.2 s
% 201 = 0.4 s
% 501 = 1.0 s

if mode == 2

% Part a) 21 CVs and C-N Method, timesteps: 0.0005, 0.001, 0.002
% Plot t = 0.008, 0.2, 0.4, 1.0
figure(1)
hold on
box on
grid on
axis([-0.2,1.2,-0.2,1.2])

plot(x_CN{1}',T_CN{1}(:,17),'ko')
plot(x_CN{2}',T_CN{2}(:,9),'k<')
plot(x_CN{3}',T_CN{3}(:,5),'kx')
plot(x_an, T_an(:,1),'k')

plot(x_CN{1}',T_CN{1}(:,401),'bo')
plot(x_CN{2}',T_CN{2}(:,201),'b<')
plot(x_CN{3}',T_CN{3}(:,101),'bx')
plot(x_an, T_an(:,2),'b')

plot(x_CN{1}',T_CN{1}(:,801),'ro')
plot(x_CN{2}',T_CN{2}(:,401),'r<')
plot(x_CN{3}',T_CN{3}(:,201),'rx')
plot(x_an, T_an(:,3),'r')

plot(x_CN{1}',T_CN{1}(:,2001),'mo')
plot(x_CN{2}',T_CN{2}(:,1001),'m<')
plot(x_CN{3}',T_CN{3}(:,501),'mx')
plot(x_an, T_an(:,4),'m')

xlabel('$$\it{x''}$$','FontSize',14,'FontName','Times New Roman','Interpreter','Latex');
ylabel('$$\it{T''}$$','FontSize',14,'FontName','Times New Roman','Interpreter','Latex');
legend('$$CV = 21, \Delta t'' = 0.0005$$', '$$CV = 21, \Delta t'' = 0.001$$', '$$CV = 21, \Delta t'' = 0.002$$', 'Analytical Solution');
title('Part a) Crank-Nicolson Method, Timestep Comparison','FontSize',14,'FontName','Times New Roman','Interpreter','Latex')

% Part b) Timestep = 0.001 and C-N Method, no. CV's: 11, 21, 101
% Plot t = 0.008, 0.2, 0.4, 1.0
figure(2)
hold on
box on
grid on
axis([-0.2,1.2,-0.2,1.2])

plot(x_CN{4}',T_CN{4}(:,9),'ko')
plot(x_CN{2}',T_CN{2}(:,9),'k<')
plot(x_CN{5}',T_CN{5}(:,9),'kx')
plot(x_an, T_an(:,1),'k')

plot(x_CN{4}',T_CN{4}(:,201),'bo')
plot(x_CN{2}',T_CN{2}(:,201),'b<')
plot(x_CN{5}',T_CN{5}(:,201),'bx')
plot(x_an, T_an(:,2),'b')

plot(x_CN{4}',T_CN{4}(:,401),'ro')
plot(x_CN{2}',T_CN{2}(:,401),'r<')
plot(x_CN{5}',T_CN{5}(:,401),'rx')
plot(x_an, T_an(:,3),'r')

plot(x_CN{4}',T_CN{4}(:,1001),'mo')
plot(x_CN{2}',T_CN{2}(:,1001),'m<')
plot(x_CN{5}',T_CN{5}(:,1001),'mx')
plot(x_an, T_an(:,4),'m')

xlabel('$$\it{x''}$$','FontSize',14,'FontName','Times New Roman','Interpreter','Latex');
ylabel('$$\it{T''}$$','FontSize',14,'FontName','Times New Roman','Interpreter','Latex');
legend('$$\Delta t'' = 0.001, CV = 11$$', '$$\Delta t'' = 0.001, CV = 21$$', '$$\Delta t'' = 0.001, CV = 101$$', 'Analytical Solution');
title('Part b) Crank-Nicolson Method, Number of CV comparison','FontSize',14,'FontName','Times New Roman','Interpreter','Latex')

% Part c) Timestep = 0.001 and no. CV's = 21, compare: C-N Method, Implicit Method 
% Plot t = 0.008, 0.2, 0.4, 1.0
figure(3)
hold on
box on
grid on
axis([-0.2,1.2,-0.2,1.2])

plot(x_CN{2}',T_CN{2}(:,9),'k<')
plot(x_FI',T_FI(:,9),'k*')
plot(x_an, T_an(:,1),'k')

plot(x_CN{2}',T_CN{2}(:,201),'b<')
plot(x_FI',T_FI(:,201),'b*')
plot(x_an, T_an(:,2),'b')

plot(x_CN{2}',T_CN{2}(:,401),'r<')
plot(x_FI',T_FI(:,401),'r*')
plot(x_an, T_an(:,3),'r')

plot(x_CN{2}',T_CN{2}(:,1001),'m<')
plot(x_FI',T_FI(:,1001),'m*')
plot(x_an, T_an(:,4),'m')

xlabel('$$\it{x''}$$','FontSize',14,'FontName','Times New Roman','Interpreter','Latex');
ylabel('$$\it{T''}$$','FontSize',14,'FontName','Times New Roman','Interpreter','Latex');
legend('$$\Delta t'' = 0.001, CV = 21, Crank-Nicolson Method$$', '$$\Delta t'' = 0.001, CV = 21, Implicit Method$$', 'Analytical Solution');
title('Part c) Crank-Nicolson Method vs Implicit Method','FontSize',14,'FontName','Times New Roman','Interpreter','Latex')

end

%%

% Optimization mode to determined the least squared error as well as
% optimize the timestep size and number of control volumes
if mode == 3

% Timestep optimization for part a)

for i = 1:length(x_CN{1})
    x_lse_t(1,i) = find(round(x_an,5) == round(x_CN{1}(1,i),5))
end

lse_t{1}(1) = sum((T_CN{1}(:,17) - T_an(x_lse_t,1)).^2)/(length(x_lse_t)-1);
lse_t{1}(2) = sum((T_CN{2}(:,9) - T_an(x_lse_t,1)).^2)/(length(x_lse_t)-1);
lse_t{1}(3) = sum((T_CN{3}(:,5) - T_an(x_lse_t,1)).^2)/(length(x_lse_t)-1);

lse_t{2}(1) = sum((T_CN{1}(:,401) - T_an(x_lse_t,2)).^2)/(length(x_lse_t)-1);
lse_t{2}(2) = sum((T_CN{2}(:,201) - T_an(x_lse_t,2)).^2)/(length(x_lse_t)-1);
lse_t{2}(3) = sum((T_CN{3}(:,101) - T_an(x_lse_t,2)).^2)/(length(x_lse_t)-1);

lse_t{3}(1) = sum((T_CN{1}(:,801) - T_an(x_lse_t,3)).^2)/(length(x_lse_t)-1);
lse_t{3}(2) = sum((T_CN{2}(:,401) - T_an(x_lse_t,3)).^2)/(length(x_lse_t)-1);
lse_t{3}(3) = sum((T_CN{3}(:,201) - T_an(x_lse_t,3)).^2)/(length(x_lse_t)-1);

lse_t{4}(1) = sum((T_CN{1}(:,2001) - T_an(x_lse_t,4)).^2)/(length(x_lse_t)-1);
lse_t{4}(2) = sum((T_CN{2}(:,1001) - T_an(x_lse_t,4)).^2)/(length(x_lse_t)-1);
lse_t{4}(3) = sum((T_CN{3}(:,501) - T_an(x_lse_t,4)).^2)/(length(x_lse_t)-1);

x_ts = [0.0005, 0.001, 0.002];

% No. CV optimization for part b)
p_array = [4,2,5];
for p = 1:3
    for i = 1:length(x_CN{p_array(p)})
        x_lse_CV{p}(1,i) = find(round(x_an,5) == round(x_CN{p_array(p)}(1,i),5));
    end
end

for p = 1:3
    i_array = [9, 202, 401, 1001];
    for i = 1:4
    lse_CV{i}(p) = sum((T_CN{p_array(p)}(:,i_array(i)) - T_an(x_lse_CV{p},i)).^2)/(length(x_lse_t)-1);
    end
end
end
x_nCV = [11, 21, 101];

% Method optimization for part c)

lse_m{1}(1) = sum((T_CN{2}(:,9) - T_an(x_lse_CV{2},1)).^2)/(length(x_lse_t)-1); 
lse_m{1}(2) = sum((T_CN{2}(:,201) - T_an(x_lse_CV{2},2)).^2)/(length(x_lse_t)-1); 
lse_m{1}(3) = sum((T_CN{2}(:,401) - T_an(x_lse_CV{2},3)).^2)/(length(x_lse_t)-1); 
lse_m{1}(4) = sum((T_CN{2}(:,1001) - T_an(x_lse_CV{2},4)).^2)/(length(x_lse_t)-1); 

lse_m{2}(1) = sum((T_FI(:,9) - T_an(x_lse_CV{2},1)).^2)/(length(x_lse_t)-1); 
lse_m{2}(2) = sum((T_FI(:,201) - T_an(x_lse_CV{2},2)).^2)/(length(x_lse_t)-1); 
lse_m{2}(3) = sum((T_FI(:,401) - T_an(x_lse_CV{2},3)).^2)/(length(x_lse_t)-1); 
lse_m{2}(4) = sum((T_FI(:,1001) - T_an(x_lse_CV{2},4)).^2)/(length(x_lse_t)-1); 

x_m = [0.008, 0.2, 0.4, 1.0];


figure(10)
hold on
box on
grid on
axis([0,0.003,10^-10,10^-2])
semilogy(x_ts, lse_t{1}, 'ko-.')
semilogy(x_ts, lse_t{2}, 'bo-.')
semilogy(x_ts, lse_t{3}, 'ro-.')
semilogy(x_ts, lse_t{4}, 'mo-.')
xlabel('$$\it{\Delta t''}$$','FontSize',14,'FontName','Times New Roman','Interpreter','Latex');
ylabel('$$\it{\sigma^2}$$','FontSize',14,'FontName','Times New Roman','Interpreter','Latex');
title('Least-squared error for $$\Delta t''$$','FontSize',14,'FontName','Times New Roman','Interpreter','Latex')

figure(11)
hold on
box on
grid on
axis([0,140,10^-10,10^-2])
semilogy(x_nCV, lse_CV{1}, 'ko-.')
semilogy(x_nCV, lse_CV{2}, 'bo-.')
semilogy(x_nCV, lse_CV{3}, 'ro-.')
semilogy(x_nCV, lse_CV{4}, 'mo-.')
ylabel('$$\it{\sigma^2}$$','FontSize',14,'FontName','Times New Roman','Interpreter','Latex');
xlabel('Number of CVs','FontSize',14,'FontName','Times New Roman','Interpreter','Latex');
title('Least-squared error for Number of CVs','FontSize',14,'FontName','Times New Roman','Interpreter','Latex')

figure(12)
hold on
box on
grid on
axis([0,1.2,10^-10,10^-2])
semilogy(x_m, lse_m{1}, 'ko-.')
semilogy(x_m, lse_m{2}, 'bo-.')
ylabel('$$\it{\sigma^2}$$','FontSize',14,'FontName','Times New Roman','Interpreter','Latex');
xlabel('$$\it{t''}$$','FontSize',14,'FontName','Times New Roman','Interpreter','Latex');
title('Least-squared error for Methods','FontSize',14,'FontName','Times New Roman','Interpreter','Latex')

