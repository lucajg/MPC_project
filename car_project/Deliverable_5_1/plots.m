load('tube_mpc_data.mat', 'Xf', 'x_safe_pos', 'K', 'U', 'U_tilde', 'X', ...
     'X_tilde', 'E', 'Q', 'R', 'W')
close all 

%% Plot the polyhedron W
figure; 
W.plot('color','r','alpha',0.1)
xlabel('x');
ylabel('V');
title('Initial Polyhedron W');

%% Plot the polyhedron E
figure;        
hold on;
E.plot('color','b','alpha',0.5);
xlabel('$x$','Interpreter', 'latex');
ylabel('$V$','Interpreter', 'latex');
title('Minimal robust positive invariant set $\mathcal{E}$','Interpreter', 'latex');
plot(0, 0, 'blackx', 'MarkerSize', 7, 'LineWidth', 1); 
axis square;
hold off;

%% Plot the polyhedrons X and X_tilde
figure;
tiledlayout(2,1)

% X
nexttile
X.plot('color','blue','alpha',0.1);
xlabel('x');
ylabel('V');
title('X')

% X_tilde
nexttile
X_tilde.plot('color','blue','alpha',0.1);
xlabel('x');
ylabel('V');
title('X tilde')

%% Plot the polyhedrons U and U_tilde
figure;
tiledlayout(2,1);

% U
nexttile
U.plot('color','g','alpha',0.1);
xlabel('u_T');
title('U')

% U_tilde
nexttile
U_tilde.plot('color','g','alpha',0.1);
xlabel('u_T');
title('U tilde');

%% Plot K*E
KE = K*E;
figure;
KE.plot('color','blue','alpha',0.1);
xlabel('u_T');
title('K*E')

%% Plot Maximal invaritant set Xf

figure;
hold on;
Xf.plot('color','r','alpha',0.3);
title('Maximal robust positive invariant set $\mathcal{X}_f$','Interpreter', 'latex');
xlabel('$x$','Interpreter', 'latex');
ylabel('$V$','Interpreter', 'latex');
grid on;
plot(0, 0, 'blackx', 'MarkerSize', 7, 'LineWidth', 1); 
axis square;
hold off;