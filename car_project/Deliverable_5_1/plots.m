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
E.plot('color','b','alpha',0.5);
xlabel('x_1');
ylabel('x_2');
title('Final Polyhedron E');

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
Xf.plot('color','r','alpha',0.1);
title('Maximal invariant set X_f');
xlabel('x');
ylabel('V');
grid on;