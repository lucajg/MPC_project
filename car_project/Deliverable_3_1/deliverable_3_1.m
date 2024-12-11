Ts = 1/10; % Sample time
car = Car(Ts);
[xs, us] = car.steady_state(120 / 3.6);
sys = car.linearize(xs, us);
[sys_lon, sys_lat] = car.decompose(sys);
% Design MPC controller
H_lon = 5; % Horizon length in seconds
mpc_lon = MpcControl_lon(sys_lon, Ts, H_lon);
% Get control input for longitudinal subsystem
%u_lon = mpc_lon.get_u(x_lon, ref_lon);

%%% In your test script:
x_lon = [0; 0] ;
ref_lon = [0; 0];
[u_lon, X_lon, U_lon] = mpc_lon.get_u(sys_lon.UserData.idx,sys_lon.UserData.idy);
% Plot and debug X lon and U lon accordingly
%%%