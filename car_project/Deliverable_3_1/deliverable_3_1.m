Ts = 1/10; % Sample time
car = Car(Ts);
[xs, us] = car.steady_state(120 / 3.6);
sys = car.linearize(xs, us);
[sys_lon, sys_lat] = car.decompose(sys);
% Design MPC controller
H_lon = 30; % Horizon length in seconds
mpc_lon = MpcControl_lon(sys_lon, Ts, H_lon);
% Get control input for longitudinal subsystem
%u_lon = mpc_lon.get_u(x_lon, ref_lon);

%%% In your test script:
x_lon = [0; 80/3.6] ;
ref_lon = 120/3.6;
[u_lon, X_lon, U_lon] = mpc_lon.get_u(x_lon,ref_lon);

t = 0:Ts:H_lon;
t_inputs = t(1:end-1);  % if U_lon has one fewer element than X_lon

figure(1); clf;

% Plot U on the first subplot
subplot(2,1,1);
plot(t_inputs, U_lon(:), 'LineWidth', 1.5);
grid on;
xlabel('Time [s]');
ylabel('Input (U) [-]');

% Plot X on the second subplot (for example, the second state)
subplot(2,1,2);
plot(t, X_lon(2,:), 'LineWidth', 1.5);
grid on;
xlabel('Time [s]');
ylabel('State (X_2) [m/s]');

sgtitle('U and X Over Time');
% Plot and debug X lon and U lon accordingly
%%%