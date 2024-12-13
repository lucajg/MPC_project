Ts = 1/10; % Sample time
car = Car(Ts);
[xs, us] = car.steady_state(120 / 3.6);
sys = car.linearize(xs, us);
[sys_lon, sys_lat] = car.decompose(sys);
% Design MPC controller
H_lon = 12; % Horizon length in seconds
H_lat = 12;
mpc_lon = MpcControl_lon(sys_lon, Ts, H_lon);
mpc_lat = MpcControl_lat(sys_lat, Ts, H_lat);
% Get control input for longitudinal subsystem
%u_lon = mpc_lon.get_u(x_lon, ref_lon);

%%% In your test script:
x_lon = [0; 80/3.6];
ref_lon = 120/3.6;
[u_lon, X_lon, U_lon] = mpc_lon.get_u(x_lon,ref_lon);

x_lat = [0; 0];
ref_lat = 3;
[u_lat, X_lat, U_lat] = mpc_lat.get_u(x_lat,ref_lat);

t = 0:Ts:H_lon;
t_inputs = t(1:end-1);  % if U_lon has one fewer element than X_lon

figure(1);

% Plot U on the first subplot
subplot(2,1,1);
plot(t_inputs, U_lon(:), 'LineWidth', 1.5);
grid on;
xlabel('Time [s]');
ylabel('Input (u_T) [-]');

% Plot X on the second subplot (for example, the second state)
subplot(2,1,2);
plot(t, X_lon(2,:), 'LineWidth', 1.5);
grid on;
xlabel('Time [s]');
ylabel('State (V) [m/s]');

sgtitle('u_T and V over time');

t = 0:Ts:H_lat;
t_inputs = t(1:end-1);  % if U_lon has one fewer element than X_lon

figure(2);

% Plot U on the first subplot
subplot(3,1,1);
plot(t_inputs, U_lat(:), 'LineWidth', 1.5);
grid on;
xlabel('Time [s]');
ylabel('Input (\delta) [rad]');

% Plot X on the second subplot (for example, the second state)
subplot(3,1,2);
plot(t, X_lat(1,:), 'LineWidth', 1.5);
grid on;
xlabel('Time [s]');
ylabel('State (y) [m]');

% Plot X on the second subplot (for example, the second state)
subplot(3,1,3);
plot(t, X_lat(2,:), 'LineWidth', 1.5);
grid on;
xlabel('Time [s]');
ylabel('State (\theta) [rad]');

sgtitle('\delta, y and \theta over time');
% Plot and debug X lon and U lon accordingly
%%%