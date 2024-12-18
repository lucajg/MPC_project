Ts = 1/10; % Sample time
car = Car(Ts);
[xs, us] = car.steady_state(120 / 3.6);
sys = car.linearize(xs, us);
[sys_lon, sys_lat] = car.decompose(sys);
% Design MPC controller
H_lon = 10; % Horizon length in seconds
H_lat = 10;
mpc_lon = MpcControl_lon(sys_lon, Ts, H_lon);
mpc_lat = MpcControl_lat(sys_lat, Ts, H_lon);
mpc = car.merge_lin_controllers(mpc_lon, mpc_lat);

estimator = LonEstimator(sys_lon, Ts);
x0 = [0 0 0 80/3.6]'; % (x, y, theta, V)
ref1 = [0 80/3.6]'; % (y_ref, V_ref)
ref2 = [3 120/3.6]'; % (y_ref, V_ref)
params = {};
params.Tf = 15;
params.myCar.model = car;
params.myCar.x0 = x0;
params.myCar.est_fcn = @estimator.estimate;
params.myCar.est_dist0 = 0;
params.myCar.u = @mpc.get_u;
params.myCar.ref = car.ref_step(ref1, ref2, 2); % delay reference step by 2s;
result = simulate(params);
visualization(car, result);