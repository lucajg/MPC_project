Ts = 1/10;
car = Car(Ts);

Tf = 2.0;        % Simulation end time

x0 = [0, 0, deg2rad(-2), 20/3.6]';

u = [deg2rad(-1), 0.7]';

params  = {};     % setup simulation parameter struct

params.Tf = Tf;
params.myCar.model = car;
params.myCar.x0 = x0;
params.myCar.u = u;
result = simulate(params);    % Simulate non liner model
visualization(car, result);

result.T          % Time at every simulation step
result.myCar.X    % state trajectory
result.myCar.U    % input trajectory