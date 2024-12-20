
Ts = 1/10; % Sample time
car = Car(Ts);
H = 10;
mpc = NmpcControl(car, H);
x0 = [0 0 0 80/3.6]';
ref = [3 100/3.6]';
u = mpc.get_u(x0, ref) % check if the open−loop prediction is reasonable
mpc.sol.value(mpc.U)
mpc.sol.value(mpc.X)
mpc.sol.stats.success
mpc.sol.stats.return_status
RK4([0;0;0;0],[0;1],car.Ts,(@car.f))
