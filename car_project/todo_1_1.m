Ts = 1/10;
car = Car(Ts);

delta = 0.01;
u_T = 0;
xx = 0;
yy = 0;
theta = 0.1;
V = 10;

u = [delta, u_T]';
x = [xx, yy, theta, V]';
x_dot = car.f(x, u)