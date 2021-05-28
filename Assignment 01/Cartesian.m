tf = 25;
dt = 0.001;
t = 0:dt:tf;
x(1) = 1;
y(1) = 0;

for i = 1 : length(t)
    xn(i) = x(i) - dt*5*y(i);
    yn(i) = y(i) + dt*5*x(i);
    x(i+1) = x(i) - (1/2)*dt*5*(y(i)+yn(i));
    y(i+1) = y(i) + (1/2)*dt*5*(x(i)+xn(i));
end

plot(x,y)