title('Trajectory of Particle in Vortex')
theta = 0:0.5:2*pi;
rho = ones(length(theta));
polarplot(theta,rho,'b--')
rlim([0 1.1])

hold on;
theta = 0:0.05:2*pi;
rho = ones(length(theta));
polarplot(theta,rho,'r-')

hold on;
theta = 0:0.005:2*pi;
rho = ones(length(theta));
polarplot(theta,rho,'g-')
%polarplot(theta,rho,'g-','LineWidth',0.5)






