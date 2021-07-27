function re = GFEM(n)

K = zeros(n, n);
for c = 1:n
    for r = 1:n
        if r == c && r == 1
            K(r, c) = n;
        elseif r == c
            K(r, c)= 2*n;
        elseif abs(r - c) == 1
            K(r, c) = -n;
        else
            K(r, c) = 0;
        end
    end
end

t = 0;
syms x
N = cell(1, n);
N{1} = triangularPulse(0, 0, 1/ n, x);

for i = 2:n
N{i} = triangularPulse(t, t+1/n, t+2/n, x);
t = t + 1/n;
end

F = zeros(n,1);

for i = 1:n
F(i) = int(N{i}.*x, x, 0,1);
end

db = K\F;

uh = 0;
for i = 1:n 
    uh = uh + db(i)*N{i};
end

u = (1 - x^3)/6; du = diff(u, x); duh = diff(uh, x);
f1 = subs(du, x, 1/(2*n)); f2 = subs(duh, x, 1/(2*n));
re = abs(f1 - f2)*2;

figure(1)
hold on
fplot(uh)
fplot(u)
hold off
grid on
ylabel('U');
xlabel('x');
xlim([0 1]);
title('Plot of Solution');
legend('Galerkin Finite Element Method', 'Exact solution' );

figure(2)
hold on
fplot(duh)
fplot(du)
hold off
grid on
ylabel('Ux');
xlabel('x');
xlim([0 1]);
title('Plot of Slope');
legend('Galerkin Finite Element Method', 'Exact solution');
end