A = [1 2 1; 3 8 1; 0 4 1];
b = [2; 12; 2];
x = rand(length(b),1);
tol = 1e-13;
fprintf('Gauss-Seidel\n');
x = Gauss_Siedel(A,b,x,tol)

