function x = Gauss_Siedel(A,b,x0,tol)

[na, ma] = size(A);
if na ~= ma
    disp('ERROR: Matrix A must be a square matrix')
    return
end

[nb, mb] = size (b);
if nb ~= na || mb~=1
   disp( 'ERROR: Matrix b must be a column matrix')
   return
end
% Separation of matrix A into lower triangular and upper triangular matrices
% A = D + L + U
D = diag(diag(A));
L = tril(A)- D;
U = triu(A)- D;
% check for convergence condition for Gauss-Seidel method
e= max(eig(-inv(D+L)*(U)));
if abs(e) >= 1
    disp ('Since the modulus of the largest Eigen value of iterative matrix is not less than 1') 
    disp ('this process is not convergent.')
    return
end
k = 1;
x(:,1) = x0;
err = 1000000000*rand(na,1);% initial error assumption for looping
i=1;
while sum(abs(err) >= tol) ~= zeros(na,1)
    
    x(:,k+1) = -inv(D+L)*(U)*x(:,k) + inv(D+L)*b; % Gauss-Seidel formula
  
    err = x(:,k+1) - x(:, k);% finding error
    k = k + 1;
    i= i+1;
end
    for i=1:10
       fprintf('\n\nIteration %2d:', i);
       fprintf('\n%f',x(:,i));
    end
fprintf('\n\nThe final answer obtained after %g iterations is  \n', k);
x = x(:,end);
