%%% Crank-Nicholson Method %%%

clc, clear, close
N =500;
u = zeros(N+1, N+1);
h = 2*pi/N; % spatial step

time_step = h^2; % time step
f = time_step/(h^2);
x = 0:h:(2*pi);
x = x';
t = 0:time_step:(N*time_step);
t = t';

% initial data
u(:,1) = sin(x);  

% boundary
u(1,:) = 0; 
u(N+1,:) = 0; 

% difference scheme
diag1 = (2+f).*ones(N-1, 1); 
diag2 = -0.5.*f.*ones(N-2, 1);
AA = diag(diag1) + diag(diag2, 1) + diag(diag2, -1); % coefficient matrix

for n = 1:N  % loop over time
    b(1) = 0.5*f*u(1, n) + 0.5*f*u(3, n) + 0.5*f*u(1, n) + (2-f)*u(2, n);
    b(N-1) = 0.5*f*u(N+1, n) + 0.5*f*u(N-1, n) + 0.5*f*u(N+1, n) + (2-f)*u(N, n);
    for i = 2:N-2
        b(i) = 0.5*f*u(i+2, n) + 0.5*f*u(i, n) + (2-f)*u(i+1, n);
    end 
    u(2:N, n+1) = AA\(b');
end

t = t(100);
exact = exp(-0.5.*t).*sin(x); % exact solution
key = u(:,100);
error = key - exact;
subplot(1,2,1),plot(x, key)
title('C-N scheme')

hold on 
subplot(1,2,1),plot(x, error, '--')
legend('Numerical result','Error')
grid on

subplot(1,2,2), plot(x, error)
xlabel('x')
ylabel('u (when tn = 100)')
legend('Error')
title('C-N scheme')
grid on

[X, T] = meshgrid(x, t);
exact = exp(-0.5.*T).*sin(X); % exact solution
exact = exact';
mesh(X, T, exact-u);
xlabel('x')
ylabel('t')
zlabel('error of u')
title('C-N scheme') 


