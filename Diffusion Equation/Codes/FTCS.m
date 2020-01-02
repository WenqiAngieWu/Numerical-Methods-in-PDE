%%% Forward Time Center Space Scheme %%%
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
for n = 1:N
    for i = 2:N
        u(i, n+1) = u(i, n) + 0.5*f*(u(i+1, n)-2*u(i, n)+u(i-1, n));
    end
end

t = t(100);
exact = exp(-0.5.*t).*sin(x); % exact solution
key = u(:,100);
error = key - exact;

subplot(1,2,1),plot(x, key)
xlabel('x')
ylabel('u (when tn = 100)')


hold on 
subplot(1,2,1),plot(x, error, '--')
legend('Numerical result','Error')
title('FTCS scheme')
grid on

subplot(1,2,2), plot(x, error)
xlabel('x')
ylabel('u (when tn = 100)')
legend('Error')
title('FTCS scheme')
grid on

[X, T] = meshgrid(x, t);
exact = exp(-0.5.*T).*sin(X); % exact solution
exact = exact';
mesh(X, T, exact-u);
xlabel('x')
ylabel('t')
zlabel('error of u')
title('FTCS scheme') 