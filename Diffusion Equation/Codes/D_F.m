%%% Dufort-Frankel %%%
% @para:
% input: order, which means the order of  the truncation error
function [] = D_F(order)

N =500;
u = zeros(N+1, N+1);
h = 2*pi/N; % spatial step
ep = order/2 + 1;
time_step = h^ep; % time step
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


% second layer of time
for i = 2:N
    u(i ,2) = u(i, 1)+0.5*f*(u(i+1, 1)-2*u(i, 1)+u(i-1, 1));
end

% difference scheme
for n = 2:N
    for i = 2:N
        u(i, n+1) = f/(1+f)*u(i+1, n) + f/(1+f)*u(i-1, n) +(1-f)/(1+f)*u(i, n-1);
    end
end

t = t(100);
exact = exp(-0.5.*t).*sin(x);
key = u(:,100);
error = key - exact;
% plot(x, key)
% hold on 
plot(x, error)
legend('Error')
%legend('Numerical result','Error')
grid on

[X, T] = meshgrid(x, t);
exact = exp(-0.5.*T).*sin(X); % exact solution
exact = exact';
mesh(X, T, exact-u);
xlabel('x')
ylabel('t')
zlabel('error of u')
title('D-F scheme') 

end


