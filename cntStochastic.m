clear all, close all
T = 1; N = 20; dt = T/N; x = zeros(1,N); x_0=0;
mu=82; Sigma=4.8;
u=rand(N,1);
z=norminv(u,mu,Sigma);
x(1)=x_0+z(1);
for j = 2:N
   x(j) = x(j-1) + z(j); % next
end
x
scatter([0:dt:T],[0,x]) % W(0) = 0
xlabel(’t’, ’FontSize’, 12), ylabel(’X(t)’, ’FontSize’, 12)
