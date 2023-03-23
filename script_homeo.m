clear all;
t0 = 0; tmax = 120; %days
N = 120;
beta = 0.4;
p_0 = 2000 % Pa ~15 mmHg
%P_int = @(x) -28.139*x +2301; %Pa, x-days
P_int = @(x) p_0*(1-exp(-beta*x)); %Pa, x-days
%Diff_axial = @(x) -3*10^(-8)*x^3+2*10^(-6)*x^2-6*10^(-5)*x+7*10^(-4); %m, x-days
%a_fit = @(x) 6.886*10^(-3) + 0.159*10^(-3)*x - 0.00141*10^(-3)*x^2 + 0.573*10^(-8)*x^3;

b_theta = 10^(-7);
a_0 = 6.71/2*10^(-3); %m
h_0 = 0.12*10^(-3); %m
%c_1 = 2*10^5;
mu = 0.8*10^6 % Pa
Jm = 0.18;
sigma_0 = 20*10^3 %Pa
dt = (tmax-t0)/N;
%k=8;
k = N+1; 
tt = zeros(k,1); 
yy = zeros(k,1); 
tt(1) = t0; 
gamma_theta(1) = 1; %gamma_theta(a0, t0)
alpha(1) = 1;
a(1) = a_0; 

function y = PEq(x, P_int, gamma_theta, a0, mu, Jm, h0)
  y = zeros (1, 1);
  dWda = 2*mu*Jm*(x(1)-x(1)^(-5))/(Jm+3-x(1)^(-4)-2*x(1)^2); % Gent, x(1) = alpha
  y(1) = P_int*gamma_theta*a0*x(1)^2-dWda*h0; %equation dlya alpha = alpha(R=a0,t)
endfunction

for i = 2:k 
  tt(i) = tt(i-1) + dt; 
  m = dt*b_theta*(P_int(tt(i-1))*a(i-1)*alpha(i-1)^2/2/h_0-sigma_0*P_int(tt(i-1))/p_0); % RHS
  gamma_theta(i) = gamma_theta(i-1) + m; 
  %P_int(tt(i))
  %[x, fval, info] = fsolve (@(x)NonlinearSistem(x,P_int(tt(i)), gamma_theta(i),a_0,c_1), [0.5; 2]);
  [x, fval, info] = fsolve (@(x)PEq(x,P_int(tt(i)), gamma_theta(i),a_0, mu, Jm, h_0), 1)
  alpha(i) = x;
  a(i) = alpha(i)*a_0*gamma_theta(i);
  end 

  a_mm = a*1000;

