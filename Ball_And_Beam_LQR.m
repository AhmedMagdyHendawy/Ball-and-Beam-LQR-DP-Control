%% Auther
% Auther : Ahmed Magdy Hendawy
% Date : 20.02.2019
clear;
close all;
clc ;
%% Linearization of the system
syms x1 x2 x3 x4 J m r J_b g tau
ps = 0.1; % The operating state x1
Ts=0.1;   % The sampling Time
tau_s=m*g*ps; % The operating control input
x_s=[ps;0;0;0]; % The operating state
f1=x2;
f2=(1/((J_b/r^2)+m))*(m*x1*(x4^2)-m*g*sin(x3));
f3=x4;
f4=(1/(m*(x1^2)+J+J_b))*(tau - m*g*x1*cos(x3)-2*m*x1*x2*x4);
% The non-linear symbolic system
f=[f1;f2;f3;f4];
% The non-linear numerical system
f_valued=subs(f,{J,m,r,J_b,g},{9.99*10^-4,0.11,0.015,9.99*10^-4,9.81});
% partial derivative of f4 w.r.t x1
df4_dx1=diff(f4,x1)
% The state space representation of the linearized system
A=[0 1 0 0; 0 0 -m*g/((J_b/r^2)+m) 0 ; 0 0 0 1; subs(df4_dx1,{x1,x2,x3,x4,tau},{ps,0,0,0,tau_s}) 0 0 0];
B=[0; 0; 0; 1/(m*ps^2+J+J_b)];
C=[1 0 0 0];
D=0;
A_valued=double(vpa(eval(subs(A,{J,m,r,J_b,g},{9.99*10^-4,0.11,0.015,9.99*10^-4,9.81}))));
B_valued=double(vpa(eval(subs(B,{J,m,r,J_b,g},{9.99*10^-4,0.11,0.015,9.99*10^-4,9.81}))));
Linearized_Sys=ss(A_valued,B_valued,C,D);
%% Discretization of the linearized system
% The state space representation of the discretized system
Discretized_Sys=c2d(Linearized_Sys,Ts);
[A_d,B_d,C_d,D_d]=ssdata(Discretized_Sys)
% Check the stability
e=eig(A_d) % The first eigenvalue is greater than 1, then the system is unstable
% Controlabilty Check
c_rank=rank(ctrb(Discretized_Sys));
controlabilty=c_rank && rank(A_d); % If the rank equals to rank of A matrix, then the system is 
                                   % controllable which is the case.
% Observability Check
O_rank=rank(obsv(Discretized_Sys));
observability=O_rank && rank(A_d); % If the rank equals to rank of A matrix, then the system is 
                                   % observable which is the case.

% The system is unstable, controllable and observable

%% solving the Difference Riccati Equation
x0=[0.25;0;0.3;0];
N=500;
P=[1 0 0 0;0 1 0 0; 0 0 1 0; 0 0 0 1];
Q=[1 0 0 0;0 1 0 0; 0 0 1 0; 0 0 0 1];
R=0.1;

Pback(:,:,1)=P;
for i = 1:N
    Pback(:,:,i+1)= (A_d')*Pback(:,:,i)*A_d + Q ...
        - (A_d')*Pback(:,:,i)*B_d*inv((B_d')*Pback(:,:,i)*B_d+R)*(B_d')*Pback(:,:,i)*A_d;
end

Pfor=[];
for i =1:N+1
    Pfor(:,:,i)=Pback(:,:,N+2-i);
end

%% Simulate both the linear and non-liear systems
x=x0-x_s; % x0 is the actual initial condition state without deviation
time=Ts*(0:1:N);
cost_to_go=[];
cost_to_go(1)=((x0-x_s)')*Pfor(:,:,1)*((x0-x_s));

cost_to_go_non=[];
cost_to_go_non(1)=((x0-x_s)')*Pfor(:,:,1)*((x0-x_s));
x_non=x0;


for i = 1:N
    % Linear system simulation
    % x is the deviation variable
    u(i)= (-inv((B_d')*Pfor(:,:,i+1)*B_d+R)*(B_d')*Pfor(:,:,i+1)*A_d)*x(:,i);
    x(:,i+1) = A_d*x(:,i)+B_d*u(i);
    cost_to_go(i+1)=(x(:,i+1)')*Pfor(:,:,i+1)*(x(:,i+1));
    
    % Non-linear system simulation
    % x_non is the actual state
    u_non(i)=(-inv((B_d')*Pfor(:,:,i+1)*B_d+R)*(B_d')*Pfor(:,:,i+1)*A_d)*(x_non(:,i)-x_s);
    x_non(:,i+1) = x_non(:,i)+Ts*double(vpa(eval(subs(f_valued,{x1,x2,x3,x4,tau},{x_non(1,i),x_non(2,i),x_non(3,i),x_non(4,i),u_non(i)+(0.11*9.81*ps)}))));
    cost_to_go_non(i+1)=((x_non(:,i+1)-x_s)')*Pfor(:,:,i+1)*((x_non(:,i+1)-x_s));
end
%% Plotting the responses, control effort and cost to go for both the linear and non-linear systems
figure;
hold on; 
title('The Linearized System States');
plot(time, x(1,:), "o-b","Linewidth",1.5);
plot(time,x(2,:),"o-r","Linewidth",1.5);
plot(time, x(3,:), "o-g","Linewidth",1.5);
plot(time,x(4,:),"o-y","Linewidth",1.5);
legend('x_1','x_2','x_3','x_4');
xlabel("Time (sec)")
ylabel('linearized states');
grid on;
hold off ;
figure;
hold on;
title('Cost To Go in Case of Linearized System');
plot(time,cost_to_go,"o-r",'Linewidth',1.5);
ylabel('Cost to go');
xlabel('Time (sec)');
legend('cost to go');
grid on;
hold off ;
figure;
hold on;
title('Control Effort in Case of Linear System');
ts=timeseries(u,time(1:end-1));
ts=setinterpmethod(ts,"zoh");
plot(ts,'o-b','Linewidth',1.5);
ylabel('Control Effort');
xlabel('Time (sec)');
legend('control effort');
grid on;
hold off ;
figure;
hold on; 
title('The Non-Linear System States');
plot(time, x_non(1,:), "o-b","Linewidth",1.5);
plot(time,x_non(2,:),"o-r","Linewidth",1.5);
plot(time, x_non(3,:), "o-g","Linewidth",1.5);
plot(time,x_non(4,:),"o-y","Linewidth",1.5);
legend('x_1','x_2','x_3','x_4');
xlabel("Time (sec)")
ylabel('linearized states');
grid on;
hold off ;
figure;
hold on;
title('Cost To Go in Case of Non-Linear System');
plot(time,cost_to_go_non,"o-r",'Linewidth',1.5);
ylabel('Cost to go');
xlabel('Time (sec)');
legend('cost to go');
grid on;
figure;
hold on;
title('Control Effort in NLS');
ts=timeseries(u_non,time(1:end-1));
ts=setinterpmethod(ts,"zoh");
plot(ts,'o-b','Linewidth',1.5);
ylabel('Control Effort');
xlabel('Time (sec)');
legend('control effort');
grid on;
hold off ;