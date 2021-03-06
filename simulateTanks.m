%% Linearization
clc
clearvars
simlength=500;
dt=0.1;
m_size=simlength/dt;

%constants
Tc=20; Th=65; Td=30;
C=0.3; a=9;
tauc=100; tau=40; 

%variables
Fc=31; Fh=20; Fd=10;
h=45.94; T=36.39;
V=C.*h.^2;
V=round(V,2);

%workpoint
u1_d=Fh;
u2_d=Fc;
v1_d=Fd;
x1_d=V;
x2_d=T;

%input range
u_jumps=0:0.5:3;

%alocating memory
x1=ones(size(u_jumps,2),m_size)*V;
x2=ones(size(u_jumps,2),m_size)*T;
y1=ones(size(u_jumps,2),m_size)*h;
y2=ones(size(u_jumps,2),m_size)*T;


x1_lin=ones(size(u_jumps,2),m_size)*V;
x2_lin=ones(size(u_jumps,2),m_size)*T;
y1_lin=ones(size(u_jumps,2),m_size)*h;
y2_lin=ones(size(u_jumps,2),m_size)*T;

legends = cell(size(u_jumps));


for i=1:size(u_jumps,2)

%setting inputs    
u1=u_jumps(i)*Fh*ones(1,m_size);
u2=u_jumps(i)*Fc*ones(1,m_size);
v1=u_jumps(i)*Fd*ones(1,m_size);

legends{i} = "u_{1} = " + num2str(u1(1)) + " u_{2} = " + num2str(u2(1)) + " v_{1} = " + num2str(v1(1));
    
for k =tauc/dt+2:m_size
    
    %Linear model

    x1_lin(i,k)=(u1_d+u2_d+v1_d-a*(x1_d./C)^0.25)*dt+x1_d+dt*(u1(k-1)-u1_d)...
        +dt*(u2(k-1)-u2_d)+dt*(v1(k-1)-v1_d)+(-a*dt*((x1_d/C)^0.25)/(4*x1_d)+1)*(x1_lin(i,k-1)-x1_d);
    x1_lin=round(x1_lin,2);
    
    x2_lin(i,k)=((u1_d*Th+u2_d*Tc+v1_d*Td-(u1_d+u2_d+v1_d)*x2_d)/x1_d)*dt+x2_d...
        +(Th*dt/x1_d-x2_d*dt/x1_d)*(u1(k-1)-u1_d)+(Tc*dt/x1_d-x2_d*dt/x1_d)*(u2(k-1)-u2_d)...
        +(Td*dt/x1_d-x2_d*dt/x1_d)*(v1(k-1)-v1_d)...
        +(-(dt/x1_d^2)*(u1_d*Th+u2_d*Tc+v1_d*Td-(u1_d+u2_d+v1_d)*x2_d))*(x1_lin(i,k-1)-x1_d)...
        +(-(dt/x1_d)*(u1_d+u2_d+v1_d)+1)*(x2_lin(i,k-1)-x2_d);
    x2_lin=round(x2_lin,2);
    
    y1_lin(i,k)=sqrt(x1_d/C)+1/(2*sqrt(x1_d*C))*(x1_lin(i,k)-x1_d);
    y1_lin(i,k)=round(y1_lin(i,k),2);
    y2_lin(i,k)=x2_lin(i,k-tau/dt);
    
    %Nonlinear model
    
    x1(i,k)=(u1(k-1)+u2(k-1-tauc/dt)+v1(k-1)-a*(x1(i,k-1)./C).^0.25)*dt+x1(i,k-1);
    x1(i,k)=round(x1(i,k),2);
    if x1(i,k)<0 %if volume < 0
        x1(i,k)=0;
    end
    
    x2(i,k)=((u1(k-1)*Th+u2(k-1-tauc/dt)*Tc+v1(k-1)*Td-(u1(k-1)+u2(k-1-tauc/dt)+v1(k-1))*x2(i,k-1))/x1(i,k-1))*dt+x2(i,k-1);
    x2(i,k)=round(x2(i,k),2);

    y1(i,k)=round(sqrt(x1(i,k)/C),2);
    y2(i,k)=x2(i,k-tau/dt);  
end

end

%% Drawing figures

for i=1:size(u_jumps,2)
figure(1)    
stairs(dt:dt:simlength,y1(i,:))
hold on
figure(2)
stairs(dt:dt:simlength,y2(i,:))
hold on
figure(3)    
stairs(dt:dt:simlength,x1(i,:))
hold on
figure(4)
stairs(dt:dt:simlength,x2(i,:))
hold on

end
figure(2)
for i=1:size(u_jumps,2)
figure(1)
stairs(dt:dt:simlength,y1_lin(i,:),'--')
hold on
figure(2)
stairs(dt:dt:simlength,y2_lin(i,:),'--')
figure(3)
stairs(dt:dt:simlength,x1_lin(i,:),'--')
hold on
figure(4)
stairs(dt:dt:simlength,x2_lin(i,:),'--')
hold on

end
figure(1)
grid on
legend(legends, 'Location', 'EastOutside')
xlabel("t [s]")
ylabel("y_{1}")
title ("Przebiegi wyj�cia y_{1} modelu linowego i nieliniowego dla r�nych skok�w sterowania")
set(gcf,'WindowState','fullscreen')
print(gcf,"./fig/y1" + "_u1=" + num2str(u1(1)) + "_u2=" + num2str(u2(1)) + "_v1=" + num2str(v1(1)) + "_dt=" + num2str(dt)+".eps" , '-depsc')

figure(2)
grid on
legend(legends, 'Location', 'EastOutside')
xlabel("t [s]")
ylabel("y_{2}")
title ("Przebiegi wyj�cia y_{2} modelu linowego i nieliniowego dla r�nych skok�w sterowania")
set(gcf,'WindowState','fullscreen')
print(gcf,"./fig/y2" + "_u1=" + num2str(u1(1)) + "_u2=" + num2str(u2(1)) + "_v1=" + num2str(v1(1)) + "_dt=" + num2str(dt)+".eps" , '-depsc')

figure(3)
grid on
legend(legends, 'Location', 'EastOutside')
xlabel("t [s]")
ylabel("x_{1}")
title ("Przebiegi zmiennej stanu x_{1} modelu linowego i nieliniowego dla r�nych skok�w sterowania")
set(gcf,'WindowState','fullscreen')
print(gcf,"./fig/x1" + "_u1=" + num2str(u1(1)) + "_u2=" + num2str(u2(1)) + "_v1=" + num2str(v1(1)) + "_dt=" + num2str(dt)+".eps" , '-depsc')

figure(4)
grid on
legend(legends, 'Location', 'EastOutside')
xlabel("t [s]")
ylabel("x_{2}")
title ("Przebiegi zmiennej stanu x_{2} modelu linowego i nieliniowego dla r�nych skok�w sterowania")
set(gcf,'WindowState','fullscreen')
print(gcf,"./fig/x2" + "_u1=" + num2str(u1(1)) + "_u2=" + num2str(u2(1)) + "_v1=" + num2str(v1(1)) + "_dt=" + num2str(dt)+".eps" , '-depsc')




%% State space and Transfer functions

%Continuous model

A_matrix=[(-a*((x1_d/C)^0.25)/(4*x1_d)),0;...
    (-(1/x1_d^2)*(u1_d*Th+u2_d*Tc+v1_d*Td-(u1_d+u2_d+v1_d)*x2_d)),(-(dt/x1_d)*(u1_d+u2_d+v1_d))];

B_matrix=[1,1,1;(Th/x1_d-x2_d/x1_d),(Tc/x1_d-x2_d/x1_d),(Td/x1_d-x2_d/x1_d)];

C_matrix=[round((1/(2*sqrt(x1_d*C))),2), 0; 0, 1];

D_matrix=zeros(2,3);
    
state_space_model=ss(A_matrix,B_matrix,C_matrix,D_matrix);

state_space_model.InputDelay=[0;tauc;0];
state_space_model.OutputDelay=[0;tau];

transfer_function = tf(state_space_model)  %#ok<NOPTS>


%Discrete model

discrete_state_space_model=c2d(state_space_model,dt);

discrete_transfer_function=tf(discrete_state_space_model)  %#ok<NOPTS>


