clear;clc;
%referance signal
load referans.mat referans Ts time 

r = referans;
r = downsample(r,100);

global n;
n = 100;

%inital values
w11 = -1.7502*ones(1,51);
w12 = -0.8314*ones(1,51);
w13 = -1.1564*ones(1,51);
w21 = -0.2857*ones(1,51);
w22 = -0.9792*ones(1,51);
w23 = -0.5336*ones(1,51);
w31 = -2.0026*ones(1,51);
w32 = 0.9642*ones(1,51);
w33 = 0.5201*ones(1,51);
W = [w11 w12 w13;
     w21 w22 w23;
     w31 w32 w33];
b = 10e-5;
mu = 0.01;
ts = 0.01;

x = zeros(3,51);
y = 0;
u = 0;
 for k=1:n
    %forward propagation    
    x(1,:) = r.*w11 + y.*w21;
    x(2,:) = r.*w12 + y.*w22;
    x(3,:) = r.*w13 + y.*w23;
    
    for a=1:3
        v(a,:) = (exp(x(a,:))-exp(-x(a,:)))./(exp(x(a,:))+exp(-x(a,:)));
    end    
    u = w31.*v(1,:) + w32.*v(2,:) + w33.*v(3,:) + b;
        
    u_prev = u; 
    y_prev = y;
    %plant
    y = (1.2*(1-0.8*exp(-0.1*k)).*y)/(1+y.^2) + u_prev;
    %tracking error
    e = r - y_prev;
    %objective function
    J = objective(e);
    
    %backward phase
        %output layer
        if round(u-u_prev) == 0
            delta_w31 = mu*e.*(1-v(1,:).^2).*x(1,:);
            delta_w32 = mu*e.*(1-v(2,:).^2).*x(2,:);
            delta_w33 = mu*e.*(1-v(3,:).^2).*x(3,:);            
        else
            delta_w31 = mu*e.*(1-v(1,:).^2).*x(1,:).*(y-y_prev)./(u-u_prev);
            delta_w32 = mu*e.*(1-v(2,:).^2).*x(2,:).*(y-y_prev)./(u-u_prev);
            delta_w33 = mu*e.*(1-v(3,:).^2).*x(3,:).*(y-y_prev)./(u-u_prev);                   
        end
        
        delta_b = -mu*e.*x(1,:);
        b = b + delta_b;
        
        w31 = w31 + delta_w31;
        w32 = w32 + delta_w32;
        w33 = w33 + delta_w33;
        
       %input layer           
%        delta_w11 = mu*(1-v(1,:).^2).*(e.*(1-v(1,:).^2).*w31).*r;
%        delta_w12 = mu*(1-v(2,:).^2).*(e.*(1-v(2,:).^2).*w32).*r;
%        delta_w13 = mu*(1-v(3,:).^2).*(e.*(1-v(3,:).^2).*w33).*r;
%        
%        delta_w21 = mu*(1-v(1,:).^2).*(e.*(1-v(1,:).^2).*w31).*y;
%        delta_w22 = mu*(1-v(2,:).^2).*(e.*(1-v(2,:).^2).*w32).*y;
%        delta_w23 = mu*(1-v(3,:).^2).*(e.*(1-v(3,:).^2).*w33).*y;
       delta_w11 = mu*(e.*(1-v(1,:).^2).*w31).*r;
       delta_w12 = mu*(e.*(1-v(2,:).^2).*w32).*r;
       delta_w13 = mu*(e.*(1-v(3,:).^2).*w33).*r;
       
       delta_w21 = mu*(e.*(1-v(1,:).^2).*w31).*y;
       delta_w22 = mu*(e.*(1-v(2,:).^2).*w32).*y;
       delta_w23 = mu*(e.*(1-v(3,:).^2).*w33).*y;
       
       w11 = w11 + delta_w11;
       w12 = w12 + delta_w12;
       w13 = w13 + delta_w13;
       
       w21 = w21 + delta_w21;
       w22 = w22 + delta_w22;
       w23 = w23 + delta_w23;
       

 end
 
figure(1)
subplot(3,3,1)
plot(w11,'g','LineWidth',1.5); hold on;
title('w11')
subplot(3,3,2)
plot(w12,'g','LineWidth',1.5); hold on;
title('w12')
subplot(3,3,3)
plot(w13,'g','LineWidth',1.5); hold on;
title('w13')
subplot(3,3,4)
plot(w21,'g','LineWidth',1.5); hold on;
title('w21')
subplot(3,3,5)
plot(w22,'g','LineWidth',1.5); hold on;
title('w22')
subplot(3,3,6)
plot(w23,'g','LineWidth',1.5); hold on;
title('w23')
subplot(3,3,7)
plot(w31,'g','LineWidth',1.5); hold on;
title('w31')
subplot(3,3,8)
plot(w32,'g','LineWidth',1.5); hold on;
title('w32')
subplot(3,3,9)
plot(w33,'g','LineWidth',1.5); hold on;
title('w33')
 
figure(2)
plot(time,referans,'linewidth',2);hold on
axis([min(time) max(time) min(referans)-0.5 max(referans)+0.5])

plot(y,'g','LineWidth',2); hold on;
xlabel('Time(sec)')
legend('r(t)','y(t)') 
figure(3)
plot(u,'m','linewidth',2);hold on
xlim([0 50])
legend('u(t) control signal')
xlabel('Time(sec)')
