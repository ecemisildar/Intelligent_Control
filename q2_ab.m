clear;clc;

syms x1 x2 x3

en = 0.1*0.15; 
en_dot = 0;
%genetic coeffs
Ke = 0.15;
Kde = 1.5;
alfa = 0.7;
beta = 1.2;

y_old = 0.5;
y = y_old;
uPID = 0;
uc_old = 0;
etr_old = 0;
ufc_ = 0;

% e=0.2;
% e_dot=0.3;
e = [-1 -0.4 0 0.4 1];
e_dot = [-1 -0.4 0 0.4 1];
u = [-1 -0.7 -0.5 -0.3 0;
    -0.7 -0.4 -0.2 0 0.3;
    -0.5 -0.2 0 0.2 0.5;
    -0.3 0 0.2 0.4 0.7;
     0 0.3 0.5 0.7 1];
 
%reference signal 
for k = 1:500
    if k<100
        r(k) = 0.4;
    elseif k>=100 && k<200
        r(k) = 0.5;
    elseif k>=200 && k<300
        r(k) = 0.6;
    elseif k>=300 && k<400
        r(k) = 0.5;
    elseif k>=400 && k<=500
        r(k) = 0.4;
    end


%fuzzy part 
for i=1:4
    for j=1:4
        %A calculation for 5 interval
        A(i) = (e(i+1)-en)/(e(i+1)-e(i));
        A(i+1) = (en-e(i))./(e(i+1)-e(i));
        %B calculation for 5 interval
        B(j) = (e(j+1)-en_dot)./(e_dot(j+1)-e_dot(j));
        B(j+1) = (en_dot-e(j))./(e_dot(j+1)-e_dot(j));
        
        ufc = A(i).*B(j)*u(i,j)+A(i+1).*B(j)*u(i+1,j)+A(i).*B(j+1)*u(i,j+1)+A(i+1).*B(j+1)*u(i+1,j+1);
    end       
end
%control signal
uc = alfa*ufc + beta*(ufc + ufc_);
uPID = uPID + uc;
ufc_ = ufc;

%runge kutta 
    if k>1
        [x1,x2,x3] = CSTR_runga_kutta_new(x1_old, x2_old, y(k-1), 3, 0.5, 1, 1, uPID, 0.1);
        y(k) = x3;
        
    else
        [x1,x2,x3] = CSTR_runga_kutta_new(0.31, 0.71, 0.5, 3, 0.5, 1, 1, uPID, 0.1);
    end
    x1_old = x1;
    x2_old = x2;


etr = r(k) - y(k);
edottr = (etr - etr_old);

en = Ke*etr;
en_dot = Kde*edottr;

%update part
uc_old = uc;   
etr_old = etr;  

end
figure(1)
plot(y,'r','LineWidth',1.5);
hold on;
plot(r,'b','LineWidth',1.5);
hold on; grid on
title('Fuzzy controller')
legend('output signal','reference signal')



