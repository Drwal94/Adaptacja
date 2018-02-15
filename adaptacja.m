clc
clear all
close all
A = poly([-0.32-0.67i -0.32+0.67i -0.18]);
B = poly([0.77i -0.77i])
C = poly(0.8);
N = 200;
y = randn(1,5)*1;
v = randn(1,505)*1;
Y_G = y;
U_G = zeros(1,5);
Y_0 = y;
U_0 = zeros(1,5);
Y_iden = y;
U_iden = zeros(1,5);

fi = [U_G(4), U_G(3), U_G(2), Y_G(4), Y_G(3), Y_G(2)];
theta = [B(1), B(2), B(3), -A(2)-0.5, -A(3), -A(4)];

for i = 1:N
    t = i+4;
    G = theta(4:6);
    U_G(t) = sum(-G(3:-1:1).*Y_G(t-2:t)) + sum(-B(3:-1:2).*U_G(t-2:t-1));
    Y_G(t+1) = sum(-A(4:-1:2).*Y_G(t-2:t)) + sum(B(3:-1:1).*U_G(t-2:t)) + sum(C(2:-1:1).*v(t:t+1));
    Y_0(t) = -sum(A(4:-1:2).*Y_0(t-3:t-1)) + sum(B(3:-1:1).*U_0(t-3:t-1)) + sum(C(2:-1:1).*v(t-1:t));
    U_0(t) = 0;
end

theta = zeros(6,1);
fi = zeros(4,6);

uf = U_iden;
yf = U_iden;

for i = 1:N
    t = i+4;
    G = theta(4:6)';
    bt = theta(1:3)';
    U_iden(t) = sum(-G(3:-1:1).*y(t-2:t)) + sum(-bt(3:-1:2).*U_iden(t-2:t-1));
    uf(t) = U_iden(t)-C(2)*uf(t-1);
    yf(t) = y(t)-C(2)*yf(t-1);
    fi = [fi; [uf(t), uf(t-1), uf(t-2), yf(t), yf(t-1), yf(t-2)]];
    if(t > 6)
        theta = lsqlin(fi(1:t-1,:),y(2:t), eye(6),ones(6,1), [], [], -2*ones(6,1), ones(6,1));
    end
    y(t+1) = sum(-A(4:-1:2).*y(t-2:t)) + sum(B(3:-1:1).*U_iden(t-2:t)) + sum(C(2:-1:1).*v(t:t+1));
   
    Y_iden = [Y_iden, fi(i,:)*theta];
end

subplot(4,1,1)
plot(U_iden)
title('Sygna³ steruj¹cy')
subplot(4,1,4)
plot(Y_0,'b')
title('Wyjœcie z obiektu przu U = 0')
subplot(4,1,3)
plot(Y_G,'r')
title('Wyjœcie z obiektu dla wyliczonego G')
subplot(4,1,2)
plot(y,'g')
title('Zidentyfikowane wyjœcie z obiektu')
[var(y(100:end)), var(Y_G(100:end)), var(Y_0(100:end))]