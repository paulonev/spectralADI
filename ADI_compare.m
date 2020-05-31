%% Двумерное уравнение теплопроводности
%% Схемы Переменных Направлений - 
%% с матрицами разностного и спектрального дифференцирования
clear

L = 2;
n = 100; Nt = 50;
tau = 0.01; h = L/(n+1);

%%%% построение сетки для разностной схемы %%%%%%%%%%%%%%%%%%%%
x_d = -L/2+h : h : L/2-h;
[Y_d,X_d] = ndgrid(x_d,x_d);

%%%% построение разностного оператора дифференцирования %%%%%
F_d = 10*sin(2*pi*X_d/L).*sin(2*pi*Y_d/L); % AU=F
U = zeros(n,n); % начальное условие
e = ones(n,1);
A_d = spdiags([e,-2*e e], -1:1, n,n)/h^2; 
E = eye(size(A_d));
Apos_d = E + 0.5*tau*A_d;
Aneg_d = E - 0.5*tau*A_d;
W_d = U;

tic
for k=1:Nt
  f = (Apos_d*(W_d.'))' + 0.5*tau*F_d;
  W_d = Aneg_d\f;       % МПН x-направление
  f = Apos_d*W_d + 0.5*tau*F_d;
  W_d = (Aneg_d\f.').'; % МПН y-направление
end
ADI_d=toc

%%%% построение сетки для спектральной схемы %%%%%%%%%%%%%%%%%%%%
N=n+2;
x_s = -cos(((1:N)-1)*pi/(N-1));
x_s = x_s(2:end-1);
[Y_s,X_s] = ndgrid(x_s,x_s);

%%%% построение спектрального оператора дифференцирования %%%%%
F_s = 10*sin(2*pi*X_s/L).*sin(2*pi*Y_s/L);
E = speye(n); 
A_s = gallery('chebspec',N); % матрица дифференцирования Чебышева
A_s = A_s*A_s;
A_s = A_s(2:end-1,2:end-1);

E = eye(size(A_s));
Apos_s = E + 0.5*tau*A_s;
Aneg_s = E - 0.5*tau*A_s;
W_s = U;

tic
for k=1:Nt
  f = (Apos_s*(W_s'))' + 0.5*tau*F_s;
  W_s = Aneg_s\f;       % МПН x-направление
  f = Apos_s*W_s + 0.5*tau*F_s;
  W_s = (Aneg_s\f').'; % МПН y-направление
end
ADI_s=toc

%%==========difference=============
subplot('position',[0.06,0.4, 0.41 0.5]), mesh(X_d,Y_d,W_d);
title('Difference ADI');
xlabel('x'); ylabel('y'); zlabel('T');
text(-3.2,-2.5,['Solution1 Time: ', num2str(ADI_d,3)]);
%%==========spectral===============
subplot('position',[0.56,0.4, 0.41 0.5]), mesh(X_s,Y_s,W_s);
title('Spectral ADI');
xlabel('x'); ylabel('y'); zlabel('T');
text(-3.0,-2.5,['Solution2 Time: ', num2str(ADI_s,3)]);

%%==========Difference===============
%subplot('position',[0.3,0.15, 0.41 0.2]),
%text( -0.10, 0.3,'Relative difference between the solutions:');
D=2*norm(ADI_s(:)-ADI_d(:))/norm(ADI_s(:)+ADI_d(:))
%text( 0.1, -0.0,['||U_1-U_2||/||U_1+U_2|| = ', num2str(D,3)]);
%axis off