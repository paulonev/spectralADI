%% Двумерное уравнение теплопроводности
%% Схемы Кранка-Николсон и Метода Переменных Направлений
clear

L = 2;
n = 145; Nt = 3;
tau = 0.05;
%%%% построение сетки %%%%%%%%%%%%%%%%%%%%
h = L/(n+1);
x = -L/2+h : h : L/2-h;
[Y,X] = ndgrid(x,x);
F = 10*sin(2*pi*X/L).*sin(2*pi*Y/L); % AU=F
U = zeros(n,n); % начальное условие
A = -gallery('poisson',n )/h^2; % Матрица Пуассона

%% Неявная схема Кранка-Николсон %%
E = eye(size(A));
Apos = E + 0.5*tau*A;
Aneg = E - 0.5*tau*A;
U1 = U;
tic
for k = 1:Nt
  f = Apos*U1(:) + tau*F(:);
  U1 = Aneg\f(:);
end
CN = toc

%% ADI %%
e = ones(n,1);
A1 = spdiags([e,-2*e e], -1:1, n,n)/h^2; 
E = eye(size(A1));
Apos = E + 0.5*tau*A1;
Aneg = E - 0.5*tau*A1;
U2 = U;
tic
for k=1:Nt
  f = (Apos*(U2.'))' + 0.5*tau*F;
  U2 = Aneg\f;       % МПН x-направление
  f = Apos*U2 + 0.5*tau*F;
  U2 = (Aneg\f.').'; % МПН y-направление
end

ADI=toc
U1=reshape(U1,n,n);

%%==========Crank-Nicolson Scheme====
subplot('position',[0.06,0.4, 0.41 0.5]), mesh(X,Y,U1);
title('Crank-Nicolson Scheme');
xlabel('x'); ylabel('y'); zlabel('z');
text(-2.99,-2.5,['Computational time CN: ', num2str(CN,3)]);
%%==========ADI Scheme===============
subplot('position',[0.56,0.4, 0.41 0.5]), mesh(X,Y,U2);
title('ADI Scheme');
xlabel('x'); ylabel('y'); zlabel('z');
text(-2.99,-2.5,['Computational time ADI: ', num2str(ADI,3)]);
%%==========Difference===============
subplot('position',[0.3,0.15, 0.41 0.2]),
text( -0.10, 0.3,'Relative difference between the solutions:');
D=2*norm(U2(:)-U1(:))/norm(U2(:)+U1(:));
text( 0.1, -0.0,['||U_1-U_2||/||U_1+U_2|| = ', num2str(D,3)]);
axis off