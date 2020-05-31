clear
clc

L = 2;
n = 10:10:100; Nt = 3;
tau = 0.005;

t = zeros(size(n),2);
%%%% построение сетки %%%%%%%%%%%%%%%%%%%%
for i=1:size(n,2)
  n1 = n(i);
  h = L/(n1+1);
  N=n1+2;
  x = -cos(((1:N)-1)*pi/(N-1));
  x = x(2:end-1);
  [Y,X] = ndgrid(x,x);

%%%% задание матриц для решения в виде AU=F %%%%%%%%%
  F = 10*sin(2*pi*X/L).*sin(2*pi*Y/L);
  U = zeros(n1,n1); % начальное условие
  C = gallery('chebspec',N); %генерирование матрицы диф-ния спектр.метода Чебышева
  C = C(2:end-1,2:end-1);
  E = speye(n1);
  A = kron(E,C)+kron(C,E); %Krank-Nicolson matrix
  A1 = C; %ADI matrix

  %% Crank-Nicolson %%
  E = eye(size(A));
  Apos = E + 0.5*tau*A;
  Aneg = E - 0.5*tau*A;
  U1 = U;

  tic
  for k = 1:Nt
    f = Apos*U1(:) + tau*F(:);
    U1 = Aneg\f(:);
  end
  t(i,1) = toc;
  
  %% ADI %%
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
  t(i,2) = toc;
  
end

U1=reshape(U1,n1,n1);
U2=reshape(U2,n1,n1);

%%==========Crank-Nicolson Scheme====
subplot('position',[0.06,0.4, 0.41 0.5]), mesh(X,Y,U1);
title('Crank-Nicolson Scheme');
xlabel('x'); ylabel('y'); zlabel('z');
text(-3.2,-2.5,['U1 time: ', num2str(t(i,1),3)]);
%%==========ADI Scheme===============
subplot('position',[0.56,0.4, 0.41 0.5]), mesh(X,Y,U2);
title('ADI Scheme');
xlabel('x'); ylabel('y'); zlabel('z');
text(-3.2,-2.5,['U2 time: ', num2str(t(i,2),3)]);
subplot('position',[0.3,0.15, 0.41 0.2]),
text( -0.10, 0.3,'Relative difference between the solutions:');
D=2*norm(U2(:)-U1(:))/norm(U2(:)+U1(:));
text( 0.1, -0.0,['||U_1-U_2||/||U_1+U_2|| = ', num2str(D,3)]);
axis off


figure
plot(n,t(:,1),"r-*",n,t(:,2),"b-*")
title("Time comparison for spectral method")
legend({'Crank-Nicolson','ADI'},'Location','northwest');
xlabel("amount of nodes"); ylabel("time");