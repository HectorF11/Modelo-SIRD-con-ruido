function []=ruido2(N,e)
X0=[1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1];  
Xn=[1; 0; 0; 0];
X=[Xn];
A=[0.95 0.04 0 0; 0.05 0.85 0 0; 0 0.1 1 0; 0 0.01 0 1];
for k=1:1:N-1
  Xn=A*Xn;
  X=[X,Xn];
endfor

  X=X+e*randn(4,N); 
  Y=[X,A*Xn];
  Y(:,1)=[];

        #Restricciones de igualdad
         function r1 = g (x)
            r1 = [ x(1)+x(2)-1;x(3);x(4);
                  x(5)+x(6)+x(7)+x(8)-1;
                  x(11)-1;x(9);x(10);x(12);
                  x(16)-1;x(13);x(14);x(15)];
          endfunction
          
          function obj = phi (x)
           obj=norm(reshape(x,4,4)*X-Y,'fro')^2; 
          endfunction
          
          #Restricciones de desigualdad
          function r2 = h (x)
           r2=[x(1); x(2); x(3); x(4); x(5); x(6); x(7); x(8); x(9); x(10); x(11); x(12); x(13); x(14); x(15); x(16);
           x(1)-x(2)-x(3)-x(4);-x(5)+x(6)-x(7)-x(8);-x(9)-x(10)+x(11)-x(12);-x(13)-x(14)-x(15)+x(16)]; 
          endfunction
          
          


[a, obj, info, iter, nf, lambda] = sqp (X0, @phi, @g, @h, zeros(1,16), ones(1,16),100,2e-10);
a=reshape(a,4,4)

iter



   disp(['||Ae X-Y||_F = ',num2str(norm(a*X-Y,'fro'))])
   disp(['||Ae-A||_F = ',num2str(norm(a-A,'fro'))])

X1=[1; 0; 0; 0];
XX=[X1];
X2=X1;
XXX=[X2];

for k=1:1:1000
  X1=a*X1;
  XX=[XX,X1];
  X2=A*X2;
  XXX=[XXX,X2];
endfor

subplot(311);
plot((XXX+e*randn(4,1001))');
title('Órbitas del modelo con ruido')
axis([0, 1000, -0.1, 1])
set(gca, 'FontSize', 14)

subplot(312);
plot(XXX');
title('Órbitas del modelo original')
axis([0, 1000, -0.1, 1])
set(gca, 'FontSize', 14)

subplot(313);
plot(XX');
title('Órbitas del modelo identificado')
axis([0, 1000, -0.1, 1])
set(gca, 'FontSize', 14)
endfunction