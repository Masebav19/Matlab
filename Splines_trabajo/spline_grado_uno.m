function [terminos] = spline_grado_uno(X,Y)
   indices = 1:length(X)-1;
   indices_f=2*indices-1;
   indices_c=2*indices;
   M_f=zeros(2,2,length(X)-1);
   Y_f=zeros(2,1,length(X)-1);
   for i=1:length(X)-1
      M_f(:,:,i)=[X(i) 1;X(i+1) 1];
      Y_f(:,:,i)=[Y(i);Y(i+1)];
   end
   M_total=zeros(2*(length(X)-1),2*(length(X)-1));
   B=zeros(2*(length(X)-1),1);
   for i=1:length(indices)
        M_total(indices_f(i):indices_c(i),indices_f(i):indices_c(i))=M_f(:,:,i);
        B(indices_f(i):indices_c(i))= Y_f(:,:,i);
   end
   
   terminos=M_total^-1*B;
end

% X=[2 7 9]
% Y=[1 3 0]
% N funciones = 2
% f1(x)= a1x+b1
% f2(x)= a2x+b2
% f1(x1)= a1*2+b1 = 1
% f1(x2)= a1*7+b1 = 3
% f2(x2)= a2*7+b2 = 3
% f2(x3)= a2*9+b2 = 0

