function [terminos,M_total,B]=spline_grado_tres(X,Y)
   M_f=zeros(2,4,length(X)-1);
   B=zeros(1,4*(length(X)-1));
   for i=1:length(X)-1
        M_f(:,:,i)=[X(i)^3 X(i)^2 X(i) 1;X(i+1)^3 X(i+1)^2 X(i+1) 1];
        B(2*i-1:2*i)=[Y(i),Y(i+1)];
   end
   M_total = zeros(4*(length(X)-1),4*(length(X)-1));
   for i=1:length(X)-1
        M_total(2*i-1:2*i,4*i-3:4*i)=M_f(:,:,i);
   end

   M_derivate=zeros(1,4*(length(X)-1),2*(length(X)-2));
   %Primera derivada
   for i=1:length(X)-2
       M_derivate(1,4*i-3:4*i+4,i)=[3*X(i+1)^2 +2*X(i+1) 1 0 -3*X(i+1)^2 -2*X(i+1) -1 0];
   end
   %Segunda derivada
   for i=length(X)-1:2*(length(X)-2)
       M_derivate(1,4*(i-(length(X)-2))-3:4*(i-(length(X)-2))+4,i)=...
           [6*X((i-(length(X)-2))+1) 2 0 0 ...
           -6*X((i-(length(X)-2))+1) -2 0 0];
   end
    k=1;
   for i=2*(length(X)-1)+1:4*(length(X)-1)-2
        M_total(i,:)=M_derivate(1,:,k);
        k=k+1;
   end
   M_total(end-1,1:4)=[6*X(1) 2 0 0];
   B(end-1)=-0.1;
   M_total(end,end-3:end)=[6*X(end) 2 0 0];
   B(end)=-0.1;
   terminos=M_total^-1*B';
end