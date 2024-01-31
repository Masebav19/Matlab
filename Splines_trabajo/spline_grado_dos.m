function [terminos]=spline_grado_dos(X,Y)
   M_f=zeros(2,3,length(X)-1);
   B=zeros(1,3*(length(X)-1)-1);
   c1=0;
   M_f(:,:,1)=[X(1)^2 X(1) 0;X(2)^2 X(2) 0];
   B(1:2)=[Y(1)-c1,Y(2)-c1];
   for i=2:length(X)-1
        M_f(:,:,i)=[X(i)^2 X(i) 1;X(i+1)^2 X(i+1) 1];
        B(2*i-1:2*i)=[Y(i),Y(i+1)];
   end
   M_total = zeros(3*(length(X)-1),3*(length(X)-1));
   for i=1:length(X)-1
        M_total(2*i-1:2*i,3*i-2:3*i)=M_f(:,:,i);
   end

   M_derivate=zeros(1,3*(length(X)-1),length(X)-2);
   M_derivate(1,1:6,1)=[2*X(2) 1 0 -2*X(2) -1 0];

   for i=2:length(X)-2
       M_derivate(1,3*i-2:3*i+3,i)=[2*X(i+1) 1 0 -2*X(i+1) -1 0];
   end
    k=1;
   for i=2*(length(X)-1)+1:3*(length(X)-1)-1
        M_total(i,:)=M_derivate(1,:,k);
        k=k+1;
   end
   % Condiciones extras
   M_total=M_total(1:end-1,:);
   M_total(:,3)=[];
    terminos=M_total^-1*B';
end