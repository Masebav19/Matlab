function [Aj]=calculo_carga_distr(data,L)
Aj=zeros(6,1);
% data-> 1x5 
%1 Enable
%2 Carga
%3 Distancia a
%4 Distancia b
%5 situación [1-8]
%Aj [xi,yi,Mi],[xf,yf,Mf]
    switch data(5)
        case 1
           Aj(2,1)= data(1)*data(2)*data(3)*data(4)^2/L^2;
           Aj(3,1)=(data(1)*data(2)*data(4)^2/L^3)*(3*data(3)+data(4));
           Aj(5,1)= -data(1)*data(2)*data(3)^2*data(4)/L^2;
           Aj(6,1)=(data(1)*data(2)*data(3)^2/L^3)*(data(3)+3*data(4));
        case 2
           Aj(2,1)= (data(1)*data(2)*data(4)/L^2)*(2*data(3)-data(4));
           Aj(3,1)=6*data(1)*data(2)*data(3)*data(4)/L^3;
           Aj(5,1)= (data(1)*data(2)*data(3)/L^2)*(2*data(4)-data(3));
           Aj(6,1)=-Aj(2,1);
        case 3
           Aj(1,1)= -data(1)*data(2)*data(4)/L;
           Aj(4,1)=-data(1)*data(2)*data(3)/L;
        case 4
           Aj(1,1)=-data(1)*data(2)*data(4)/L;
           Aj(4,1)=-data(1)*data(2)*data(3)/L;
        case 5
           Aj(2,1)= data(1)*data(2);
           Aj(3,1)=(data(1)*data(2)*data(3)/L)*(L-data(3));
           Aj(5,1)= Aj(2,1);
           Aj(6,1)=-Aj(3,1); 
        case 6
           Aj(2,1)= data(1)*data(2)*L/2;
           Aj(3,1)= data(1)*data(2)*L^2/12;
           Aj(5,1)= Aj(2,1);
           Aj(6,1)=-Aj(3,1); 
        case 7
           Aj(2,1)= (data(1)*data(2)*data(3)/(2*L^3))*(2*L^3-2*data(3)^2*L+3*data(3)^3);
           Aj(3,1)=(data(1)*data(2)*data(3)^2/(12*L^2))*(6*L^2-8*data(3)*L+3*data(3)^2);
           Aj(5,1)= (data(1)*data(2)*data(3)^3/(2*L^3))*(2*L-data(3));
           Aj(6,1)=-(data(1)*data(2)*data(3)^3/(12*L^2))*(4*L-3*data(3));  
        case 8
           Aj(2,1)= data(1)*3*data(2)*L/20;
           Aj(3,1)= data(1)*data(2)*L^2/30;
           Aj(5,1)= data(1)*7*data(2)*L/20;
           Aj(6,1)=-data(1)*data(2)*L^2/30; 
    end
end