function [Sj,alpha,T,L]=calculo_barra_param(coord_nodo_i,coord_nodo_f,E,A)
    V=coord_nodo_f-coord_nodo_i;
    L = sqrt(V(1)^2+V(2)^2);
    alpha = atan2d(V(2),V(1));
    T=[cosd(alpha) sind(alpha) 0 0 0 0;-sind(alpha) cosd(alpha) 0 0 0 0;0 0 1 0 0 0;...
        0 0 0 cosd(alpha) sind(alpha) 0;0 0 0 -sind(alpha) cosd(alpha) 0;0 0 0 0 0 1];
    Sm=E*A/L*[1 0 0 -1 0 0;0 12 6 0 -12 6;0 6 4 0 -6 2;-1 0 0 1 0 0;0 -12 -6 0 12 -6;0 6 2 0 -6 4];
    Sj=T'*Sm*T;
end