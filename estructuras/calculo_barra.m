function [Sj,alpha,T,L]=calculo_barra(coord_nodo_i,coord_nodo_f,E,A)
    V=coord_nodo_f-coord_nodo_i;
    L = sqrt(V(1)^2+V(2)^2);
    alpha = atan2d(V(2),V(1));
    Sj=E*A/L*[(cosd(alpha))^2  (cosd(alpha))*(sind(alpha))   -(cosd(alpha))^2 -(cosd(alpha))*(sind(alpha));...
       (cosd(alpha))*(sind(alpha))  (sind(alpha))^2   -(cosd(alpha))*(sind(alpha)) -(sind(alpha))^2; -(cosd(alpha))^2 -(cosd(alpha))*(sind(alpha)) ...
       (cosd(alpha))^2  (cosd(alpha))*(sind(alpha));-(cosd(alpha))*(sind(alpha)) -(sind(alpha))^2    (cosd(alpha))*(sind(alpha))  (sind(alpha))^2];
     T=[cosd(alpha) sind(alpha) 0 0;0 0 cosd(alpha) sind(alpha)];
end