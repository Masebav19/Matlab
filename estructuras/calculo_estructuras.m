% Calculo general de celocias
clc
clear variables

% Determinamos datos generales del modelo
E=200*10^9; %[Pa] - Acero
A=(8/100*8/100); %[m^2]
grados_liber=2; % Numero de grados de libertad por nudo

% Deteminamos los nudos
% Node2D(numero,coordenada,load,restrain)

% solicito el numero de nodos
Nnodos = input("Ingrese el numero de nodos: ");
% Pido ingresar la información de los nodos 
nodo_info = {};
for i=1:Nnodos
   msg=sprintf("Ingrese la información del nodo %d con el siguiente formato:\n{[Cooordx,Cooordy],[cargax,cargay],[restriccx,resticcy]} 0 para restringido y 1 para libre: ",i);
   nodo_info=[nodo_info;input(msg)];
end
% Pido ingresar el numero de Barras 
Nbarra = input("Ingrese el numero de barras: ");
barra_info=[];
% Pido la información de cada barra
for i = 1:Nbarra
   msg=sprintf("Ingrese la información de la barra %d con el siguiente formato:\n[nodo_i,nodo_f]: ",i);
   barra_info=[barra_info;input(msg)];
 
end
% claculo de los parámetros del elemento
SJ=zeros(4,4,Nbarra);%matriz de rigidez global
T=zeros(2,4,Nbarra);%Matrz de tranformación
alpha=zeros(1,Nbarra);%angulo de cada barra
for i=1:Nbarra
   [SJ(:,:,i),alpha(i),T(:,:,i),~]=calculo_barra(nodo_info{barra_info(i,1),1},nodo_info{barra_info(i,2),1},E,A);
end
%calculo de los índices para el calculo de la matriz de rigidez global de
%sistema
indices = zeros(Nbarra,4);
for i=1:Nbarra
   indices(i,:)=[grados_liber*barra_info(i,1)-1,grados_liber*barra_info(i,1),grados_liber*barra_info(i,2)-1,grados_liber*barra_info(i,2)];
end
%creo la matriz global del sistema
M_rigidez_global=zeros(grados_liber*Nnodos,grados_liber*Nnodos);
%con cada indice creo la matriz de cada barra dependiendo del nodo inicial
%y final
for i=1:Nbarra
    matriz_local=SJ(:,:,i);
    M_rigidez_global(indices(i,:),indices(i,:))= M_rigidez_global(indices(i,:),indices(i,:))+matriz_local;
end


% vector de cargas 
% vector de desplazamiento
indices_restringidos=zeros(1,grados_liber*Nnodos);
v_carga=zeros(grados_liber*Nnodos,1);
v_dezpla=zeros(grados_liber*Nnodos,1);
indice_total=1:grados_liber*Nnodos;
for i=1:Nnodos
    v_carga(grados_liber*i-1:grados_liber*i,1)=nodo_info{i,2}';
    indices_restringidos(1,grados_liber*i-1:grados_liber*i)=nodo_info{i,3};
end

index_free=find(indices_restringidos==1);
indice_libre=indice_total(index_free);

% indices correspondeinetes a los nudos restringidos
indice_nofree=find(indices_restringidos==0);
indice_restrin=indice_total(indice_nofree);
% matriz Sff
S_2free=M_rigidez_global(indice_libre,indice_libre);
S_free=M_rigidez_global(indice_libre,indice_restrin);
%vector de carga de los nudos libres
Af=v_carga(indice_libre');
%vector de desplazamientos para los nudos restringidos
Dr=v_dezpla(indice_restrin');
%sistema de ecuaciones para los delsplazamientos
Df=S_2free^-1*(Af-S_free*Dr)
%vector de desplazamiento de los nudos
Dj=zeros(grados_liber*Nnodos,1);
Dj(indice_libre')=Df;
Dj(indice_restrin')=v_dezpla(indice_restrin')
%acciones de nudos
reacciones=M_rigidez_global*Dj