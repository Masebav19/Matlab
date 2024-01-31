% Calculo general de celocias
clc
clear variables

% Determinamos datos generales del modelo
E=200*10^9; %[Pa] - Acero
%A=(8/100*8/100); %[m^2]
grados_liber=3; % Numero de grados de libertad por nudo

% Deteminamos los nudos
% Node2D(numero,coordenada,load,restrain)

% solicito el numero de nodos
Nnodos = input("Ingrese el numero de nodos: ");
% Pido ingresar la información de los nodos 
nodo_info = {};
for i=1:Nnodos
   msg=sprintf("Ingrese la información del nodo %d con el siguiente formato:\n{[Cooordx,Cooordy],[cargax,cargay],[restriccx,resticcy,resticcz]} 0 para restringido y 1 para libre: ",i);
   nodo_info=[nodo_info;input(msg)];
end
% Pido ingresar el numero de Barras 
Nbarra = input("Ingrese el numero de barras: ");
barra_info=[];
% Pido la información de cada barra
for i = 1:Nbarra
   msg=sprintf("Ingrese la información del elemento %d con el siguiente formato:\n[nodo_i,nodo_f]: ",i);
   barra_info=[barra_info;input(msg)];
 
end
% Ingreso del area del elemento
A=zeros(1,Nbarra);
for i = 1:Nbarra
   msg=sprintf("Ingrese el area del elemento %d: ",i);
   A(1,i)=input(msg);
 
end

% claculo de los parámetros del elemento
SJ=zeros(6,6,Nbarra);%matriz de rigidez global
T=zeros(6,6,Nbarra);%Matrz de tranformación
alpha=zeros(1,Nbarra);%angulo de cada barra
L=zeros(1,Nbarra);%angulo de cada barra
for i=1:Nbarra
   [SJ(:,:,i),alpha(i),T(:,:,i),L(i)]=calculo_barra_param(nodo_info{barra_info(i,1),1},nodo_info{barra_info(i,2),1},E,A(1,i));
end

%Ingreso de la carga distribuida
Aj_fe=zeros(6,Nbarra);
msg=sprintf("Se pide el ingreso de las cargas distribuidas en base a las tablas del 1 al 8 con el siguiente formato:\n[Y(1)/N(0),P,a,b]"); 
disp(msg);
for i=1:Nbarra
     msg=sprintf("Viga %d: ",i);
     disp(msg);
    for k=1:8
    msg=sprintf("Tabla %d: ",k);
    data=input(msg);
    Aj_fe(:,i)= Aj_fe(:,i)+calculo_carga_distr([data,k],L(i));
    end
end

%calculo de los índices para el calculo de la matriz de rigidez global de
%sistema
indices = zeros(Nbarra,6);
for i=1:Nbarra 
   indices(i,:)=[grados_liber*barra_info(i,1)-2,grados_liber*barra_info(i,1)-1,grados_liber*barra_info(i,1)...
       ,grados_liber*barra_info(i,2)-2,grados_liber*barra_info(i,2)-1,grados_liber*barra_info(i,2)];
end
%creo la matriz global del sistema
M_rigidez_global=zeros(grados_liber*Nnodos,grados_liber*Nnodos);
%con cada indice creo la matriz de cada barra dependiendo del nodo inicial
%y final
for i=1:Nbarra
    matriz_local=SJ(:,:,i);
    M_rigidez_global(indices(i,:),indices(i,:))= M_rigidez_global(indices(i,:),indices(i,:))+matriz_local;
end
v_carga_fe=zeros(grados_liber*Nnodos,1);
for i=1:Nbarra
    matriz_local=Aj_fe(:,i);
    v_carga_fe(indices(i,:),1)= v_carga_fe(indices(i,:),1)+matriz_local;
end
% vector de cargas 
% vector de desplazamiento
indices_restringidos=zeros(1,grados_liber*Nnodos);
v_carga=zeros(grados_liber*Nnodos,1);
v_carga_fe=zeros(grados_liber*Nnodos,1);
v_dezpla=zeros(grados_liber*Nnodos,1);
indice_total=1:grados_liber*Nnodos;
for i=1:Nnodos
    v_carga(grados_liber*i-1:grados_liber*i,1)=nodo_info{i,2}';
    indices_restringidos(1,grados_liber*i-2:grados_liber*i)=nodo_info{i,3};
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
Af_fe=v_carga_fe(indice_libre');
%vector de desplazamientos para los nudos restringidos
Dr=v_dezpla(indice_restrin');
%sistema de ecuaciones para los delsplazamientos
Df=S_2free^-1*(Af-Af_fe-S_free*Dr)
%vector de desplazamiento de los nudos
Dj=zeros(grados_liber*Nnodos,1);
Dj(indice_libre')=Df;
Dj(indice_restrin')=v_dezpla(indice_restrin')

Aj=zeros(grados_liber*Nnodos,1);
Aj(indice_libre')=Af_fe;
%acciones de nudos
reacciones=M_rigidez_global*Dj+Aj