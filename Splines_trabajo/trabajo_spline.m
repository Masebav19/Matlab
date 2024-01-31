clear
clc
img=imread("imag_spline.jpeg");
%% Spline de grado uno
%puntos de la imagen em mm
X=[20.738,39.461,43.316,70.53,97.997,102.469,102.762,102.883,108.498,122.275,133.123,169.069];
Y=149.754-[106.426,95.086,92.138,56.986,42.522,37.706,27.46,23.295,23.838,37.80,48.142,106.426];
%Calculo los términos del spline
terminos=spline_grado_uno(X,Y)
funciones ={};
k=1;
%armo las funciones
for i=1:2:length(terminos)
    funciones(k)={@(x) terminos(i)*x+terminos(i+1)}
    k=k+1;
end
%creo los vectores para graficar
x_graph=zeros(1,2*length(funciones))
y_graph=zeros(1,2*length(funciones))
k=1;
for i=1:2:(length(x_graph))
      x_graph(i:i+1)=[X(k),X(k+1)];
      y_graph(i:i+1)=[funciones{k}(X(k)),funciones{k}(X(k+1))]
      k=k+1;
end
figure(1)
subplot(2,1,1)
plot(x_graph,y_graph)
subplot(2,1,2)
imshow(img)
%% Spline de segundo grado
%puntos de la imagen em mm
X=[20.738,39.461,43.316,70.53,97.997,102.469,102.762,102.883,108.498,122.275,133.123,169.069];
Y=149.754-[106.426,95.086,92.138,56.986,42.522,37.706,27.46,23.295,23.838,37.80,48.142,106.426];
%Calculo los términos del spline
terminos=spline_grado_dos(X,Y)
funciones ={};
k=2;
%armo las funciones
funciones(1)={@(x) terminos(1).*x.^2+terminos(2).*x}
for i=3:3:length(terminos)
    funciones(k)={@(x) terminos(i).*x.^2+terminos(i+1).*x+terminos(i+2)}
    k=k+1;
end
%creo los vectores para graficar
x_graph=zeros(1,10*length(funciones))
y_graph=zeros(1,10*length(funciones))
k=1;
for i=1:10:(length(x_graph))
      x_graph(i:i+9)=X(k):(X(k+1)-X(k))/9:X(k+1);
      y_graph(i:i+9)=[funciones{k}(x_graph(i:i+9))]
      k=k+1;
end
figure(2)
subplot(2,1,1)
plot(x_graph,y_graph)
subplot(2,1,2)
imshow(img)

%% Spline de tercer grado
%puntos de la imagen em mm
X=[20.738,39.461,43.316,70.53,97.997,102.469,102.762,102.883,108.498,122.275,133.123,169.069];
Y=149.754-[106.426,95.086,92.138,56.986,42.522,37.706,27.46,23.295,23.838,37.80,48.142,106.426];
%Calculo los términos del spline
terminos=spline_grado_tres(X,Y)
funciones ={};
k=1;
%armo las funciones
for i=1:4:length(terminos)
    funciones(k)={@(x) terminos(i).*x.^3+terminos(i+1).*x.^2+terminos(i+2).*x+terminos(i+3)};
    k=k+1;
end
%creo los vectores para graficar
x_graph=zeros(1,10*length(funciones));
y_graph=zeros(1,10*length(funciones));
k=1;
for i=1:10:(length(x_graph))
      x_graph(i:i+9)=X(k):(X(k+1)-X(k))/9:X(k+1);
      y_graph(i:i+9)=[funciones{k}(x_graph(i:i+9))];
      k=k+1;
end
figure(3)
subplot(2,1,1)
plot(x_graph,y_graph)
subplot(2,1,2)
imshow(img)

