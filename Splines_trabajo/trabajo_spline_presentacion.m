clear
clc
figure('Position',[500,200,400,300],'Name','Splines')
uicontrol('units','normalized','position',[0.2,0.2,0.5,0.05],...
    'style','popup','string',{'Spline Lineal','Spline cuadrático','Spline cúbico'},...
    'value',1,'callback',@changetype);
uicontrol('units','normalized','position',[0.01,0.9,0.8,0.09],...
    'style','text','String','Universidad Politécnica Salesiana','fontname','Calibri Ligth','FontSize',14);
uicontrol('units','normalized','position',[-0.05,0.8,0.8,0.09],...
    'style','text','String','Análisis Numérico y Programación','fontname','Calibri Ligth','FontSize',12);
uicontrol('units','normalized','position',[-0.01,0.7,0.5,0.09],...
    'style','text','String','M. Sc. Paulina Morillo','fontname','Calibri Ligth','FontSize',12);
uicontrol('units','normalized','position',[-0.14,0.6,0.5,0.09],...
    'style','text','String','Grupo1','fontname','Calibri Ligth','FontSize',12);
uicontrol('units','normalized','position',[-0.05,0.5,0.4,0.09],...
    'style','text','String','Período 62','fontname','Calibri Ligth','FontSize',12);
uicontrol('units','normalized','position',[-0.26,0.4,0.8,0.09],...
    'style','text','String','Año 2023','fontname','Calibri Ligth','FontSize',12);
logo=imread("logo.png");
ax=axes('Position',[0.5,0.35,0.4,0.4]);
imshow(logo);
function changetype(hObj,~)
    img=imread("imag_spline.jpeg");
    %puntos de la imagen em mm
    X=[20.738,39.461,43.316,70.53,97.997,102.469,102.762,102.883,108.498,122.275,133.123,169.069];
    Y=149.754-[106.426,95.086,92.138,56.986,42.522,37.706,27.46,23.295,23.838,37.80,48.142,106.426];
          switch get(hObj,'value')
              case 1    
                %Calculo los términos del spline
                terminos=spline_grado_uno(X,Y);
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
                figure('Name','Spline lineal')
                subplot(2,1,1)
                plot(x_graph,y_graph,"-k",'LineWidth',3)
                hold on
                plot(X,Y,'o')
                hold off
                subplot(2,1,2)
                imshow(img)
                %terminos = [a1,b1,a2,b2]
                a=[];
                k=1;
                for i=1:2:length(terminos)
                a(k)=terminos(i);
                k=k+1;
                end
                k=1;
                b=[];
                for i=2:2:length(terminos)
                b(k)=terminos(i);
                k=k+1;
                end
                figure('Name','Tabla de valores',Position=[100,100,280,200])
                uitable("Data",[a',b'],'ColumnName',{'Termino a','Termino b'},'Units','normalized',Position=[0,0.4,0.72,0.56])
              case 2
                %Calculo los términos del spline
                terminos=spline_grado_dos(X,Y);
                funciones ={};
                k=2;
                %armo las funciones
                funciones(1)={@(x) terminos(1).*x.^2+terminos(2).*x};
                for i=3:3:length(terminos)
                    funciones(k)={@(x) terminos(i).*x.^2+terminos(i+1).*x+terminos(i+2)};
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
                figure('Name','Spline Cuadrático');
                subplot(2,1,1)
                plot(x_graph,y_graph,"-k",'LineWidth',3)
                hold on
                plot(X,Y,'o')
                hold off
                title("Relieve con Spline grado uno")
                subplot(2,1,2)
                imshow(img)

                a=[];
                k=1;
                for i=1:3:length(terminos)
                a(k)=terminos(i);
                k=k+1;
                end
                
                k=1;
                b=[];
                for i=2:3:length(terminos)
                b(k)=terminos(i);
                k=k+1;
                end

                k=2;
                c=0;
                for i=3:3:length(terminos)
                c(k)=terminos(i);
                k=k+1;
                end
    
                figure('Name','Tabla de valores',Position=[100,100,290,200]);
                uitable("Data",[a',b',c'],'ColumnName',{'Termino a','Termino b','Termino c'},'Units','normalized',Position=[0,0.4,0.98,0.57]);
              case 3    
                %Calculo los términos del spline
                terminos=spline_grado_tres(X,Y);
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
                figure('Name','Spline Cúbico')
                subplot(2,1,1)
                plot(x_graph,y_graph,"-k",'LineWidth',3)
                hold on
                plot(X,Y,'o')
                hold off
                title("Relieve con Spline grado tres")
                title("Relieve con Spline grado dos")
                subplot(2,1,2)
                imshow(img)

                a=[];
                k=1;
                for i=1:4:length(terminos)
                a(k)=terminos(i);
                k=k+1;
                end
                
                k=1;
                b=[];
                for i=2:4:length(terminos)
                b(k)=terminos(i);
                k=k+1;
                end

                k=1;
                c=[];
                for i=3:4:length(terminos)
                c(k)=terminos(i);
                k=k+1;
                end

                k=1;
                d=[];
                for i=4:4:length(terminos)
                d(k)=terminos(i);
                k=k+1;
                end
    
                figure('Name','Tabla de valores',Position=[100,100,400,200]);
                uitable("Data",[a',b',c',d'],'ColumnName',{'Termino a','Termino b','Termino c','Termino d'},'Units','normalized',Position=[0,0.4,0.92,0.56]);
          end
end