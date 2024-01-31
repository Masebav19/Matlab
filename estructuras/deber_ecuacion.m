clear
clc
value=input("Ingrese los valores de a,b,c y d enrre corchetes y separados por comas: ");
syms x
y=value(1)*x^3+value(2)*x^2+value(3)*x+value(4)
yp=diff(y)
ypp=diff(yp)
eq= yp==0;
t_eval=double(solve(eq));
t_graph=[];
x_eval=[];
y_eval=[];
n=0;
for K=1:length(t_eval)
    if isreal(t_eval(K))
        t_graph=[t_graph,t_eval(K)];
        aux=ypp;
        x_eval=[x_eval,double(subs(aux,t_eval(K)))];
        if x_eval(K) > 0 
            aux2=y;
            y_eval=[y_eval,double(subs(aux2,t_eval(K)))];
           msg=sprintf("Valor mínimo: %.2f",double(subs(aux2,t_eval(K))));
           disp(msg)
        else
           aux2=y;
           y_eval=[y_eval,double(subs(aux2,t_eval(K)))];
           msg=sprintf("Valor máximo: %.2f",double(subs(aux2,t_eval(K))));
           disp(msg) 
        end
    else
        n=n+1;
    end
    if n==length(t_eval)
       disp("La función no tiene valores mínimos ni máximos") 
    end
end
if length(t_graph)>0
fplot(y,t_graph+[-1,1])
title("Gráfica de la función ax^3+bx^2+cx+d");
ylabel("y")
xlabel("x")
hold on 
for j=1:length(t_graph)
    msg=sprintf("punto %d:(%.2f,%.2f)",j,t_graph(j),y_eval(j));
    text(t_graph(j),y_eval(j),msg);
end
hold off
else
  fplot(y,[-2,2])
title("Gráfica de la función ax^3+bx^2+cx+d");
ylabel("y")
xlabel("x")  
end
