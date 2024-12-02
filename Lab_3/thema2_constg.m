clear;
syms x y
f=(1/3)*x^2+3*y^2;

[sol1,k1,x1,y1,g1]=steepest(f,[1 1],0.001,x,y,0.1);
plotting(k1,x1,y1,f,g1);
[sol2,k2,x2,y2,g2]=steepest(f,[1 1],0.001,x,y,0.3);
plotting(k2,x2,y2,f,g2);
[sol3,k3,x3,y3,g3]=steepest(f,[1 1],0.001,x,y,3);
plotting(k3,x3,y3,f,g3);
[sol4,k4,x4,y4,g4]=steepest(f,[1 1],0.001,x,y,5);
plotting(k4,x4,y4,f,g4);

function plotting(k,x1,y1,f,g)

    [x_vals, y_vals] = meshgrid(-5:0.01:5, -5:0.01:5); 


    f_numeric = matlabFunction(f);
    z_vals = f_numeric(x_vals, y_vals); 

    x0_str = num2str(x1(1), '%.2f');
    y0_str = num2str(y1(1), '%.2f');

    figure;
    hold on;
    title(['x_0 = [',num2str(x1),',',num2str(y1),']']);
    contour(x_vals, y_vals, z_vals,30, 'LineWidth', 1.5);
    scatter(x1, y1, 'filled');  


    hold on; 
    plot(x1, y1, '-o', 'Color', 'r', 'LineWidth', 1.5);
    xlabel('x');
    ylabel('y');
    title(['Ισοβαρείς καμπύλες της συνάρτησης f(x, y) για ','x_0 = [',num2str(x1(1)),',',num2str(y1(1)),']']);
    grid on;
    saveas(gcf, ['contour_plot_x0_', x0_str, '_y0_', y0_str,' γ_k' , num2str(g(1)),'.png']);

    f_array=zeros(1,length(x1));
    for i=1:length(x1)
        f_array(i)=f_numeric(x1(i),y1(i));
    end
    figure();
    plot(1:k,f_array,'LineWidth',1);
    legend('Σύκλιση f');
    xlabel('Αριθμός k επαναλήψεων');
    grid on;
    saveas(gcf, ['convergence_plot_x0_', x0_str, '_y0_', y0_str, ' γ_k' , num2str(g(1)),'.png']);


    figure();
    plot(1:k,x1,1:k,y1,'LineWidth',1);
    legend('x1','x2');
    xlabel('Αριθμός k επαναλήψεων');
    grid on;
    saveas(gcf, ['Συγκλιση_x0_', x0_str, '_y0_', y0_str, ' γ_k' , num2str(g(1)),'.png']);
end

function [x_k,k,x_array,y_array,g_array]=steepest(f,x_0,e,x,y,g)
    grad=gradient(f,[x,y]);

    k=1;

    x_k(1)=x_0(1);
    x_k(2)=x_0(2);

    grad_f=subs(grad,{x,y},{x_k(1),x_k(2)});
    x_array=[x_k(1)];
    y_array=[x_k(2)];

    g_array=[];
    m=0;
    while norm(grad_f)>e

        grad_f=subs(grad,{x,y},{double(x_k(1)),double(x_k(2))});
        dk=-grad_f;

        
        x_k(1)=x_k(1)+g*dk(1);
        x_k(2)=x_k(2)+g*dk(2);

        g_array=[g_array g];
        x_array=[x_array x_k(1)];
        y_array=[y_array x_k(2)];
        k=k+1;
        if k>200
            break;
        end
    end

end
function [i, a, b,a_list,b_list] = dixotomos_der(l, a, b, f)
    n = floor(log(l / (b - a)) / log(0.5));
    df = matlabFunction(diff(f));
    a_list=[a];
    b_list=[b];
    for i = 1:n
        x = (a + b) / 2;
        if df(x) == 0
            break;
        elseif df(x) < 0
            a = x;
            a_list=[a_list a];
            b_list=[b_list b];
        elseif df(x) > 0
            b = x;
            a_list=[a_list a];
            b_list=[b_list b];
        end
    end
end