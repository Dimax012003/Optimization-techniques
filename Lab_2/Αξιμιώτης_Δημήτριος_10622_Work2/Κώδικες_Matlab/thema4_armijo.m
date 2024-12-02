clear;

syms x y
f=x^5*exp(-x^2-y^2);
[sol1,k1,x1,y1,g1]=Levenberg(f,[1,-1],0.0001,x,y);
plotting(k1,x1,y1,f,g1);

[sol2,k2,x2,y2,g2]=Levenberg(f,[-1,1],0.0001,x,y);
plotting(k2,x2,y2,f,g2);

[sol3,k3,x3,y3,g3]=Levenberg(f,[0,0],0.0001,x,y);
plotting(k3,x3,y3,f,g3);

function plotting(k,x1,y1,f,g)

    [x_vals, y_vals] = meshgrid(-4:0.01:4, -4:0.01:4); 


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
    saveas(gcf, ['contour_plot_x0_', x0_str, '_y0_', y0_str, '.png']);

    f_array=zeros(1,length(x1));
    for i=1:length(x1)
        f_array(i)=f_numeric(x1(i),y1(i));
    end
    figure();
    plot(1:k,f_array);

    legend('Σύκλιση f');
    xlabel('Αριθμός k επαναλήψεων')
    grid on;
    saveas(gcf, ['convergence_plot_x0_', x0_str, '_y0_', y0_str, '.png']);

    figure();
    plot(1:(k-1),g);
    xlabel('Αριθμός k επαναλήψεων');
    ylabel('Μεταβολή βήματος γ');
    grid on;
    saveas(gcf, ['step_size_plot_x0_', x0_str, '_y0_', y0_str, '.png']);
end




function [x_k,k,x_array,y_array,g_array]=Levenberg(f,x_0,e,x,y)
    grad=gradient(f,[x,y]);
    H=hessian(f);

    k=1;
    x_k(1)=x_0(1);
    x_k(2)=x_0(2);
    grad_f=subs(grad,{x,y},{x_k(1),x_k(2)});

    x_array=[x_k(1)];
    y_array=[x_k(2)];
    g_array=[];
    m=0;
    while abs(grad_f)>e

        grad_f=subs(grad,{x,y},{double(x_k(1)),double(x_k(2))});

        h=subs(H,{x,y},{double(x_k(1)),double(x_k(2))});

        eigenvalues = eig(h);
        isPositiveDefinite = all(eigenvalues > 0);

        if isPositiveDefinite
            dk=-inv(h)*grad_f;
            
        else 
            min_e=abs(min(eigenvalues))+1;
            dk=-inv(h+min_e*[1 0;0 1])*grad_f;
        end
        
        [g,m]=armijo(f, grad_f, dk, 0.001, 0.2, 8, x_k, x, y,m);

        g_array=[g_array g];
        x_k(1)=x_k(1)+g*double(dk(1));
        x_k(2)=x_k(2)+g*double(dk(2));
        x_array=[x_array x_k(1)];
        y_array=[y_array x_k(2)];
        k=k+1;
        if k>200
            break;
        end
    end
    disp(double(subs(H,{x,y},{x_k(1),x_k(2)})));

end

function [g,m] = armijo(f, gradf, dk, a, b, g_0, xk, x, y,m)

    while double(subs(f, {x, y}, {xk(1), xk(2)})) + a * b^m * g_0 * double(dk.'* gradf) < ...
          double(subs(f, {x, y}, {xk(1) + g_0 * b^m * double(dk(1)), xk(2) + g_0 * b^m * double(dk(2))}))
        m = m + 1;
    end
    g = g_0 * b^m;
end