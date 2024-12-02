clear;
syms x y
f=x^5*exp(-x^2-y^2);

[sol1,k1,x1,y1,g1]=newton(f,[1,-1],0.0001,x,y);
plotting(k1,x1,y1,f,g1);

[sol2,k2,x2,y2,g2]=newton(f,[-1,1],0.0001,x,y);
plotting(k2,x2,y2,f,g2);

[sol3,k3,x3,y3,g3]=newton(f,[0,0],0.0001,x,y);
plotting(k3,x3,y3,f,g3);

[sol4,k4,x4,y4,g4]=newton(f,[-1.2,-0.1],0.0001,x,y);
plotting(k4,x4,y4,f,g4);

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


function [x_k,k,x_array,y_array,g_array]=newton(f,x_0,e,x,y)
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
    while norm(grad_f)>e

        grad_f=subs(grad,{x,y},{double(x_k(1)),double(x_k(2))});

        h=subs(H,{x,y},{x_k(1),x_k(2)});

        eigenvalues = eig(h);
        isPositiveDefinite = all(eigenvalues > 0);
        
        
        if isPositiveDefinite
            dk=-inv(h)*(grad_f);

        else 
            disp('Δεν ορίζεται ο εσσιανος');
            disp(double(h));
            break;
        end

        g=0.8;

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

end

function [num_iterate, a1, b1,a,b] = dixotomos(x1, x2, l, e, a, b, f_x, str)
    fx = [];
    num_iterate = 0;
    while abs(b(length(b)) - a(length(a))) > l
        x1 = [x1 ((b(length(b)) + a(length(a))) / 2) - e];
        x2 = [x2 ((b(length(b)) + a(length(a))) / 2) + e];
        if f_x(x1(length(x1))) < f_x(x2(length(x2)))
            a = [a a(length(a))];
            b = [b x2(length(x2))];
        elseif f_x(x1(length(x1))) > f_x(x2(length(x2)))
            a = [a x1(length(x1))];
            b = [b b(length(b))];
        end 
        num_iterate = num_iterate + 2;
    end
    a1=a(length(a));
    b1=b(length(b));
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

function [num_count,a1,b1,a_list,b_list] =golden(x1,x2,a,b,l,f_x,g)
    num_count=1;
    a_list=a;
    b_list=b;
    
    
    fx1=f_x(x1);
    fx2=f_x(x2);


    while abs(b-a)>l
        if fx1 > fx2
            a=x1;
            x1=x2;
            x2=a+g*(b-a);
            fx1=fx2;
            fx2=f_x(x2);
        else
            b=x2;
            x2=x1;
            x1=a+(1-g)*(b-a);
            fx2=fx1;
            fx1=f_x(x1);
        end

        a_list = [a_list a];
        b_list = [b_list b];
        num_count = num_count + 1;
    end
    a1 = a;
    b1 = b;
end

function [k, a1, b1,a_total,b_total] = fibo(e, l, a_init, b_init, f_x)
    a_total = [];
    b_total = [];
    k = 1;
    a = a_init;
    b = b_init;
    a_total = [a_total a];
    b_total = [b_total b];
    F = [1, 1];
    
    while F(end) < (b - a) / l + 3
        F = [F, F(end) + F(end - 1)];
    end
    n = length(F); 
    x1 =a+(F(n-2)/F(n))*(b-a);
    x2 =a+(F(n-1)/F(n))*(b-a);
    fx1=f_x(x1);
    fx2=f_x(x2);

    while k <= n - 2
        if fx1 < fx2
            b = x2;
            x2 = x1;
            x1 = a + (F(n-k-2)/F(n-k))*(b-a);
            b_total = [b_total b];
            a_total = [a_total a_total(length(a_total))];
            fx2=fx1;
            fx1=f_x(x1);
            
        else
            a = x1;
            x1 = x2;
            x2 = a + (F(n - k - 1) / F(n - k)) * (b - a);
            a_total = [a_total a];
            b_total = [b_total b_total(length(b_total))];
            fx1=fx2;
            fx2=f_x(x2);
           
        end
        
        k = k + 1;
        
        if k==n-2
            x1=(a + b)/2;  
            x2=x1+e;       
            if fx1>fx2
                a = x1;
                a_total = [a_total a];
                b_total = [b_total b];
            else
                b = x2;
                b_total = [b_total b];
                a_total = [a_total a];
            end
            break;
        end
    end
    a1=a_total(length(a_total));
    b1=b_total(length(b_total));
end

function [g,m] = armijo(f, gradf, dk, a, b, g_0, xk, x, y,m)

    while double(subs(f, {x, y}, {xk(1), xk(2)})) + a * b^m * g_0 * double(dk.'* gradf) < ...
          double(subs(f, {x, y}, {xk(1) + g_0 * b^m * double(dk(1)), xk(2) + g_0 * b^m * double(dk(2))}))
        m = m + 1;
    end
    g = g_0 * b^m;
end