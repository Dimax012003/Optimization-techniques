clear;
a = -1;
b = 3;
l = [1e-7,1e-6,1e-5,1e-4,1e-3];
e = 0.00000002;
kol1 = [];
kol2 = [];
kol3 = [];
a_k1=[];
b_k1=[];
a_k2=[];
b_k2=[];
a_k3=[];
b_k3=[];
f1=[];
f2=[];
f3=[];
for j = 1:length(l)
    [k_1, a1, b1,a1_list,b1_list] = fibo(e, l(j), a, b, @f_x1);
    [k_2, a2, b2,a2_list,b2_list] = fibo(e, l(j), a, b, @f_x2);
    [k_3, a3, b3,a3_list,b3_list] = fibo(e, l(j), a, b, @f_x3);

    kol1 = [kol1 k_1];
    a_k1=[a_k1 a1];
    b_k1=[b_k1 b1];

    kol2=[kol2 k_2];
    a_k2=[a_k2 a2];
    b_k2=[b_k2 b2];

    kol3=[kol3 k_3];
    a_k3=[a_k3 a3];
    b_k3=[b_k3 b3];


   
        figure(1);
        plot(1:(k_1)+1,b1_list,'DisplayName',['bk ',num2str(l(j))]);
        hold on;
        plot(1:(k_1)+1,a1_list,'DisplayName',['ak ',num2str(l(j))]);
        xlabel('Αριθμός επαναλήψεων k');
        ylabel('Μεταβολή υποδιαστημάτων');
        legend show;
        grid on;

        figure(2);
        plot(1:(k_2)+1,b2_list,'DisplayName',['bk ',num2str(l(j))]);
        hold on;
        plot(1:(k_2)+1,a2_list,'DisplayName',['ak ',num2str(l(j))]);
        xlabel('Αριθμός επαναλήψεων k');
        ylabel('Μεταβολή υποδιαστημάτων');
        legend show;
        grid on;
        
        figure(3);
        plot(1:(k_3)+1,b3_list,'DisplayName',['bk ',num2str(l(j))]);
        hold on;
        plot(1:(k_3)+1,a3_list,'DisplayName',['ak ',num2str(l(j))]);
        xlabel('Αριθμός επαναλήψεων k');
        ylabel('Μεταβολή υποδιαστημάτων');
        legend show;
        grid on;

    
    f1=[f1 (a1+b1)/2];
    f2=[f2 (a2+b2)/2];
    f3=[f3 (a3+b3)/2];

end
%%Γραφική παράσταση για το τελικό διάστημα αναζήτησης.
figure();
plot(kol1, a_k1,kol1,b_k1);
title('f(x)=(x-2)^2 + x*log(x+3)');
xlabel('k τιμές επαναλήψεων');
ylabel('Μεταβολές ανω και κατω διαστήματος');
legend('Τιμές a_k','Τιμές b_k');
grid on;
saveas(gcf, 'lamda_plot_f1_ab.jpg');

figure();
plot(kol2, a_k2,kol2,b_k2);
title('f(x)=e^{-2*x} + (x-2)^2')
xlabel('k τιμές επαναλήψεων');
ylabel('Μεταβολές ανω και κατω διαστήματος');
legend('Τιμές a_k','Τιμές b_k');
grid on;
saveas(gcf, 'lamda_plot_f2_ab.jpg');

figure();
plot(kol3, a_k3,kol3,b_k3);
title('f(x)=e^x * (x^3 - 1) + (x - 1) * sin(x)');
xlabel('k τιμές επαναλήψεων');
ylabel('Μεταβολές ανω και κατω διαστήματος');
legend('Τιμές a_k','Τιμές b_k');
grid on;
saveas(gcf, 'lamda_plot_f3_ab.jpg');


figure();
plot(l, kol1);
xlabel('l values');
ylabel('Αριθμός υπολογισμών της f1(x)');
grid on;
saveas(gcf, 'lamda_plot_f1.jpg');

figure();
plot(l, kol2);
xlabel('l values');
ylabel('Αριθμός υπολογισμών της f2(x)');
grid on;
saveas(gcf, 'lamda_plot_f2.jpg');

figure();
plot(l, kol3);
xlabel('l values');
ylabel('Αριθμός υπολογισμών της f3(x)');
grid on;
saveas(gcf, 'lamda_plot_f3.jpg');

% Ορισμός της μεθόδου Fibonacci
function [k, a1, b1,a_total,b_total] = fibo(e, l, a_init, b_init, f_x)
    a_total = [];
    b_total = [];
    k = 1;
    a = a_init;
    b = b_init;
    a_total = [a_total a];
    b_total = [b_total b];
    F = [1, 1];
    
    while F(end) < (b - a) / l
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

function g = f_x1(x)
    g = (x - 2)^2 + x * log(x + 3);
end

function g = f_x2(x)
    g = exp(-2 * x) + (x - 2)^2;
end

function g = f_x3(x)
    g = exp(x) * (x^3 - 1) + (x - 1) * sin(x);
end
