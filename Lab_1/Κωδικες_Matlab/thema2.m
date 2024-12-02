clear;
clc;
g=0.618;
a=-1;
b=3;
x1=a+(1-g)*(b-a);
x2=a+g*(b-a);
l=[1e-7,1e-6,1e-5,1e-4,1e-3];
kol1=[];
kol2=[];
kol3=[];
a_k1=[];
a_k2=[];
b_k1=[];
b_k2=[];
a_k3=[];
b_k3=[];
f1=[];
f2=[];
f3=[];

for j = 1:length(l)
    [k_1,a1,b1,a1_list,b1_list]=golden(x1,x2,a,b,l(j),@f_x1,g);
    [k_2,a2,b2,a2_list,b2_list]=golden(x1,x2,a,b,l(j),@f_x2,g);
    [k_3,a3,b3,a3_list,b3_list]=golden(x1,x2,a,b,l(j),@f_x3,g);

    kol1 = [kol1 k_1];
    a_k1 = [a_k1 a1];
    b_k1 = [b_k1 b1];

    kol2 = [kol2 k_2];
    a_k2 = [a_k2 a2];
    b_k2 = [b_k2 b2];

    kol3 = [kol3 k_3];
    a_k3 = [a_k3 a3];
    b_k3 = [b_k3 b3];


    figure(1);
    plot(1:k_1, b1_list, 'DisplayName', ['bk ', num2str(l(j))]);
    hold on;
    plot(1:k_1, a1_list, 'DisplayName', ['ak ', num2str(l(j))]);
    xlabel('Αριθμός επαναλήψεων k');
    ylabel('Μεταβολή υποδιαστημάτων');
    legend show;
    grid on;

    figure(2);
    plot(1:k_2, b2_list, 'DisplayName', ['bk ', num2str(l(j))]);
    hold on;
    plot(1:k_2, a2_list, 'DisplayName', ['ak ', num2str(l(j))]);
    xlabel('Αριθμός επαναλήψεων k');
    ylabel('Μεταβολή υποδιαστημάτων');
    legend show;
    grid on;

    figure(3);
    plot(1:k_3, b3_list, 'DisplayName', ['bk ', num2str(l(j))]);
    hold on;
    plot(1:k_3, a3_list, 'DisplayName', ['ak ', num2str(l(j))]);
    xlabel('Αριθμός επαναλήψεων k');
    ylabel('Μεταβολή υποδιαστημάτων');
    legend show;
    grid on;

    f1 = [f1 (a1+b1)/2];
    f2 = [f2 (a2+b2)/2];
    f3 = [f3 (a3+b3)/2];
end

%%Γραφική παράσταση για το τελικό διάστημα αναζήτησης.
figure();
plot(kol1, a_k1, kol1, b_k1);
title('f(x)=(x-2)^2 + x*log(x+3)');
xlabel('k τιμές επαναλήψεων');
ylabel('Μεταβολές ανω και κατω διαστήματος');
legend('Τιμές a_k', 'Τιμές b_k');
grid on;
saveas(gcf, 'lamda_plot_f1_ab.jpg');

figure();
plot(kol2, a_k2, kol2, b_k2);
title('f(x)=e^{-2*x} + (x-2)^2');
xlabel('k τιμές επαναλήψεων');
ylabel('Μεταβολές ανω και κατω διαστήματος');
legend('Τιμές a_k', 'Τιμές b_k');
grid on;
saveas(gcf, 'lamda_plot_f2_ab.jpg');

figure();
plot(kol3, a_k3, kol3, b_k3);
title('f(x)=e^x * (x^3 - 1) + (x - 1) * sin(x)');
xlabel('k τιμές επαναλήψεων');
ylabel('Μεταβολές ανω και κατω διαστήματος');
legend('Τιμές a_k', 'Τιμές b_k');
grid on;
saveas(gcf, 'lamda_plot_f3_ab.jpg');

%%Γραφική παράσταση αριθμών επαναλήψεων συναρτήσει του l

figure();
plot(l, kol1);
xlabel('l τιμές');
ylabel('Αριθμός υπολογισμών της f1(x)');
grid on;
saveas(gcf, 'lamda_plot_f1.jpg');

figure();
plot(l, kol2);
xlabel('l τιμές');
ylabel('Αριθμός υπολογισμών της f2(x)');
grid on;
saveas(gcf, 'lamda_plot_f2.jpg');

figure();
plot(l, kol3);
xlabel('l τιμές');
ylabel('Αριθμός υπολογισμών της f3(x)');
grid on;
saveas(gcf, 'lamda_plot_f3.jpg');


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

function g=f_x1(x)
    g=(x-2)^2+x*log(x + 3);
end

function g=f_x2(x)
    g=exp(-2*x)+(x-2)^2;
end

function g = f_x3(x)
    g=exp(x)*(x^3-1)+(x-1)*sin(x);
end
