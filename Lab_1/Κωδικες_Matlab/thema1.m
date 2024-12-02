clear;
clc;
x1 = [];
x2 = [];
l = 0.01;
e = 0.000005:0.0002:0.0049;
a = -1;
b = 3;
k1 = [];
k2 = [];
k3 = [];


for i = [1:length(e)]
    fx1 = dixotomos(x1, x2, l, e(i), a, b, @f_x1, '(x-2)^2 + x*log(x+3)');
    fx2 = dixotomos(x1, x2, l, e(i), a, b, @f_x2, 'e^{-2*x} + (x-2)^2');
    fx3 = dixotomos(x1, x2, l, e(i), a, b, @f_x3, 'e^x * (x^3 - 1) + (x - 1) * sin(x)');
    
    k1 = [k1 fx1];
    k2 = [k2 fx2];
    k3 = [k3 fx3];
end

figure();
plot(e, k1);
xlabel('e values');
ylabel('Αριθμός υπολογισμών της f1(x)');
grid on;
saveas(gcf, 'plot_examplef1.jpg');

figure();
plot(e, k2);
xlabel('e values');
ylabel('Αριθμός υπολογισμών της f2(x)');
grid on;
saveas(gcf, 'plot_examplef2.jpg');

figure();
plot(e, k3);
xlabel('e values');
ylabel('Αριθμός υπολογισμών της f3(x)');
grid on;
saveas(gcf, 'plot_examplef3.jpg'); 

epsilon = 0.001;
lamda = [0.0021:0.001:0.01];

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
for j = [1:length(lamda)]
    [k_1, a1, b1,a1_list,b1_list] = dixotomos(x1, x2, lamda(j), epsilon, a, b, @f_x1, '(x-2)^2 + x*log(x+3)');
    [k_2, a2, b2,a2_list,b2_list] = dixotomos(x1, x2, lamda(j), epsilon, a, b, @f_x2, 'e^{-2*x} + (x-2)^2');
    [k_3, a3, b3,a3_list,b3_list] = dixotomos(x1, x2, lamda(j), epsilon, a, b, @f_x3, 'e^x * (x^3 - 1) + (x - 1) * sin(x)');

    kol1 = [kol1 k_1];
    a_k1=[a_k1 a1];
    b_k1=[b_k1 b1];

    kol2=[kol2 k_2];
    a_k2=[a_k2 a2];
    b_k2=[b_k2 b2];

    kol3=[kol3 k_3];
    a_k3=[a_k3 a3];
    b_k3=[b_k3 b3];
    
   
        figure(4);
        plot(1:(k_1/2)+1,b1_list,'DisplayName',['bk ',num2str(lamda(j))]);
        hold on;
        plot(1:(k_1/2)+1,a1_list,'DisplayName',['ak ',num2str(lamda(j))]);
        xlabel('Αριθμός επαναλήψεων k');
        ylabel('Μεταβολή υποδιαστημάτων');
        legend show;
        grid on;

        figure(5);
        plot(1:(k_2/2)+1,b2_list,'DisplayName',['bk ',num2str(lamda(j))]);
        hold on;
        plot(1:(k_2/2)+1,a2_list,'DisplayName',['ak ',num2str(lamda(j))]);
        xlabel('Αριθμός επαναλήψεων k');
        ylabel('Μεταβολή υποδιαστημάτων');
        legend show;
        grid on;
        
        figure(6);
        plot(1:(k_3/2)+1,b3_list,'DisplayName',['bk ',num2str(lamda(j))]);
        hold on;
        plot(1:(k_3/2)+1,a3_list,'DisplayName',['ak ',num2str(lamda(j))]);
        xlabel('Αριθμός επαναλήψεων k');
        ylabel('Μεταβολή υποδιαστημάτων');
        legend show;
        grid on;

    
    % υπολογισμοί ελαχίστου για κάθε l τιμή
    f1=[f1 (a1+b1)/2];
    f2=[f2 (a2+b2)/2];
    f3=[f3 (a3+b3)/2];
end

%%Γραφική παράσταση για το τελικό διάστημα αναζήτησης.
figure();
plot(kol1/2, a_k1,kol1/2,b_k1);
title('f(x)=(x-2)^2 + x*log(x+3)');
xlabel('k τιμές επαναλήψεων');
ylabel('Μεταβολές ανω και κατω διαστήματος');
legend('Τιμές a_k','Τιμές b_k');
grid on;
saveas(gcf, 'lamda_plot_f1_ab.jpg');

figure();
plot(kol2/2, a_k2,kol2/2,b_k2);
title('f(x)=e^{-2*x} + (x-2)^2')
xlabel('k τιμές επαναλήψεων');
ylabel('Μεταβολές ανω και κατω διαστήματος');
legend('Τιμές a_k','Τιμές b_k');
grid on;
saveas(gcf, 'lamda_plot_f2_ab.jpg');

figure();
plot(kol3/2, a_k3,kol3/2,b_k3);
title('f(x)=e^x * (x^3 - 1) + (x - 1) * sin(x)');
xlabel('k τιμές επαναλήψεων');
ylabel('Μεταβολές ανω και κατω διαστήματος');
legend('Τιμές a_k','Τιμές b_k');
grid on;
saveas(gcf, 'lamda_plot_f3_ab.jpg');






figure();
plot(lamda, kol1);
xlabel('l τιμές');
ylabel('Αριθμός υπολογισμών της f1(x)');
grid on;
saveas(gcf, 'lamda_plot_f1.jpg');

figure();
plot(lamda, kol2);
xlabel('l τιμές');
ylabel('Αριθμός υπολογισμών της f2(x)');
grid on;
saveas(gcf, 'lamda_plot_f2.jpg');

figure();
plot(lamda, kol3);
xlabel('l τιμές');
ylabel('Αριθμός υπολογισμών της f3(x)');
grid on;
saveas(gcf, 'lamda_plot_f3.jpg');

% Ορισμός συναρτήσεων
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

function g = f_x1(x)
    g = (x - 2)^2 + x * log(x + 3);
end

function g = f_x2(x)
    g = exp(-2 * x) + (x - 2)^2;
end

function g = f_x3(x)
    g = exp(x) * (x^3 - 1) + (x - 1) * sin(x);
end
