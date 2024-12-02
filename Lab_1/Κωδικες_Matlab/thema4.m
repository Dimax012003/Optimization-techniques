clear;
a = -1;
b = 3;
syms f1(x) f2(x) f3(x)
f1(x) = (x - 2)^(2) + x * log(x + 3);
f2(x) = exp(-2 * x) + (x - 2)^2;
f3(x) = exp(x) * (x^3 - 1) + (x - 1) * sin(x);
l = [1e-7,1e-6,1e-5,1e-4,1e-3];
kol1 = [];
kol2 = [];
kol3 = [];
a_k1=[];
b_k1=[];
a_k2=[];
b_k2=[];
a_k3=[];
b_k3=[];
f_1=[];
f_2=[];
f_3=[];
for j = 1:length(l)
    [k_1, a1, b1,a1_list,b1_list] = dixotomos_der(l(j), a, b, f1);
    [k_2, a2, b2,a2_list,b2_list] = dixotomos_der(l(j), a, b, f2);
    [k_3, a3, b3,a3_list,b3_list] = dixotomos_der(l(j), a, b, f3);

    kol1 = [kol1 k_1];
    a_k1=[a_k1 a1];
    b_k1=[b_k1 b1];

    kol2=[kol2 k_2];
    a_k2=[a_k2 a2];
    b_k2=[b_k2 b2];

    kol3=[kol3 k_3];
    a_k3=[a_k3 a3];
    b_k3=[b_k3 b3];
    % Σχεδίαση και αποθήκευση γραφήματος f1
    

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

    
    
    f_1=[f_1 (a1+b1)/2];
    f_2=[f_2 (a2+b2)/2];
    f_3=[f_3 (a3+b3)/2];

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

% Συνάρτηση διχοτόμησης με παραγώγους
function [i, a, b,a_list,b_list] = dixotomos_der(l, a, b, f)
    n = floor(log(l / (b - a)) / log(0.5) );
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
