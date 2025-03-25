clear;

% αυτό το κομμάτι δεν είναι απαραίτητο εάν θέλετε κάνετε comment μέχρι το
% τέλος της loop.Χρησιμοποίησα αυτό το κομμάτι του κώδικα από 6-14 μονο στο
% script.m
samples=100;
y_pred=zeros(1,samples);
V_values=zeros(1,samples);
for i=1:samples
    V=100+30*(rand-0.5);
    V_values(i)=V;
    y_pred(i)=main(V,1);
end
% μέχρι εδώ δλδ.

main(V,0);
main(100,0);

function y_hat=main(V,hide)
c = [54.13, 21.56, 34.08, 49.19, 33.03, 21.84, 29.96, 24.87, 47.24, 33.97, 26.89, 32.76, 39.98, 37.12, 58.83, 61.65, 59.73]; 

X = zeros(17, 100); 
for i = 1:100
    X(:, i) = FirstGen(c,V);  
end
y = zeros(1, 100); 
g = zeros(1, 100); 
for i = 1:100
    [y(i), g(i)] = FitnessFun(X(:, i), c); 
end
M=1000;
y1=zeros(1,M);
means=zeros(1,M);
for j = 1:M
    X = NextGen(X, c, y, 100, g,V);  
    for i = 1:100
        [y(i), g(i)] = FitnessFun(X(:, i), c); 
    end
    %disp(j);
    y1(j)=min(y);
    means(j)=mean(y);
end
y_hat=min(y);
disp(['Minimum fitness value: ', num2str(min(y)),' V:',num2str(V)]);
if (hide~=1)
    figure();
    plot(1:M,y1,LineWidth=1);
    title(['V=',num2str(V)]);
    grid on;

    figure();
    plot(1:M,means,LineWidth=1);
    title(['V=',num2str(V)]);
    ylabel('Average Fitness value through generation');
    xlabel('Generations');
    grid on;
end
end