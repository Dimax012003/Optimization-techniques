% Μπορείτε αν θέλετε αφού τρέξετε την Main και πάρετε τους πίνακες V_values
% και y_pred να τρέξετε αυτό το script.

[b0_1,b1_1]=coefficients(V_values,y_pred);
y_hat=b0_1+b1_1*V_values;

figure();
scatter(V_values,y_pred,LineWidth=1);
hold on;
xlabel('V values');
ylabel('Optimal value for fitness function');
plot(V_values,y_hat,LineWidth=1);
grid on;

n2=length(V_values);
adjR2=1-((n2-1)/(n2-2))*(sum((y_pred-y_hat).^(2)))/(sum((y_pred-mean(y_pred)).^(2)));


[B,y_hat2]=ls(y_pred,V_values,2);

[V_values, idx] = sort(V_values);
y_pred = y_pred(idx);
y_hat = y_hat(idx);
y_hat2 = y_hat2(idx);

figure();
scatter(V_values,y_pred,LineWidth=1);
hold on;
xlabel('V values');
ylabel('Optimal value for fitness function');
plot(V_values, y_hat2,LineWidth=2,Color=[0.3 0.3 0.3]);
legend('Πραγματικά δεδομένα','Πολυωνυμικό μοντέλο');
grid on;


function [b,y_hat]=ls(y,x,n)
    X=zeros(length(x),n+1);
    for i=1:n+1
        if i==1
            X(:,i)=ones(length(x),1);
        else
            X(:,i)=x.^(i-1);
        end
    end

    b=inv(X'*X)*X'*y';
    y_hat=X*b;
end

function [b0,b1]=coefficients(x,y)
    b1=(sum((x-mean(x)).*(y-mean(y))))/(sum((x-mean(x)).^(2)));
    b0=(sum(y)-b1*sum(x))/length(x);
end