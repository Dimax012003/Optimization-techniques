[X,Y]=meshgrid(-5:0.1:5);
graph=X.^(5).*exp(-X.^(2)-Y.^(2));
figure(1);
surf(X,Y,graph);