function X=CrossOver(x,c,V,fitfun)
n=75;
    X=zeros(17,n);
    for i=1:n
        x_star=zeros(17,1);

        while Check(x_star,c,V)~=1 
            parent1=randi(100);                    %Roulette(x,fitfun,1);
            parent2=randi(100);                    %Roulette(x,fitfun,1);
            while parent2==parent1
                parent2=randi(100);                           %Roulette(x,fitfun,1);
            end
            x_star=(x(:,parent1)+x(:,parent2))/2;
        end
        X(:,i)=x_star;
    end
end
