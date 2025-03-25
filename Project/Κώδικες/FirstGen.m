function x=FirstGen(c,V)
    e=10^(-3);
    x=zeros(17,1);
    while abs(x(1)+x(2)+x(3)+x(4)-V)>e
        x(1)=rand*c(1);
        x(2)=rand*c(2);
        x(3)=rand*c(3);
        x(4)=rand*c(4);
    end

    k=1;
    while x(10)<=0 || k==1  
        x(9)=rand*c(9);
        x(10)=-x(9)+x(4);
        if x(10)>c(10)
            k=1;
        else 
            k=0;
        end
    end

    k=1;
    while x(8)<=0 || k==1  
        x(7)=rand*c(7);
        x(8)=-x(7)+x(2);
        if x(8)>c(8)
            k=1;
        else 
            k=0;
        end
    end

    k=1;
    while x(6)<=0 || k==1  
        x(5)=rand*c(5);
        x(6)=-x(5)+x(1);
        if x(6)>c(6)
            k=1;
        else 
            k=0;
        end
    end

    while x(15)<=0 || x(16)>c(16) || x(12)<=0 || x(13)<=0 || x(17)<=0 || x(17)>c(17) || abs(x(17)+x(12)+x(15)+x(16)-V)>e || x(13)>c(13) || x(15)>c(15) || x(16)<=0 
        x(12)=rand*c(12);
        x(14)=rand*c(14);
        x(16)=x(14)+x(5);
        x(11)=rand*c(11);
        x(17)=x(11)+x(10);
        x(13)=x(9)+x(3)+x(8)-x(11)-x(12);
        x(15)=x(13)+x(7)+x(6)-x(14);
    end
end
