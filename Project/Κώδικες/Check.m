function state=Check(x,c,V)
  
    count=0;
    a=0;
    b=0;
    e=0.002;

    for i=1:17
        if x(i)<c(i)
            count=count+1;
        end
    end
    if count==17
        a=1;
    end

    if abs(x(1)+x(2)+x(3)+x(4)-V)<e
        if abs(x(10)+x(9)-x(4))<e
            if abs(x(8)+x(7)-x(2))<e
                if abs(x(6)+x(5)-x(1))<e
                    if  abs(x(16)-x(14)-x(5))<e
                        if abs(x(17)-x(11)-x(10))<e && abs(x(13)-x(9)-x(3)-x(8)+x(11)+x(12))<e && abs(x(15)-x(13)-x(7)-x(6)+x(14))<e && abs(x(17)+x(12)+x(15)+x(16)-V)<e
                            b=1;
                        end
                    end
                end
            end
        end
    end
    if a==1 && b==1
        state=1;
    else 
        state=0;
    end
end
