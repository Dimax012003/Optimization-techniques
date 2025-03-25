function [p,select_index]=Roulette(x,fitfun,num_index)
    n=length(x);
    p=zeros(100,1);
    total=sum(fitfun);
    for i=1:n
        p(i)=fitfun(i)/total;
    end

    CDF=cumsum(p);
    select_index=zeros(num_index,1);
    for i=1:length(select_index)
        r = rand;  
        select_index(i) = find(CDF >= r, 1);
    end
end
