function [Ttotal,Ability]=FitnessFun(x,c)
    a=[1.25,1.25,1.25,1.25,1.25,1.5,1.5,1.5,1.5,1.5,1,1,1,1,1,1,1];
    t=0*ones(17,1);
    Ttotal=0;
    for i=1:17
        T=t(i)+a(i)*x(i)/(1-x(i)/c(i));
        Ttotal=Ttotal+T;
    end
    Ability=1/(Ttotal);
end
