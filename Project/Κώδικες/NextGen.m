function X=NextGen(x,c,fitfun,num_index,ability,V)
    X1=CrossOver(x,c,V,fitfun);
    X2=Mutation(x,fitfun,c,V);
    [~,indexes]=sort(ability);
    X3=x(:,indexes(91:100));
    [~,index]=Roulette(x,ability,num_index);
    X4=x(:,index);
    X=[X1,X2,X3,X4];
    X=X(:,randperm(100));
end
