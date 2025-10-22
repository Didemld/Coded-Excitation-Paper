function [c_best,C_best,trinv] = C_random_rounding_tr(C,A_var,A_tr,nof_it)



[U,S,V] = svd(C);
Wch = U*sqrt(S)*V';

best_tri = inf;
trinv = zeros(nof_it,1);

for it = 1:nof_it
    
    c_random = Wch*randn(size(C,1),1); 

    C_random = c_random*c_random';


    sumcA_new=zeros(size(A_var{1,1},1));

    for i =1:size(A_tr,1)
        for j =1:size(A_tr,1)
    
            sumcA_new =sumcA_new+C_random(i,j)*double(A_var{i,j});
        end
    end
           
    
    trinv(it) = real(trace(inv(sumcA_new)));
    
    if trinv(it) < best_tri
        best_tri = trinv(it);
        c_best = c_random;
        C_best = C_random;
    end
    
        
    
end

