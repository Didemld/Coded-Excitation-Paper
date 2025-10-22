%solve via cvx
% Optimization problem via determinant
% max lambda
% s.t. 
%




%
% addpath('/Users/didemdogan/surfdrive/encoding_matrix_code_MacOS/cvx')
% cvx_setup;


%%
clear all
close all
clc

% % 
load('A_var_256_1x1_code2')
load('A_tr_256_1x1_code2')



maxIter = 300;
beta = 5*10^6;

c_init = rand(size(A_var,1),1);
c_init = c_init/norm(c_init);

c_new = c_init;
for iter = 1: maxIter    
    C = c_new*c_new'; 
    result = cellfun(@(x,y) x.*y, num2cell(C), A_var, 'UniformOutput', false);
    sumcA = sum(cat(3,result{:}),3);
    
    

    for i = 1:size(A_tr,1)
    
        middle_term = 0;
        for j = 1:size(A_tr,1)  
    
            middle_term =middle_term+c_new(j)'*(A_var{i,j});
        end
        sumcA_derivative = -trace(inv(sumcA)* middle_term *inv(sumcA));
    
        c_new(i) = c_new(i) - beta * sumcA_derivative; 

    end 
    values(iter) =real(trace(inv(sumcA)));
    plot((values))
    hold on
    drawnow
%     if iter >1 && values(iter) > values(iter-1)
%         break;
%     end
end


c_new = c_new/norm(c_new);

objective = trace(inv(sumcA))




