%solve via cvx
%
% addpath('/Users/didemdogan/surfdrive/encoding_matrix_code_MacOS/cvx')
% cvx_setup;

%%
clear all
close all

%
addpath('/mnt/Data1/didem/encoding_matrix_code_MacOS/cvx_linux')
cvx_setup;



load('A_var_192_6x6_code3')
load('A_tr_192_6x6_code3')


n = 10;
m =5;
rng(5)
A_tr = zeros(n,n);
A_var = cell(n);

A_mat = randn(n*m);
% A_mat =5* eye(n*m);
% A_mat = diag(randn(1,m*n)');
A_mat2=A_mat*A_mat';

for i =1:size(A_tr,1)
    for j =1:size(A_tr,1)  

        A_var{i,j} = A_mat2((i-1)*m+1:i*m,(j-1)*m+1:j*m);
        A_tr(i,j)=trace(A_var{i,j});
        % A_var{i,j} = eye(n);
    end
end


for i =1:size(A_tr,1)
    for j =1:size(A_tr,1)  

        A_var{i,j} = real(A_var{i,j});
        % A_var{i,j} = eye(n);
    end
end
A_tr = real(A_tr);

for i =1:size(A_tr,1)
    for j =1:size(A_tr,1)
        A_var{i,j} =double(A_var{i,j});
    end
end


cvx_solver
cvx_begin
disp('Initializing cvx cost function');
%     disp(['Total number of iterations is ',num2str(size(A,3))]);
%     disp('Iteration no.: ');

% Variables

variable c(size(A_var,1),1);
variable C(length(c),length(c));


tic
result = cellfun(@(x,y) x.*y, num2cell(C), A_var, 'UniformOutput', false);
sumcA = sum(cat(3,result{:}),3);
toc

obj = log_det(sumcA);


maximize obj

% Constraints
disp('Initializing cvx constraints');
subject to

C == hermitian_semidefinite(size(C,1));
sum(diag(C)) == 1;

% Typical relaxation for discrete problems

disp('Finished cvx init');
cvx_end

c = sqrt(diag(C)).*sign(C(:,1));


[V,D]=eig((A_tr));
[lambda, idx] = max(diag(D));
v1=V(:,idx);

v1'*A_tr*v1;

C_new =v1*v1';

sumcA_new=zeros(size(A_var{1,1},1));

for i =1:size(A_tr,1)
    for j =1:size(A_tr,1)

        sumcA_new =sumcA_new+C_new(i,j)*double(A_var{i,j});
    end
end
log(det(sumcA_new))


% sumcA_l=zeros(size(A_var{1,1},1));
% 
% for i =1:size(A_tr,1)
%     for j =1:size(A_tr,1)
% 
%         sumcA_l =sumcA_l+C_lambda(i,j)*double(A_var{i,j});
%     end
% end
% 
% trace(inv(sumcA_l))


