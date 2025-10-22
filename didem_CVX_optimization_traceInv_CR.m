%solve via cvx
%
% addpath('/Users/didemdogan/surfdrive/encoding_matrix_code_MacOS/cvx')
% cvx_setup;

%%
clear all
close all

addpath('/mnt/Data1/didem/encoding_matrix_code_MacOS/cvx_linux')
cvx_setup;

load('A_var_256_1x1_code2')
load('A_tr_256_1x1_code2')

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

cvx_solver
cvx_begin
disp('Initializing cvx cost function');
%     disp(['Total number of iterations is ',num2str(size(A,3))]);
%     disp('Iteration no.: ');

% Variables

variable lambda;
variable c(size(A_var,1),1);
variable C(length(c),length(c));


sumcA=zeros(size(A_var{1,1},1));

variable R(size(sumcA));


result = cellfun(@(x,y) x.*y, num2cell(C), A_var, 'UniformOutput', false);
sumcA = sum(cat(3,result{:}),3);


minimize 10^10*trace(R)

% Constraints
disp('Initializing cvx constraints');
subject to

[R eye(size(R));
eye(size(R)) sumcA]== hermitian_semidefinite(2*size(sumcA,1));

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

%%
sumcA_new=zeros(size(A_var{1,1},1));

for i =1:size(A_tr,1)
    for j =1:size(A_tr,1)

        sumcA_new =sumcA_new+C_new(i,j)*double(A_var{i,j});
    end
end
trace(inv(sumcA_new))
trace(inv(sumcA))


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


