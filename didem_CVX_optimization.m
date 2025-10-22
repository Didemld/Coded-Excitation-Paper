%solve via cvx
% Optimization problem via determinant
% max lambda
% s.t. 
%

%%
clear all
close all

timer_main = tic
%
addpath('/Users/didemdogan/surfdrive/Didem PhD Tez Kodlari/Chapter7_TCI_paper_didem_lexi/encoding_matrix_code_MacOS/cvx_mac')
cvx_setup;




load('A_var_256_3x3_code2')
load('A_tr_256_3x3_code2')

% load('A_var_located_else2')
% load('A_tr_located_else_2')

for i =1:size(A_tr,1)
    for j =1:size(A_tr,1)
        A_var{i,j} =double(A_var{i,j});
    end
end

sz= size(A_var{1,1},1);
cvx_solver
cvx_begin
disp('Initializing cvx cost function');
%     disp(['Total number of iterations is ',num2str(size(A,3))]);
%     disp('Iteration no.: ');

% Variables

variable lambda;
variable c(size(A_var,1),1);
variable C(length(c),length(c));

tic
result = cellfun(@(x,y) x.*y, num2cell(C), A_var, 'UniformOutput', false);
sumcA = sum(cat(3,result{:}),3);
toc

tic
maximize lambda;

% Constraints
disp('Initializing cvx constraints');
subject to


C == hermitian_semidefinite(length(c));
% Opt constraints

% Typical relaxation for discrete problems
sumcA-lambda *eye(sz) == hermitian_semidefinite(sz);

sum(diag(C)) == 1;

disp('Finished cvx init');
cvx_end
toc

time_main =toc(timer_main)

c = sqrt(diag(C)).*sign(C(:,1));

%%
sumcA_new=zeros(size(A_var{1,1},1));

for i =1:size(A_tr,1)
    for j =1:size(A_tr,1)

        sumcA_new =sumcA_new+C_new2(i,j)*double(A_var{i,j});
    end
end
max(eig(sumcA_new))
max(eig(sumcA))




