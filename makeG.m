
function G  = makeG(g_receive,g_transmit,O_c,p,L,numelements,M)
% tic
G=(zeros(numelements,M,L,'single'));
for l =1:L
%     G(:,:,l) =g_receive(:,:,l)*diag(g_transmit(:,:,l)*O_1*p(l,:).'); %31.8s
    G(:,:,l) = bsxfun (@times,g_receive(:,:,l), (g_transmit(:,:,l)*O_c*p(l,:).').'); %1.18s
end
% toc
end