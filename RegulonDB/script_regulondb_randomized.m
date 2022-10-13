% script_regulondb_randomized
%
% This script generates statistics on randomized networks that "look like"
% the e-coli network.
%
% We will focus on only the nTFs part, since the other output genes could
% cloud the network.
%

load Mat/M_cell

M_out1 = M_out(:,1:nTFs);
n = full(sum(~~M_out1(:)));
N = 1e6;

I_rand = randi([1 nTFs],n,N);
J_rand = randi([1 nTFs],n,N);
n_FBL = zeros(N,1);
for ii = 1:N
	M1 = sparse(I_rand(:,ii),J_rand(:,ii),1,nTFs,nTFs,n);	
	M2 = M1';
	
	M33 = M1 & M2;
	M3 = triu(M33,1); % combine all, and take only upper triangle
	n_FBL(ii) = sum(M3(:));
	
end








