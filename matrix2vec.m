function N_vec = matrix2vec( N )
%UNTITLED �˴���ʾ�йش˺�����ժҪ
%   �˴���ʾ��ϸ˵��
[rows,cols] = size(N);
dim = cols;
N_vec = zeros(rows*cols,1);
for i = 1:1:rows
    N_vec(dim*(i-1)+1 : 1 : dim*i) = N(i,:)';
end
end

