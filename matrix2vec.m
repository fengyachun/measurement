function N_vec = matrix2vec( N )
%UNTITLED 此处显示有关此函数的摘要
%   此处显示详细说明
[rows,cols] = size(N);
dim = cols;
N_vec = zeros(rows*cols,1);
for i = 1:1:rows
    N_vec(dim*(i-1)+1 : 1 : dim*i) = N(i,:)';
end
end

