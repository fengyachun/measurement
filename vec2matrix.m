function [ N_matrix ] = vec2matrix( N_vec, dim )
%UNTITLED2 �˴���ʾ�йش˺�����ժҪ
%   �˴���ʾ��ϸ˵��
[rows,cols] = size(N_vec);
nums = rows/dim;
N_matrix = zeros(nums,dim);
for i = 1:1:nums
    N_matrix( i,: ) = N_vec( dim*(i-1)+1 : 1 : dim*i)';  
end
end

