function [ P_present_knn ] = func_get_knn_point( P_present, knn_index )
%UNTITLED23 �˴���ʾ�йش˺�����ժҪ
%   �˴���ʾ��ϸ˵��
[rows,cols] = size(P_present); 
[ r , k ] = size( knn_index );
P_present_knn = zeros(rows,cols*k);
for i = 1:1:rows
        tempindex = knn_index(i,:)';
        tempnomal = P_present(tempindex,:);
        tempnormalrow = zeros(1,k*cols);
        for j = 1:1:k
            tempnormalrow(1,(cols*(j-1)+1):1:cols*j) = tempnomal(j,:);
        end
        P_present_knn(i,:) = tempnormalrow;
end
end

