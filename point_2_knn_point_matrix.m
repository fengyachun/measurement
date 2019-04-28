function [ knn_point_matrix ] = point_2_knn_point_matrix( point_cloud, knn_index, k )
%UNTITLED3 此处显示有关此函数的摘要
%   此处显示详细说明
[rows,cols] = size(point_cloud);
knn_point_matrix = zeros(rows,cols*k);
for i = 1:1:rows
    tempindex = knn_index(i,:);
    tempindex = tempindex';
    tempnomal = point_cloud(tempindex,:);
    tempnormalrow = zeros(1,k*cols);
    for j = 1:1:k
        tempnormalrow(1,(cols*(j-1)+1):1:cols*j) = tempnomal(j,:);
    end
    knn_point_matrix(i,:) = tempnormalrow;
end

end

