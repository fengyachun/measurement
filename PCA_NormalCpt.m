function [ normal_vector_eig] = PCA_NormalCpt( point_cloud )
%UNTITLED 此处显示有关此函数的摘要
%   point_cloud为输入的点云，二维平面点或者三维空间点都可以，二维平面点点数不小于两个，三维空间点的个数不小于3个
%默认输入的点都是按行输入，所以行数为输入点个数，列数为输入点的维数，如果输入不规范，需要转置
[rows,cols] = size(point_cloud);
if(rows < cols)
    point_cloud = point_cloud';
end

[rows,cols] = size(point_cloud);%rows，输入点的个数，cols，点的维数
average_point = sum(point_cloud)/rows;
point_cloud = point_cloud-average_point;

%%特征值方法
covar_matrix = point_cloud'*point_cloud / rows; %协方差矩阵
[eig_vector,eig_value ] = eig(covar_matrix);
[sort_eig_value, sort_index] = sort(diag(eig_value),'descend');
normal_vector_eig = eig_vector(:,sort_index(end))/norm(eig_vector(:,sort_index(end)));

%%SVD方法
% [U,S,V] = svd(point_cloud);%%point_cloud是N*3的，point_cloud = U*S*V'
% normal_vector_svd = V(:,end)/norm(V(:,end));
end

