function [ normal_vector_eig] = PCA_NormalCpt( point_cloud )
%UNTITLED �˴���ʾ�йش˺�����ժҪ
%   point_cloudΪ����ĵ��ƣ���άƽ��������ά�ռ�㶼���ԣ���άƽ��������С����������ά�ռ��ĸ�����С��3��
%Ĭ������ĵ㶼�ǰ������룬��������Ϊ��������������Ϊ������ά����������벻�淶����Ҫת��
[rows,cols] = size(point_cloud);
if(rows < cols)
    point_cloud = point_cloud';
end

[rows,cols] = size(point_cloud);%rows�������ĸ�����cols�����ά��
average_point = sum(point_cloud)/rows;
point_cloud = point_cloud-average_point;

%%����ֵ����
covar_matrix = point_cloud'*point_cloud / rows; %Э�������
[eig_vector,eig_value ] = eig(covar_matrix);
[sort_eig_value, sort_index] = sort(diag(eig_value),'descend');
normal_vector_eig = eig_vector(:,sort_index(end))/norm(eig_vector(:,sort_index(end)));

%%SVD����
% [U,S,V] = svd(point_cloud);%%point_cloud��N*3�ģ�point_cloud = U*S*V'
% normal_vector_svd = V(:,end)/norm(V(:,end));
end

