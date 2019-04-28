function [ res_point_cloud, res_dis, res_index ] = func_k_nearest_point_and_dis( point_cloud,centre_point,k )

[ rows , cols ] = size(point_cloud);
diff_vector =  point_cloud-centre_point;
% diff_dot = diff_vector*diff_vector';
diff_norm2 = zeros( rows , 1 );
for i = 1:1:rows
    diff_norm2( i ) = norm( diff_vector(i , :));
end

[sort_res,sort_index] = sort(diff_norm2,'ascend');
% res_point_cloud = point_cloud(sort_index(2:1:k+1),:);
% res_index = sort_index(2:1:k+1);
res_point_cloud = point_cloud(sort_index(1:1:k),:);
res_index       = sort_index(1:1:k);
res_dis         = sqrt( sort_res(1:1:k) );
end

