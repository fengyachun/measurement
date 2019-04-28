function [ normal_P_present, knn_index, knn_P_present_dis ] = func_knn_normal_and_index_dis_cmpt( P_present, k )

[rows, cols] = size( P_present );
normal_P_present     = zeros( rows , cols );
knn_index            = zeros( rows , k );
knn_P_present_dis    = zeros( rows , k );

for i = 1:1:rows

        [knn_point_i,res_dis, knn_index_i] = func_k_nearest_point_and_dis( P_present,P_present(i,:),k);
        knn_index(i,:) = knn_index_i';%%%
        knn_P_present_dis(i,:) = res_dis';
        %%%%%
        normal_eig_i = PCA_NormalCpt( knn_point_i ); 
        normal_P_present(i,:) = normal_eig_i';  
end

end

