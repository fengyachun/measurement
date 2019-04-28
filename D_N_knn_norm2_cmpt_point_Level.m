function [ D_N_knn_norm2 ] = D_N_knn_norm2_cmpt_point_Level( D_N,k )

[rows,cols] = size(D_N);
dim = cols/k;
D_N_knn_norm2 = zeros(rows,k);
for i = 1:1:rows
    for j = 1:1:k
        D_N_knn_norm2(i,j) = norm(D_N(i,dim*(j-1)+1:1:dim*j))*norm(D_N(i,dim*(j-1)+1:1:dim*j));
%         D_N_knn_norm2(i,j) = norm(D_N(i,dim*(j-1)+1:1:dim*j));
    end
end

end

