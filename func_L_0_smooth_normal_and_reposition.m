function [P_resosition, N_smooth] = func_L_0_smooth_normal_and_reposition( P_init, k, step_max)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% <<<<<<<<  main initial >>>>>>>>>> %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
P_present = P_init;
[rows,cols] = size(P_present);  
%%%%%%%%%%%%parameter setting of normal vector cluster 
beta_max = 1e5;
kapa_normal = 1.2;
yita = 0.3;
beta_init = 1;  % D_N_ik_j ? yita/beta;  Update: beta = kapa_normal*beta

%%%%%%%%%%%%parameter setting of point reposition 
lamda_init = 1;
delta = 10;
kapa_point = 1.5;
lamda_max = 1e4;  % W_ik_j^2 ? delta/lamda;  Update: lamda = kapa_point * lamda

for step = 1:1:step_max
    beta  = beta_init * 1.1^(step-1); %%yita/beta:  beta = kapa_normal * beta
    lamda = lamda_init * 1.2^(step-1); %%delta/lamda:  lamda = kapa_point * lamda    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% <<<<<<<<  main step0 -  calculate knn & PCA normal >>>>>>>>>>  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%k½üÁÚ
    [ knn_index, normal_point_cloud ] = func_knn_normal_and_index_cmpt( P_present, k );
    %%%%%%%%%%%%%%%
    knn_PCA_normal= func_get_knn_point( normal_point_cloud, knn_index );    
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% <<<<<<<<  main step1 -- normal smooth>>>>>>>>>>  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
    D_N = zeros(rows,cols*k);
    for i = 1:1:k
        disp(i);
        D_N(:,(cols*(i-1)+1):1:(cols*i)) = normal_point_cloud - knn_PCA_normal(:,(cols*(i-1)+1):1:(cols*i));%% D(N)_ik+j=N_i-N_M(i,j)
    end
   
    % N_vec_j_knn = C_j * N_vec;   
    C_0 = zeros(rows*cols,rows*cols);
    C_k = repmat(C_0,[1,1,k]);
    for j = 1:1:k
        for i = 1:1:rows
            C_k((cols*(i-1)+1) : 1 : (cols*i) , (cols * (knn_index(i,j) - 1) + 1) : 1 : cols * knn_index(i,j), j )  = eye(cols);
        end
    end
     % A_j = I - C_j;
    A_j_T_A_j_K = repmat(C_0,[1,1,k]); 
    A_N_k = repmat(C_0,[1,1,k]);
    I = eye(rows*cols);
    for j = 1:1:k
        A_N_k(:,:,j) = ( I - C_k(:,:,j) );
        A_j_T_A_j_K(:,:,j) = ( I - C_k(:,:,j) )' * ( I - C_k(:,:,j));%%A_T_j = (I-C_k(:,:,j))
    end
    A_N_K_sum = C_0;
    for i = 1:1:rows*cols
        for j = 1:1:rows*cols
            A_N_K_sum(i,j) = sum(A_j_T_A_j_K(i,j,:));
        end
    end
    N_init = matrix2vec(normal_point_cloud);
    D_N_new = D_N;
    %%%%normal smooth 
    while beta < beta_max    
        %disp(beta);
        %update theta£¬objective function£ºmin_theta { beta * || D_N_new-theta ||_2 + yita * || theta ||_0 }
        %D_N_knn_norm2 = D_N_knn_norm2_cmpt_row_vector_level( D_N_new , k );%sigma_j_k|D(N)_ik+j| page7
        D_N_knn_norm2 = D_N_knn_norm2_cmpt_point_Level( D_N_new , k );
        theta = D_N_new;
        ratio = yita / beta;
        disp(ratio);
        [r,c] = find( D_N_knn_norm2 < ratio );
        for i = 1:1:length(r)
            theta(r(i), cols*(c(i)-1) + 1 : 1 : cols*c(i) ) = 0;
        end    
        %update N_new£¬objective function£ºmin_N_new {  || N - N_new ||_2 + beta * || D_N_new - theta ||_2 }
       
        %%N_new = (I-beta * sigma_j_k ( A_T_j * A_j ) )_-1 * ( N_init + beta * sigma_j_k( A_T_j * theta_:k_j ) )
        A_T_j_theta_k_j = zeros(rows*cols , 1);
        for j = 1:1:k
            A_T_j_theta_k_j = A_T_j_theta_k_j + A_N_k(:,:,j)' * matrix2vec( theta( : , cols*(j-1)+1:1:cols*j ));
        end 
        b_vec = N_init + beta *A_T_j_theta_k_j;
        A_inverse = inv(I + beta * A_N_K_sum);    
        % N_vec_new = (I + beta*A_K_sum)\b_vec;
        N_vec_new = A_inverse * b_vec;
        N_matrix =vec2matrix( N_vec_new, cols );
        N_matrix =  N_to_unit_N( N_matrix,1 );%%
        knn_point_matrix = point_2_knn_point_matrix( N_matrix, knn_index, k );      
        for i = 1:1:k        
            D_N_new(:,(cols*(i-1)+1):1:(cols*i)) = N_matrix - knn_point_matrix(:,(cols*(i-1)+1):1:(cols*i));%% D(N)_ik+j=N_i-N_M(i,j)
        end    
        %update beta
        beta = kapa_normal*beta;    
    end
   
    N_smooth = N_matrix;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% <<<<<<<<  main step2 >>>>>>>>>>  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%reposition£»    
    %
    N_vec = matrix2vec( N_matrix );
    N_vec_diag = diag(N_vec, 0); %
    N_nn_diag = N_nn_diag_cmpt( N_matrix );%%%Ni' * Ni;£¨rows*cols£©x (rows*cols);
    N_vec_diag_T_mul_N_vec_diag = N_vec_diag*N_vec_diag;

    P_init_vec = matrix2vec( P_init );
    P_present_vec = matrix2vec( P_present );
    %C_k;
    A_P_k = repmat(C_0,[1,1,k]);
    
    for j = 1:1:k
        A_P_k( : , : , j ) = N_vec_diag -  diag( C_k( : , : , j ) *  N_vec , 0 ) * C_k( : , : , j );
    end;
    %%%( S1 + lamda * S2 ) * alpha = (b1 + b2) 
    S1 = N_vec_diag_T_mul_N_vec_diag;
    S2 = zeros( rows * cols , rows * cols);
    for j = 1 : 1 : k
        S2 = S2 + A_P_k( : , : , j )' * N_nn_diag * A_P_k( : , : , j );
    end
    I_P = eye(rows*cols);
    B_P_k = zeros(rows*cols,k);
    % P_present_vec_init = P_present_vec;
    %%%%point reposition  
    while lamda < lamda_max
        disp('lamda:');
        disp(lamda);
        % update W_k, by D_P_2 < delta/lamda, not D_P
        %update by norm 1
%         D_P = D_P_cmpt(P_present,N_matrix,C_k);
%         W_k = D_P;
%         [r,c] = find( abs(D_P) < delta/lamda );
%         W_k( r , c ) = 0;   
        %update by norm 2
        D_P = D_P_cmpt(P_present,N_matrix,C_k);
        D_P_2 = D_P .* D_P;
        W_k = D_P;
        [r,c] = find( D_P_2 < delta/lamda );
        W_k( r , c ) = 0;   
        %substep2: update alpha
        for j = 1: 1 : k
            B_P_k( : , j ) = ( I_P -C_k( : , : , j ) ) * P_present_vec - N_vec_diag * W_k_j_vec_cmpt( W_k(:,j), cols );
        end
    
        b1 = zeros( rows * cols , 1 );
        for j = 1: 1 : k
            b1 = b1 + lamda * A_P_k( : , : , j )' * N_nn_diag * B_P_k( : , j );
        end

        b2 = N_vec_diag * (P_present_vec - P_init_vec);
        A_inverse = inv(S1 + lamda * S2);
        alpha_vec = A_inverse * ( -b1 - b2 ); 
        %substep3: update point position    
        P_present_vec = P_present_vec + N_vec_diag * alpha_vec;
        P_present =  vec2matrix( P_present_vec, cols );
        %%update lamda
        lamda = kapa_point*lamda;
    end
    %%%%point reposition    
    P_resosition = P_present;     
end

end

