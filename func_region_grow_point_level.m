function [ region_grow_group_index, each_region_point_num ] = func_region_grow_point_level( P_present , k , theta_thresh, thresh_index, multiple )

disp('region growing');

[rows,cols] = size(P_present); 
[ normal_P_present, knn_index, knn_P_present_dis ] = func_knn_normal_and_index_dis_cmpt( P_present, k );

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%region grow
process_status = ones( rows , 1 );
region_grow_group_index = zeros( rows , 1 );
each_region_point_num = [];
for i = 1:1:rows    
    if process_status(i)==0
        continue;
    end
    index_seed = [ ];
    index_seed = [ index_seed ; i ] ;    
    head_index = 1;
%     P_seed =        P_present( i , : );
%     N_seed = normal_P_present( i , : );
    while head_index - 1 < length( index_seed )
        seed_i =        index_seed( head_index ) ;
        if process_status( seed_i )==0
%             index_seed( head_index ) = 0;
            head_index = head_index + 1 ; 
            continue;
        end
        
        P_seed =        P_present( seed_i , : );
        N_seed = normal_P_present( seed_i , : );
        dis_thresh_point_plane = knn_P_present_dis( seed_i , thresh_index ) ;%
        dis_thresh_point_point = multiple * dis_thresh_point_plane;   %%
        for j = 2:1:k
            seed_knn_i_j = knn_index( seed_i , j );
            
            if process_status( seed_knn_i_j ) == 0
                continue;
            end
            
            P_seed_j =                       P_present( seed_knn_i_j , : );
            N_seed_j =                normal_P_present( seed_knn_i_j , : );
            
            if acos( abs( N_seed_j*N_seed' ) ) > theta_thresh
                continue;
            end
            
            if abs( ( P_seed_j - P_seed )*N_seed' ) > dis_thresh_point_plane
                continue;
            end
            
            if norm(  P_seed_j - P_seed  ) > dis_thresh_point_point
                continue;
            end
            [ r , c ] = find( ( index_seed - seed_knn_i_j ) ==0 );
            if length( r ) > 0
                continue;
            end
            index_seed = [index_seed; seed_knn_i_j ] ;            
        end
        process_status( seed_i ) = 0 ;
        head_index = head_index + 1 ;
    end
    disp('cluster');
    disp( size( index_seed ) );

    new_region_num = length( index_seed ) ;
    region_grow_group_index(sum( each_region_point_num ) + 1 : 1 : sum( each_region_point_num ) + new_region_num) = index_seed;    
    each_region_point_num = [each_region_point_num; new_region_num];
    
end


end

