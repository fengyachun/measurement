function region_point = get_region_pointfrom_P_present(P_present,region_grow_group_index,each_region_point_num,index_sub_region )
%UNTITLED15 此处显示有关此函数的摘要
%   取出region index_sub_region 对应的区域的点云
if index_sub_region==1
        s = 0;      
        L = each_region_point_num(index_sub_region);
        tempIndex = region_grow_group_index( s + 1 : 1: s + L , 1 );
        region_point = P_present(tempIndex,:);
        return;
end
if index_sub_region > 1
        s = sum( each_region_point_num( 1 : 1 : index_sub_region - 1 ) );
        L = each_region_point_num(index_sub_region);
        tempIndex = region_grow_group_index( s + 1 : 1: s + L , 1 );
        region_point = P_present(tempIndex,:);
end 

end

