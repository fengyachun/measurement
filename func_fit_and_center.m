function [ n, average_center ] = func_fit_and_center( point )

[ rows , cols ] = size( point );
n = zeros( 1 , cols );
average_center = zeros( 1 , cols );
if rows == 1
    average_center = point;
    n = 0;
    return;
end
if rows < cols && rows > 1 
    average_center = sum( point ) / rows;
    n = 0;
    return;
end
average_center = sum( point ) / rows;
point = point - average_center;
[U,S,V] = svd(point);%%point_cloud = U*S*V'
n = ( V(:,end)/norm(V(:,end)) )';
end

