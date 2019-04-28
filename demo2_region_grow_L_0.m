clc;
clear all;
close all;

src = load('..\P_present_100_27.mat');
P_init = src.P_present;


% x_theta = (-28)*pi/180;%%%%%%%%%%%%rotate 
% rotateMatrix = [1,         0,             0;
%                 0,    cos(x_theta),  -sin(x_theta);
%                 0,    sin(x_theta),   cos(x_theta)];
% P_init = P_init*rotateMatrix;%
% P_init(:,2) = -P_init(:,2);


[ rows, cols ] = size(P_init);
%%%%%%%%%%%%%%%%%region growing
 k = 9;
theta_thresh = 10*pi/180;
thresh_index = 3;
multiple  = 10;
[ region_grow_group_index , each_region_point_num ] = func_region_grow_point_level( P_init , k , theta_thresh , thresh_index , multiple );
[region_num, c] = size(each_region_point_num);


figure(1);
figure(2);

n_arr = [];
center_arr = [];
%%%%%%%%%%%%%%%%%
for j = 1:1:region_num 
    %%%%%%%%get the point of sub-region
    region_point_temp = get_region_pointfrom_P_present(P_init , region_grow_group_index , each_region_point_num , j ); 
    %%%plot each sub-region
    figure(1);    
    plot3(region_point_temp(:,1) , region_point_temp(:,2) , region_point_temp(:,3),'.');  
    hold on;
    %%%
    %N_thresh
    [r_num, c_num] = size( region_point_temp );
    if  r_num<30
        continue
    end
    %%%%%%%%%%%%%%%%% L-0 smooth
    step_max = 1;
    k = 9;
    [region_P_resosition, region_N_smooth] = func_L_0_smooth_normal_and_reposition( region_point_temp, k, step_max); 
    %%%%%%%%%%%%%%%%%
    
    %%%%%%%%%%%%%%%%%
    [ n, average_center ] = func_fit_and_center( region_P_resosition );     
    
    %%%%%%%%%%%%%%%%%
    min_x = min(region_P_resosition(:,1));
    max_x = max(region_P_resosition(:,1));
    min_y = min(region_P_resosition(:,2));
    max_y = max(region_P_resosition(:,2));  
    min_z = min(region_P_resosition(:,3));
    max_z = max(region_P_resosition(:,3));  
    d_x = (max_x-min_x)/10;
    d_y = (max_y-min_y)/10;
    d_z = (max_z-min_z)/10;
    
    if d_y<d_z
        figure(2);        
        [X,Z]=meshgrid(min_x:d_x:max_x,min_z:d_z:max_z);
        Y = -1/n(2) * ( n(1) * X + n(3) * Z - n * average_center');
        mesh(X,Y,Z);   
        alpha(0);
        hold on;            
    end
    
    if d_y>d_z
        figure(2);        
        [X,Y]=meshgrid(min_x:d_x:max_x,min_y:d_y:max_y);
        Z = -1/n(3) * ( n(1) * X + n(2) * Y - n * average_center'); 
        mesh(X,Y,Z);   
        alpha(0);
        hold on;         
    end   
    
    %%%%%%%%%%%%%%%%%     
    
    n_arr = [n_arr;n];
    center_arr = [center_arr;average_center];
    
    %%%%%%%%%%%%%%%%%
    figure(2);    
    plot3(region_P_resosition(:,1) , region_P_resosition(:,2),region_P_resosition(:,3),'.k');  
    hold on;
    figure(5);
    hold on;
    plot(region_P_resosition(:,3),region_P_resosition(:,2),'.k');  
    hold off;
    %%%%%%%%%%%%%%%%%

end
%%sort
sort_y = center_arr(:,2);
[sort_res,sort_index] = sort(sort_y,'ascend');
n_arr = n_arr(sort_index,:);
center_arr = center_arr(sort_index,:);
%%%%%%%

figure(2);  
xlabel('x-axis');%
ylabel('y-axis');%
zlabel('z-axis');%
axis equal; 

figure(1);  
xlabel('x-axis');%
ylabel('y-axis');%
zlabel('z-axis');%
axis equal;

%%%%%%%%%%%%
[ r , c ] = size(n_arr);
across = [];
project_n_ = [];
for i = 1:1:r-1
    for j = i+1:1:r
       nn = cross(n_arr(i,:),n_arr(j,:)); 
       if (norm(nn) > 0.7)
           if nn(1)<0
               nn =-nn;
           end
           project_n_ = [project_n_;nn];
       end
    end
end
%
[r0,c0] = size(project_n_);
project_n = sum(project_n_)/r0;
project_n = project_n/norm(project_n);


%%%%%%%%%%%%%%%%%
% min_y = min(P_init(:,2));
% max_y = max(P_init(:,2));
% min_z = min(P_init(:,3));
% max_z = max(P_init(:,3));    
% d_y = (max_y-min_y)/20;
% d_z = (max_z-min_z)/20;
% average_center = sum( P_init ) / rows;
% figure(2);
% [Y,Z]=meshgrid(min_y:d_y:max_y,min_z:d_z:max_z);
% 
% % X = -1/project_n(1) * ( project_n(2) * Y + project_n(3) * Z - project_n * average_center');
% X = -1/project_n(1) * ( project_n(2) * Y + project_n(3) * Z );
% mesh(X,Y,Z);   
% alpha(0)
% hold on;

% figure(3);
% quiver3(0 , 0 , 0 ,project_n(1) , project_n(2) ,project_n(3) , 20 , 'r');
% hold on;
% xlabel('x-axis');%
% ylabel('y-axis');%
% zlabel('z-axis');%
% axis equal;

%%%%
n_status = ones(r,1);
n_index = zeros(r,1);
n_num = [];
for i = 1:1:r-1
    if n_status(i)==0
        continue;
    end
    itemindex = [];
    n_seed = n_arr(i,:);
    itemindex = [itemindex;i];
    for j = i+1:1:r
        if n_status(j)==0
            continue;
        end
        if abs(n_seed*n_arr(j,:)') < 0.7
           continue;
        end
        if norm(n_seed - n_arr(j,:))>1.8
            n_arr(j,:) = -n_arr(j,:);
        end
        itemindex = [itemindex;j];
        n_status(j) = 0;
    end
    item_num = length( itemindex ) ;
    n_index( sum( n_num ) + 1:1:sum( n_num ) + item_num ) = itemindex;
    n_num = [n_num ; item_num];
end

cross_p_arr = [];
for i = 1:1:r-1    
        A = [ project_n ; n_arr(i,:) ; n_arr(i+1,:)  ];
        B =[0 ; n_arr(i,:)*center_arr(i,:)' ; n_arr(i+1,:)*center_arr(i+1,:)'];
        P = inv(A)*B; 
        cross_p_arr = [cross_p_arr;P'];
end

figure(2);
hold on;
plot3(cross_p_arr(:,1),cross_p_arr(:,2),cross_p_arr(:,3),'*r','LineWidth',5);
hold off;

dis_step = [];
[rr,cc] = size(cross_p_arr);
for i = 1:1:rr-1
    distemp = norm(cross_p_arr(i,:)-cross_p_arr(i+1,:));
    dis_step = [dis_step ; distemp]; %the stairs height and width
end
    
figure(4)
plot3(P_init(:,1) , P_init(:,2) , P_init(:,3),'.b');  
xlabel('x-axis');
ylabel('y-axis');
zlabel('z-axis');
axis equal;

figure(5)
hold on;
plot(cross_p_arr(:,3),cross_p_arr(:,2),'*r','LineWidth',5);
xlabel('x-axis');
ylabel('y-axis');


