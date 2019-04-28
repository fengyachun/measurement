clc;
clear all;
close all;
src = load('1_00_m_27_degree_points.txt');
src = src(:,1:1:3);
rows_depth = 640;
cols_depth = 480;
subsample_rate = 6;
P_src = [];
for i = 1:subsample_rate:rows_depth
    for j = 1:subsample_rate:cols_depth
        P_src = [P_src ; src(cols_depth*(i-1)+j , :)];
    end
end

P_src = P_src*1000;
[r,c] = size(P_src);
P_present_init = [];
for i = 1:1:r
    if norm(P_src(i,:))<50
        continue
    end
    if abs(P_src(i,1))>100
        continue;
    end
    if norm(P_src(i,:))>1500
        continue;
    end
    P_present_init = [P_present_init;P_src(i,:)];
    
end
% plot3(P_present_init(:,1),P_present_init(:,2),P_present_init(:,3),'.' );
% axis equal;
P_present = P_present_init(1:1:end, :);
% figure(2);
% plot3(P_present(:,1),P_present(:,2),P_present(:,3),'.' );
% axis equal;
save('1_00_m_27_degree_points.mat','P_present');
plot3(P_present(:,1),P_present(:,2),P_present(:,3),'.');
axis equal
