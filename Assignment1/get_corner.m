 [r,face] = read_ply('teabox.ply');
  %% Label the vertex of teabox
   v_1 = r(1,:);
   v_2 = r(2,:);
   v_3 = r(3,:);
   v_4 = r(4,:);
   v_5 = r(5,:);
   v_6 = r(6,:);
   v_7 = r(7,:);
   v_8 = r(8,:); 
   %% the order of the corner
   world_point_1 = [v_1 ; v_2; v_3; v_4; v_5; v_6];
   world_point_2 = [v_1 ; v_2; v_3; v_4; v_5; v_6; v_8];
   world_point_3 = [v_1 ; v_2; v_3; v_4; v_5; v_8];
   world_point_4 = [v_1 ; v_2; v_3; v_4; v_5; v_7; v_8];
   world_point_5 = [v_1 ; v_2; v_3; v_4; v_7; v_8];
   world_point_6 = [v_1 ; v_2; v_3; v_4; v_6; v_7; v_8];
   world_point_7 = [v_1 ; v_2; v_3; v_4; v_6; v_7];
   world_point_8 = [v_1 ; v_2; v_3; v_4; v_5; v_6; v_7];
 
   face = face + 1 ;%% Matlab array start from 1
   
   
%% create intrinst matrix 
intrinst_m = [2960.37845 0 0;
                0  2960.37845 0;
                1841.68855 1235.23369 1];
cameraParams = cameraParameters('IntrinsicMatrix',intrinst_m); 

   n = 1 ;
   images{n} = imread(sprintf('DSC_97%d.JPG',n+42));
   figure(1)
   imshow(images{n});
   [x,y]= getpts;%% select the corner in 2D image
   C(1:size(x),:,n) =[x,y];
   hold on;
   plot(C(:,1,n),C(:,2,n),'r*');

[worldOrientation,worldLocation] = estimateWorldCameraPose(C(1:size(x),:,n),world_point_1,cameraParams,'MaxReprojectionError' ,2);

figure(2);%% Show the estimateWorldCameraPose result
 col = [0; 6; 4; 3; 4; 6;0;4];
 patch('Faces',face,'Vertices',r,'FaceVertexCData',col,'FaceColor','interp');
 view(3);
 hold on
 plotCamera('Size',0.01,'Orientation',worldOrientation,'Location',worldLocation);
