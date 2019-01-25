


intrinst_m = [2960.37845 0 0;
                0  2960.37845 0;
                1841.68855 1235.23369 1];
cameraParams = cameraParameters('IntrinsicMatrix',intrinst_m);
            
image = imread('DSC_9775.JPG');            
[rotationMatrix,translationVector] = Ransac( image,500,1 );


%%
% rotationMatrix = [-0.863437229252907 -0.0109928051817204 -0.504336504106434;
%                    0.504453465762010 -0.0221623476498083 -0.863154407523559;
%                   -0.00168879269284626 -0.999693947456363 0.0246811547208555];
% translationVector = [0.0274285961791067 0.104845265182136 0.573097471583233];
% rotationMatrix = [-0.968953305627949,-0.0412070078860021,-0.243785713309358;0.240343640108527,-0.388304026566212,-0.889637520348512;-0.0580036737806554,-0.920609561895568,0.386152312403276];
% translationVector = [0.0225297990117465,0.104282195302060,0.503091570809018];
peak_thresh = 2;
% result_pose = [] ;
result_pose = {} ;
[orientation location] = extrinsicsToCameraPose(rotationMatrix,translationVector);
result_pose{1} = [orientation location'];


for j = 775 : 820

image_pre = imread(sprintf('DSC_9%d.JPG',j));
X = rgb2gray(image_pre);
[fa,da] = vl_sift(single(X),'PeakThresh', peak_thresh) ;

[ result D s] = image2world_ex2(image_pre,fa,da, rotationMatrix,translationVector, cameraParams);

image_curr = imread(sprintf('DSC_9%d.JPG',j+1));
X2 = rgb2gray(image_curr);
[fb,db] = vl_sift(single(X2),'PeakThresh', peak_thresh) ;

[matches, scores] = vl_ubcmatch(D, db,2) ;

image_point_cur = fb(1:2,matches(2,:))';
image_point_pre = s(:,matches(1,:))';

% figure(j);
% showMatchedFeatures(image_pre,image_curr,image_point_pre,image_point_cur,'montage','PlotOptions',{'ro','go','y--'})
% saveas(figure(j),sprintf('trackind feature_%d.png',j-773));
% figure(1);
% t = worldToImage(cameraParams,rotationMatrix,translationVector,result(matches(1,:),:));
% 
%         imshow(image_curr)
%         hold on ;
%         plot(t(:,1),t(:,2),'rx','Color','g');

rotationVector = rotationMatrixToVector(rotationMatrix);
theta = [rotationVector translationVector];
N = 100 ;
lamda = 0.001 ;
c = 4.685 ; 

M = result(matches(1,:),:) ; 
m = image_point_cur ;

for i = 1 : N
    
    e1 = projection_err(cameraParams,theta,M,m);
    
%     alpha =  1.48257968 * median(e1)  ;
%     wi = e1/alpha ;
%     remove = find((wi) > c );
%     v = (1-(wi/(c.^2)).^2).^2 ;
%     v(remove) = 0 ;
% 
%     W = diag(v);
 


    err = [e1(1:size(M,1)) e1(size(M,1)+1 :2*size(M,1))];
    err = abs(err(:,1)) + abs(err(:,2));
    alpha =  1.48257968 * median(err)  ;
    remove = find((err/alpha) > c );
    v = (1-(err/alpha/(c.^2)).^2).^2 ;
    v(remove) = 0 ;
    v = [v ; v];
    W = diag(v);
    e1(find(v == 0 )) = 0;
    er = e1' * e1 ;
    if i == 1 
        fprintf('initail %d err %d\n' ,j-773, norm(e1));
        
    end
    
    e = 0.000001;
    J = zeros(size(e1,1),6);
    J(:, 1) = ( projection_err(cameraParams,theta+[e 0 0 0 0 0],M,m) - e1) / e;
    J(:, 2) = ( projection_err(cameraParams,theta+[0 e 0 0 0 0],M,m)- e1) / e;
    J(:, 3) = ( projection_err(cameraParams,theta+[0 0 e 0 0 0],M,m)- e1) / e;
    J(:, 4) = ( projection_err(cameraParams,theta+[0 0 0 e 0 0],M,m)- e1) / e;
    J(:, 5) = ( projection_err(cameraParams,theta+[0 0 0 0 e 0],M,m)- e1) / e;
    J(:, 6) = ( projection_err(cameraParams,theta+[0 0 0 0 0 e],M,m)- e1) / e;
%     J1  = Jacob( theta,M,m );
%     delta = -inv(J'* J + lamda * eye(6)) * (J' * e1) ;
    delta = -inv(J'* W* J + lamda * eye(6)) * (J' *W * e1) ;
    e_new =  projection_err(cameraParams,theta+delta',M,m);
    e_new(find(v == 0 )) = 0 ;
    err_new = e_new'* e_new ;
    if norm(err_new) > norm(er) 
        lamda = 10 * lamda ;
%         disp('gradient');
    else 
        lamda = 0.1 * lamda ;
        theta = theta + delta' ;
%         fprintf('delta %d\n',delta);
      
        

%         disp('newton');        
    end
    
    if norm(delta) < 1e-10
        break ; 
    end
        
    
    
end
        
fprintf('After %d err %d\n',j-773,norm(e_new));

rotationMatrix = rotationVectorToMatrix(theta(1:3));
translationVector = theta(4:6) ;
[r,face] = read_ply('teabox.ply');
point = worldToImage(cameraParams,rotationMatrix,theta(4:6),M);
line_point = worldToImage(cameraParams,rotationMatrix,theta(4:6),r);
figure(j-773);
        imshow(image_curr)
       
        hold on ;
plot(point(:,1), point(:,2),'rx','Color','g') ;        
plot([line_point(1,1),line_point(2,1)],[line_point(1,2),line_point(2,2)],'LineWidth', 1, 'Color','b');
plot([line_point(3,1),line_point(4,1)],[line_point(3,2),line_point(4,2)],'LineWidth', 1, 'Color','b');
plot([line_point(5,1),line_point(6,1)],[line_point(5,2),line_point(6,2)],'LineWidth', 1, 'Color','b');
plot([line_point(7,1),line_point(8,1)],[line_point(7,2),line_point(8,2)],'LineWidth', 1, 'Color','b');
plot([line_point(2,1),line_point(3,1)],[line_point(2,2),line_point(3,2)],'LineWidth', 1, 'Color','b');
plot([line_point(1,1),line_point(4,1)],[line_point(1,2),line_point(4,2)],'LineWidth', 1, 'Color','b');
plot([line_point(6,1),line_point(7,1)],[line_point(6,2),line_point(7,2)],'LineWidth', 1, 'Color','b');
plot([line_point(5,1),line_point(8,1)],[line_point(5,2),line_point(8,2)],'LineWidth', 1, 'Color','b');
plot([line_point(4,1),line_point(8,1)],[line_point(4,2),line_point(8,2)],'LineWidth', 1, 'Color','b');
plot([line_point(3,1),line_point(7,1)],[line_point(3,2),line_point(7,2)],'LineWidth', 1, 'Color','b');
plot([line_point(2,1),line_point(6,1)],[line_point(2,2),line_point(6,2)],'LineWidth', 1, 'Color','b');
plot([line_point(1,1),line_point(5,1)],[line_point(1,2),line_point(5,2)],'LineWidth', 1, 'Color','b');

saveas(figure(j-773),sprintf('trackind_%d.png',j-773));
[orientation,location] = extrinsicsToCameraPose(rotationMatrix,translationVector);
result_pose{j-773} = [orientation location'];
end
%%
[r,face] = read_ply('teabox.ply');

 face = face + 1 ;
 grid on;
 col = [0; 6; 4; 3; 4; 6;0;4];
 figure(1);
 patch('Faces',face,'Vertices',r,'FaceVertexCData',col,'FaceColor','interp');
 view(3);
 hold on
 
 a = [1:size(result_pose,2)]'; b = num2str(a); c = cellstr(b);
 for i = 1 : size(result_pose,2)
 result = result_pose{i};
 orientation = result(1:3,1:3);
 location = result(1:3,4);
 
  
text( location(1),location(2),location(3), c(i),'Color','r','FontSize',20);
 plotCamera('Size',0.01,'Orientation',orientation,'Location',location,'Color','b');
 hold on ;
 end
 %%
% 
% syms u0 v0 f 
% syms tx ty tz wx wy wz 
% syms X Y Z x y 
% 
% 
% K=[f 0 u0;
%     0 f v0;
%     0 0 1];
%  
% % Expression for the rotation matrix based on the Rodrigues formula
%  
% theta=sqrt(wx^2+wy^2+wz^2);
% omega=  [0 -wz wy;
%          wz 0 -wx;
%         -wy wx 0;];
% R = eye(3) + (sin(theta)/theta)*omega + ((1-cos(theta))/theta^2)*(omega*omega);
% % Expression for the translation vector
% 
% t=[tx;ty;tz];
% % perspective projection of the model point (X,Y,Z)
%  
% uvs=K*[R t]*[X; Y; Z; 1];
% u=uvs(1)/uvs(3);
% v=uvs(2)/uvs(3);
% % calculate the geometric distance in x and y direction
%  
% % u,v = the x and y positions of the projection of the corresponding model point
%  
% dx=u;
% dy=v;
% e = sqrt((x - dx)^2 + (y - dy)^2) ;
% % Evaluate the symbolic expression of the Jacobian w.r.t. the estimated parameters
% Jx=jacobian(dx,[wx wy wz tx ty tz]);
% Jy=jacobian(dy,[wx wy wz tx ty tz]);
% diary ('myVariable.txt');
% Jx 
% Jy
% % your symbolic variable
% diary off
%  
% %%
%     tx = translationVector(1);
%     ty = translationVector(2);
%     tz = translationVector(3);
%     rotationVector = rotationMatrixToVector(rotationMatrix);
%     wx =  rotationVector(1);
%     wy =  rotationVector(2);
%     wz =  rotationVector(3);
%     f = 2960.37845;
%     u0 = 1841.68855;
%     v0 = 1235.23369;
%     X =1 ;
%    
%     Y =2;
%     Z =3;
%     Jx = [ - (X*(f*((2*wx*(cos((wx^2 + wy^2 + wz^2)^(1/2)) - 1)*(wy^2 + wz^2))/(wx^2 + wy^2 + wz^2)^2 + (wx*sin((wx^2 + wy^2 + wz^2)^(1/2))*(wy^2 + wz^2))/(wx^2 + wy^2 + wz^2)^(3/2)) - u0*((wx*wy*sin((wx^2 + wy^2 + wz^2)^(1/2)))/(wx^2 + wy^2 + wz^2)^(3/2) - (wx*wy*cos((wx^2 + wy^2 + wz^2)^(1/2)))/(wx^2 + wy^2 + wz^2) - (wz*(cos((wx^2 + wy^2 + wz^2)^(1/2)) - 1))/(wx^2 + wy^2 + wz^2) + (wx^2*wz*sin((wx^2 + wy^2 + wz^2)^(1/2)))/(wx^2 + wy^2 + wz^2)^(3/2) + (2*wx^2*wz*(cos((wx^2 + wy^2 + wz^2)^(1/2)) - 1))/(wx^2 + wy^2 + wz^2)^2)) - Y*(f*((wx*wz*sin((wx^2 + wy^2 + wz^2)^(1/2)))/(wx^2 + wy^2 + wz^2)^(3/2) - (wx*wz*cos((wx^2 + wy^2 + wz^2)^(1/2)))/(wx^2 + wy^2 + wz^2) - (wy*(cos((wx^2 + wy^2 + wz^2)^(1/2)) - 1))/(wx^2 + wy^2 + wz^2) + (wx^2*wy*sin((wx^2 + wy^2 + wz^2)^(1/2)))/(wx^2 + wy^2 + wz^2)^(3/2) + (2*wx^2*wy*(cos((wx^2 + wy^2 + wz^2)^(1/2)) - 1))/(wx^2 + wy^2 + wz^2)^2) + u0*(sin((wx^2 + wy^2 + wz^2)^(1/2))/(wx^2 + wy^2 + wz^2)^(1/2) + (wx^2*cos((wx^2 + wy^2 + wz^2)^(1/2)))/(wx^2 + wy^2 + wz^2) - (wx^2*sin((wx^2 + wy^2 + wz^2)^(1/2)))/(wx^2 + wy^2 + wz^2)^(3/2) + (2*wx*wy*wz*(cos((wx^2 + wy^2 + wz^2)^(1/2)) - 1))/(wx^2 + wy^2 + wz^2)^2 + (wx*wy*wz*sin((wx^2 + wy^2 + wz^2)^(1/2)))/(wx^2 + wy^2 + wz^2)^(3/2))) + Z*(u0*((2*wx*(cos((wx^2 + wy^2 + wz^2)^(1/2)) - 1)*(wx^2 + wy^2))/(wx^2 + wy^2 + wz^2)^2 - (2*wx*(cos((wx^2 + wy^2 + wz^2)^(1/2)) - 1))/(wx^2 + wy^2 + wz^2) + (wx*sin((wx^2 + wy^2 + wz^2)^(1/2))*(wx^2 + wy^2))/(wx^2 + wy^2 + wz^2)^(3/2)) - f*((wx*wy*cos((wx^2 + wy^2 + wz^2)^(1/2)))/(wx^2 + wy^2 + wz^2) - (wz*(cos((wx^2 + wy^2 + wz^2)^(1/2)) - 1))/(wx^2 + wy^2 + wz^2) - (wx*wy*sin((wx^2 + wy^2 + wz^2)^(1/2)))/(wx^2 + wy^2 + wz^2)^(3/2) + (wx^2*wz*sin((wx^2 + wy^2 + wz^2)^(1/2)))/(wx^2 + wy^2 + wz^2)^(3/2) + (2*wx^2*wz*(cos((wx^2 + wy^2 + wz^2)^(1/2)) - 1))/(wx^2 + wy^2 + wz^2)^2)))/(tz - X*((wy*sin((wx^2 + wy^2 + wz^2)^(1/2)))/(wx^2 + wy^2 + wz^2)^(1/2) + (wx*wz*(cos((wx^2 + wy^2 + wz^2)^(1/2)) - 1))/(wx^2 + wy^2 + wz^2)) + Y*((wx*sin((wx^2 + wy^2 + wz^2)^(1/2)))/(wx^2 + wy^2 + wz^2)^(1/2) - (wy*wz*(cos((wx^2 + wy^2 + wz^2)^(1/2)) - 1))/(wx^2 + wy^2 + wz^2)) + Z*(((cos((wx^2 + wy^2 + wz^2)^(1/2)) - 1)*(wx^2 + wy^2))/(wx^2 + wy^2 + wz^2) + 1)) - ((X*((wx*wy*sin((wx^2 + wy^2 + wz^2)^(1/2)))/(wx^2 + wy^2 + wz^2)^(3/2) - (wx*wy*cos((wx^2 + wy^2 + wz^2)^(1/2)))/(wx^2 + wy^2 + wz^2) - (wz*(cos((wx^2 + wy^2 + wz^2)^(1/2)) - 1))/(wx^2 + wy^2 + wz^2) + (wx^2*wz*sin((wx^2 + wy^2 + wz^2)^(1/2)))/(wx^2 + wy^2 + wz^2)^(3/2) + (2*wx^2*wz*(cos((wx^2 + wy^2 + wz^2)^(1/2)) - 1))/(wx^2 + wy^2 + wz^2)^2) - Z*((2*wx*(cos((wx^2 + wy^2 + wz^2)^(1/2)) - 1)*(wx^2 + wy^2))/(wx^2 + wy^2 + wz^2)^2 - (2*wx*(cos((wx^2 + wy^2 + wz^2)^(1/2)) - 1))/(wx^2 + wy^2 + wz^2) + (wx*sin((wx^2 + wy^2 + wz^2)^(1/2))*(wx^2 + wy^2))/(wx^2 + wy^2 + wz^2)^(3/2)) + Y*(sin((wx^2 + wy^2 + wz^2)^(1/2))/(wx^2 + wy^2 + wz^2)^(1/2) + (wx^2*cos((wx^2 + wy^2 + wz^2)^(1/2)))/(wx^2 + wy^2 + wz^2) - (wx^2*sin((wx^2 + wy^2 + wz^2)^(1/2)))/(wx^2 + wy^2 + wz^2)^(3/2) + (2*wx*wy*wz*(cos((wx^2 + wy^2 + wz^2)^(1/2)) - 1))/(wx^2 + wy^2 + wz^2)^2 + (wx*wy*wz*sin((wx^2 + wy^2 + wz^2)^(1/2)))/(wx^2 + wy^2 + wz^2)^(3/2)))*(f*tx - Y*(f*((wz*sin((wx^2 + wy^2 + wz^2)^(1/2)))/(wx^2 + wy^2 + wz^2)^(1/2) + (wx*wy*(cos((wx^2 + wy^2 + wz^2)^(1/2)) - 1))/(wx^2 + wy^2 + wz^2)) - u0*((wx*sin((wx^2 + wy^2 + wz^2)^(1/2)))/(wx^2 + wy^2 + wz^2)^(1/2) - (wy*wz*(cos((wx^2 + wy^2 + wz^2)^(1/2)) - 1))/(wx^2 + wy^2 + wz^2))) + tz*u0 - X*(u0*((wy*sin((wx^2 + wy^2 + wz^2)^(1/2)))/(wx^2 + wy^2 + wz^2)^(1/2) + (wx*wz*(cos((wx^2 + wy^2 + wz^2)^(1/2)) - 1))/(wx^2 + wy^2 + wz^2)) - f*(((cos((wx^2 + wy^2 + wz^2)^(1/2)) - 1)*(wy^2 + wz^2))/(wx^2 + wy^2 + wz^2) + 1)) + Z*(f*((wy*sin((wx^2 + wy^2 + wz^2)^(1/2)))/(wx^2 + wy^2 + wz^2)^(1/2) - (wx*wz*(cos((wx^2 + wy^2 + wz^2)^(1/2)) - 1))/(wx^2 + wy^2 + wz^2)) + u0*(((cos((wx^2 + wy^2 + wz^2)^(1/2)) - 1)*(wx^2 + wy^2))/(wx^2 + wy^2 + wz^2) + 1))))/(tz - X*((wy*sin((wx^2 + wy^2 + wz^2)^(1/2)))/(wx^2 + wy^2 + wz^2)^(1/2) + (wx*wz*(cos((wx^2 + wy^2 + wz^2)^(1/2)) - 1))/(wx^2 + wy^2 + wz^2)) + Y*((wx*sin((wx^2 + wy^2 + wz^2)^(1/2)))/(wx^2 + wy^2 + wz^2)^(1/2) - (wy*wz*(cos((wx^2 + wy^2 + wz^2)^(1/2)) - 1))/(wx^2 + wy^2 + wz^2)) + Z*(((cos((wx^2 + wy^2 + wz^2)^(1/2)) - 1)*(wx^2 + wy^2))/(wx^2 + wy^2 + wz^2) + 1))^2, - (X*(f*((2*wy*(cos((wx^2 + wy^2 + wz^2)^(1/2)) - 1)*(wy^2 + wz^2))/(wx^2 + wy^2 + wz^2)^2 - (2*wy*(cos((wx^2 + wy^2 + wz^2)^(1/2)) - 1))/(wx^2 + wy^2 + wz^2) + (wy*sin((wx^2 + wy^2 + wz^2)^(1/2))*(wy^2 + wz^2))/(wx^2 + wy^2 + wz^2)^(3/2)) - u0*((wy^2*sin((wx^2 + wy^2 + wz^2)^(1/2)))/(wx^2 + wy^2 + wz^2)^(3/2) - (wy^2*cos((wx^2 + wy^2 + wz^2)^(1/2)))/(wx^2 + wy^2 + wz^2) - sin((wx^2 + wy^2 + wz^2)^(1/2))/(wx^2 + wy^2 + wz^2)^(1/2) + (2*wx*wy*wz*(cos((wx^2 + wy^2 + wz^2)^(1/2)) - 1))/(wx^2 + wy^2 + wz^2)^2 + (wx*wy*wz*sin((wx^2 + wy^2 + wz^2)^(1/2)))/(wx^2 + wy^2 + wz^2)^(3/2))) + Z*(u0*((2*wy*(cos((wx^2 + wy^2 + wz^2)^(1/2)) - 1)*(wx^2 + wy^2))/(wx^2 + wy^2 + wz^2)^2 - (2*wy*(cos((wx^2 + wy^2 + wz^2)^(1/2)) - 1))/(wx^2 + wy^2 + wz^2) + (wy*sin((wx^2 + wy^2 + wz^2)^(1/2))*(wx^2 + wy^2))/(wx^2 + wy^2 + wz^2)^(3/2)) - f*(sin((wx^2 + wy^2 + wz^2)^(1/2))/(wx^2 + wy^2 + wz^2)^(1/2) + (wy^2*cos((wx^2 + wy^2 + wz^2)^(1/2)))/(wx^2 + wy^2 + wz^2) - (wy^2*sin((wx^2 + wy^2 + wz^2)^(1/2)))/(wx^2 + wy^2 + wz^2)^(3/2) + (2*wx*wy*wz*(cos((wx^2 + wy^2 + wz^2)^(1/2)) - 1))/(wx^2 + wy^2 + wz^2)^2 + (wx*wy*wz*sin((wx^2 + wy^2 + wz^2)^(1/2)))/(wx^2 + wy^2 + wz^2)^(3/2))) - Y*(f*((wy*wz*sin((wx^2 + wy^2 + wz^2)^(1/2)))/(wx^2 + wy^2 + wz^2)^(3/2) - (wy*wz*cos((wx^2 + wy^2 + wz^2)^(1/2)))/(wx^2 + wy^2 + wz^2) - (wx*(cos((wx^2 + wy^2 + wz^2)^(1/2)) - 1))/(wx^2 + wy^2 + wz^2) + (wx*wy^2*sin((wx^2 + wy^2 + wz^2)^(1/2)))/(wx^2 + wy^2 + wz^2)^(3/2) + (2*wx*wy^2*(cos((wx^2 + wy^2 + wz^2)^(1/2)) - 1))/(wx^2 + wy^2 + wz^2)^2) + u0*((wx*wy*cos((wx^2 + wy^2 + wz^2)^(1/2)))/(wx^2 + wy^2 + wz^2) - (wz*(cos((wx^2 + wy^2 + wz^2)^(1/2)) - 1))/(wx^2 + wy^2 + wz^2) - (wx*wy*sin((wx^2 + wy^2 + wz^2)^(1/2)))/(wx^2 + wy^2 + wz^2)^(3/2) + (wy^2*wz*sin((wx^2 + wy^2 + wz^2)^(1/2)))/(wx^2 + wy^2 + wz^2)^(3/2) + (2*wy^2*wz*(cos((wx^2 + wy^2 + wz^2)^(1/2)) - 1))/(wx^2 + wy^2 + wz^2)^2)))/(tz - X*((wy*sin((wx^2 + wy^2 + wz^2)^(1/2)))/(wx^2 + wy^2 + wz^2)^(1/2) + (wx*wz*(cos((wx^2 + wy^2 + wz^2)^(1/2)) - 1))/(wx^2 + wy^2 + wz^2)) + Y*((wx*sin((wx^2 + wy^2 + wz^2)^(1/2)))/(wx^2 + wy^2 + wz^2)^(1/2) - (wy*wz*(cos((wx^2 + wy^2 + wz^2)^(1/2)) - 1))/(wx^2 + wy^2 + wz^2)) + Z*(((cos((wx^2 + wy^2 + wz^2)^(1/2)) - 1)*(wx^2 + wy^2))/(wx^2 + wy^2 + wz^2) + 1)) - ((Y*((wx*wy*cos((wx^2 + wy^2 + wz^2)^(1/2)))/(wx^2 + wy^2 + wz^2) - (wz*(cos((wx^2 + wy^2 + wz^2)^(1/2)) - 1))/(wx^2 + wy^2 + wz^2) - (wx*wy*sin((wx^2 + wy^2 + wz^2)^(1/2)))/(wx^2 + wy^2 + wz^2)^(3/2) + (wy^2*wz*sin((wx^2 + wy^2 + wz^2)^(1/2)))/(wx^2 + wy^2 + wz^2)^(3/2) + (2*wy^2*wz*(cos((wx^2 + wy^2 + wz^2)^(1/2)) - 1))/(wx^2 + wy^2 + wz^2)^2) - Z*((2*wy*(cos((wx^2 + wy^2 + wz^2)^(1/2)) - 1)*(wx^2 + wy^2))/(wx^2 + wy^2 + wz^2)^2 - (2*wy*(cos((wx^2 + wy^2 + wz^2)^(1/2)) - 1))/(wx^2 + wy^2 + wz^2) + (wy*sin((wx^2 + wy^2 + wz^2)^(1/2))*(wx^2 + wy^2))/(wx^2 + wy^2 + wz^2)^(3/2)) + X*((wy^2*sin((wx^2 + wy^2 + wz^2)^(1/2)))/(wx^2 + wy^2 + wz^2)^(3/2) - (wy^2*cos((wx^2 + wy^2 + wz^2)^(1/2)))/(wx^2 + wy^2 + wz^2) - sin((wx^2 + wy^2 + wz^2)^(1/2))/(wx^2 + wy^2 + wz^2)^(1/2) + (2*wx*wy*wz*(cos((wx^2 + wy^2 + wz^2)^(1/2)) - 1))/(wx^2 + wy^2 + wz^2)^2 + (wx*wy*wz*sin((wx^2 + wy^2 + wz^2)^(1/2)))/(wx^2 + wy^2 + wz^2)^(3/2)))*(f*tx - Y*(f*((wz*sin((wx^2 + wy^2 + wz^2)^(1/2)))/(wx^2 + wy^2 + wz^2)^(1/2) + (wx*wy*(cos((wx^2 + wy^2 + wz^2)^(1/2)) - 1))/(wx^2 + wy^2 + wz^2)) - u0*((wx*sin((wx^2 + wy^2 + wz^2)^(1/2)))/(wx^2 + wy^2 + wz^2)^(1/2) - (wy*wz*(cos((wx^2 + wy^2 + wz^2)^(1/2)) - 1))/(wx^2 + wy^2 + wz^2))) + tz*u0 - X*(u0*((wy*sin((wx^2 + wy^2 + wz^2)^(1/2)))/(wx^2 + wy^2 + wz^2)^(1/2) + (wx*wz*(cos((wx^2 + wy^2 + wz^2)^(1/2)) - 1))/(wx^2 + wy^2 + wz^2)) - f*(((cos((wx^2 + wy^2 + wz^2)^(1/2)) - 1)*(wy^2 + wz^2))/(wx^2 + wy^2 + wz^2) + 1)) + Z*(f*((wy*sin((wx^2 + wy^2 + wz^2)^(1/2)))/(wx^2 + wy^2 + wz^2)^(1/2) - (wx*wz*(cos((wx^2 + wy^2 + wz^2)^(1/2)) - 1))/(wx^2 + wy^2 + wz^2)) + u0*(((cos((wx^2 + wy^2 + wz^2)^(1/2)) - 1)*(wx^2 + wy^2))/(wx^2 + wy^2 + wz^2) + 1))))/(tz - X*((wy*sin((wx^2 + wy^2 + wz^2)^(1/2)))/(wx^2 + wy^2 + wz^2)^(1/2) + (wx*wz*(cos((wx^2 + wy^2 + wz^2)^(1/2)) - 1))/(wx^2 + wy^2 + wz^2)) + Y*((wx*sin((wx^2 + wy^2 + wz^2)^(1/2)))/(wx^2 + wy^2 + wz^2)^(1/2) - (wy*wz*(cos((wx^2 + wy^2 + wz^2)^(1/2)) - 1))/(wx^2 + wy^2 + wz^2)) + Z*(((cos((wx^2 + wy^2 + wz^2)^(1/2)) - 1)*(wx^2 + wy^2))/(wx^2 + wy^2 + wz^2) + 1))^2, - (Z*(u0*((2*wz*(cos((wx^2 + wy^2 + wz^2)^(1/2)) - 1)*(wx^2 + wy^2))/(wx^2 + wy^2 + wz^2)^2 + (wz*sin((wx^2 + wy^2 + wz^2)^(1/2))*(wx^2 + wy^2))/(wx^2 + wy^2 + wz^2)^(3/2)) - f*((wy*wz*cos((wx^2 + wy^2 + wz^2)^(1/2)))/(wx^2 + wy^2 + wz^2) - (wx*(cos((wx^2 + wy^2 + wz^2)^(1/2)) - 1))/(wx^2 + wy^2 + wz^2) - (wy*wz*sin((wx^2 + wy^2 + wz^2)^(1/2)))/(wx^2 + wy^2 + wz^2)^(3/2) + (wx*wz^2*sin((wx^2 + wy^2 + wz^2)^(1/2)))/(wx^2 + wy^2 + wz^2)^(3/2) + (2*wx*wz^2*(cos((wx^2 + wy^2 + wz^2)^(1/2)) - 1))/(wx^2 + wy^2 + wz^2)^2)) - Y*(u0*((wx*wz*cos((wx^2 + wy^2 + wz^2)^(1/2)))/(wx^2 + wy^2 + wz^2) - (wy*(cos((wx^2 + wy^2 + wz^2)^(1/2)) - 1))/(wx^2 + wy^2 + wz^2) - (wx*wz*sin((wx^2 + wy^2 + wz^2)^(1/2)))/(wx^2 + wy^2 + wz^2)^(3/2) + (wy*wz^2*sin((wx^2 + wy^2 + wz^2)^(1/2)))/(wx^2 + wy^2 + wz^2)^(3/2) + (2*wy*wz^2*(cos((wx^2 + wy^2 + wz^2)^(1/2)) - 1))/(wx^2 + wy^2 + wz^2)^2) + f*((wz^2*sin((wx^2 + wy^2 + wz^2)^(1/2)))/(wx^2 + wy^2 + wz^2)^(3/2) - (wz^2*cos((wx^2 + wy^2 + wz^2)^(1/2)))/(wx^2 + wy^2 + wz^2) - sin((wx^2 + wy^2 + wz^2)^(1/2))/(wx^2 + wy^2 + wz^2)^(1/2) + (2*wx*wy*wz*(cos((wx^2 + wy^2 + wz^2)^(1/2)) - 1))/(wx^2 + wy^2 + wz^2)^2 + (wx*wy*wz*sin((wx^2 + wy^2 + wz^2)^(1/2)))/(wx^2 + wy^2 + wz^2)^(3/2))) + X*(f*((2*wz*(cos((wx^2 + wy^2 + wz^2)^(1/2)) - 1)*(wy^2 + wz^2))/(wx^2 + wy^2 + wz^2)^2 - (2*wz*(cos((wx^2 + wy^2 + wz^2)^(1/2)) - 1))/(wx^2 + wy^2 + wz^2) + (wz*sin((wx^2 + wy^2 + wz^2)^(1/2))*(wy^2 + wz^2))/(wx^2 + wy^2 + wz^2)^(3/2)) - u0*((wy*wz*sin((wx^2 + wy^2 + wz^2)^(1/2)))/(wx^2 + wy^2 + wz^2)^(3/2) - (wy*wz*cos((wx^2 + wy^2 + wz^2)^(1/2)))/(wx^2 + wy^2 + wz^2) - (wx*(cos((wx^2 + wy^2 + wz^2)^(1/2)) - 1))/(wx^2 + wy^2 + wz^2) + (wx*wz^2*sin((wx^2 + wy^2 + wz^2)^(1/2)))/(wx^2 + wy^2 + wz^2)^(3/2) + (2*wx*wz^2*(cos((wx^2 + wy^2 + wz^2)^(1/2)) - 1))/(wx^2 + wy^2 + wz^2)^2)))/(tz - X*((wy*sin((wx^2 + wy^2 + wz^2)^(1/2)))/(wx^2 + wy^2 + wz^2)^(1/2) + (wx*wz*(cos((wx^2 + wy^2 + wz^2)^(1/2)) - 1))/(wx^2 + wy^2 + wz^2)) + Y*((wx*sin((wx^2 + wy^2 + wz^2)^(1/2)))/(wx^2 + wy^2 + wz^2)^(1/2) - (wy*wz*(cos((wx^2 + wy^2 + wz^2)^(1/2)) - 1))/(wx^2 + wy^2 + wz^2)) + Z*(((cos((wx^2 + wy^2 + wz^2)^(1/2)) - 1)*(wx^2 + wy^2))/(wx^2 + wy^2 + wz^2) + 1)) - ((X*((wy*wz*sin((wx^2 + wy^2 + wz^2)^(1/2)))/(wx^2 + wy^2 + wz^2)^(3/2) - (wy*wz*cos((wx^2 + wy^2 + wz^2)^(1/2)))/(wx^2 + wy^2 + wz^2) - (wx*(cos((wx^2 + wy^2 + wz^2)^(1/2)) - 1))/(wx^2 + wy^2 + wz^2) + (wx*wz^2*sin((wx^2 + wy^2 + wz^2)^(1/2)))/(wx^2 + wy^2 + wz^2)^(3/2) + (2*wx*wz^2*(cos((wx^2 + wy^2 + wz^2)^(1/2)) - 1))/(wx^2 + wy^2 + wz^2)^2) - Z*((2*wz*(cos((wx^2 + wy^2 + wz^2)^(1/2)) - 1)*(wx^2 + wy^2))/(wx^2 + wy^2 + wz^2)^2 + (wz*sin((wx^2 + wy^2 + wz^2)^(1/2))*(wx^2 + wy^2))/(wx^2 + wy^2 + wz^2)^(3/2)) + Y*((wx*wz*cos((wx^2 + wy^2 + wz^2)^(1/2)))/(wx^2 + wy^2 + wz^2) - (wy*(cos((wx^2 + wy^2 + wz^2)^(1/2)) - 1))/(wx^2 + wy^2 + wz^2) - (wx*wz*sin((wx^2 + wy^2 + wz^2)^(1/2)))/(wx^2 + wy^2 + wz^2)^(3/2) + (wy*wz^2*sin((wx^2 + wy^2 + wz^2)^(1/2)))/(wx^2 + wy^2 + wz^2)^(3/2) + (2*wy*wz^2*(cos((wx^2 + wy^2 + wz^2)^(1/2)) - 1))/(wx^2 + wy^2 + wz^2)^2))*(f*tx - Y*(f*((wz*sin((wx^2 + wy^2 + wz^2)^(1/2)))/(wx^2 + wy^2 + wz^2)^(1/2) + (wx*wy*(cos((wx^2 + wy^2 + wz^2)^(1/2)) - 1))/(wx^2 + wy^2 + wz^2)) - u0*((wx*sin((wx^2 + wy^2 + wz^2)^(1/2)))/(wx^2 + wy^2 + wz^2)^(1/2) - (wy*wz*(cos((wx^2 + wy^2 + wz^2)^(1/2)) - 1))/(wx^2 + wy^2 + wz^2))) + tz*u0 - X*(u0*((wy*sin((wx^2 + wy^2 + wz^2)^(1/2)))/(wx^2 + wy^2 + wz^2)^(1/2) + (wx*wz*(cos((wx^2 + wy^2 + wz^2)^(1/2)) - 1))/(wx^2 + wy^2 + wz^2)) - f*(((cos((wx^2 + wy^2 + wz^2)^(1/2)) - 1)*(wy^2 + wz^2))/(wx^2 + wy^2 + wz^2) + 1)) + Z*(f*((wy*sin((wx^2 + wy^2 + wz^2)^(1/2)))/(wx^2 + wy^2 + wz^2)^(1/2) - (wx*wz*(cos((wx^2 + wy^2 + wz^2)^(1/2)) - 1))/(wx^2 + wy^2 + wz^2)) + u0*(((cos((wx^2 + wy^2 + wz^2)^(1/2)) - 1)*(wx^2 + wy^2))/(wx^2 + wy^2 + wz^2) + 1))))/(tz - X*((wy*sin((wx^2 + wy^2 + wz^2)^(1/2)))/(wx^2 + wy^2 + wz^2)^(1/2) + (wx*wz*(cos((wx^2 + wy^2 + wz^2)^(1/2)) - 1))/(wx^2 + wy^2 + wz^2)) + Y*((wx*sin((wx^2 + wy^2 + wz^2)^(1/2)))/(wx^2 + wy^2 + wz^2)^(1/2) - (wy*wz*(cos((wx^2 + wy^2 + wz^2)^(1/2)) - 1))/(wx^2 + wy^2 + wz^2)) + Z*(((cos((wx^2 + wy^2 + wz^2)^(1/2)) - 1)*(wx^2 + wy^2))/(wx^2 + wy^2 + wz^2) + 1))^2, f/(tz - X*((wy*sin((wx^2 + wy^2 + wz^2)^(1/2)))/(wx^2 + wy^2 + wz^2)^(1/2) + (wx*wz*(cos((wx^2 + wy^2 + wz^2)^(1/2)) - 1))/(wx^2 + wy^2 + wz^2)) + Y*((wx*sin((wx^2 + wy^2 + wz^2)^(1/2)))/(wx^2 + wy^2 + wz^2)^(1/2) - (wy*wz*(cos((wx^2 + wy^2 + wz^2)^(1/2)) - 1))/(wx^2 + wy^2 + wz^2)) + Z*(((cos((wx^2 + wy^2 + wz^2)^(1/2)) - 1)*(wx^2 + wy^2))/(wx^2 + wy^2 + wz^2) + 1)), 0, u0/(tz - X*((wy*sin((wx^2 + wy^2 + wz^2)^(1/2)))/(wx^2 + wy^2 + wz^2)^(1/2) + (wx*wz*(cos((wx^2 + wy^2 + wz^2)^(1/2)) - 1))/(wx^2 + wy^2 + wz^2)) + Y*((wx*sin((wx^2 + wy^2 + wz^2)^(1/2)))/(wx^2 + wy^2 + wz^2)^(1/2) - (wy*wz*(cos((wx^2 + wy^2 + wz^2)^(1/2)) - 1))/(wx^2 + wy^2 + wz^2)) + Z*(((cos((wx^2 + wy^2 + wz^2)^(1/2)) - 1)*(wx^2 + wy^2))/(wx^2 + wy^2 + wz^2) + 1)) - (f*tx - Y*(f*((wz*sin((wx^2 + wy^2 + wz^2)^(1/2)))/(wx^2 + wy^2 + wz^2)^(1/2) + (wx*wy*(cos((wx^2 + wy^2 + wz^2)^(1/2)) - 1))/(wx^2 + wy^2 + wz^2)) - u0*((wx*sin((wx^2 + wy^2 + wz^2)^(1/2)))/(wx^2 + wy^2 + wz^2)^(1/2) - (wy*wz*(cos((wx^2 + wy^2 + wz^2)^(1/2)) - 1))/(wx^2 + wy^2 + wz^2))) + tz*u0 - X*(u0*((wy*sin((wx^2 + wy^2 + wz^2)^(1/2)))/(wx^2 + wy^2 + wz^2)^(1/2) + (wx*wz*(cos((wx^2 + wy^2 + wz^2)^(1/2)) - 1))/(wx^2 + wy^2 + wz^2)) - f*(((cos((wx^2 + wy^2 + wz^2)^(1/2)) - 1)*(wy^2 + wz^2))/(wx^2 + wy^2 + wz^2) + 1)) + Z*(f*((wy*sin((wx^2 + wy^2 + wz^2)^(1/2)))/(wx^2 + wy^2 + wz^2)^(1/2) - (wx*wz*(cos((wx^2 + wy^2 + wz^2)^(1/2)) - 1))/(wx^2 + wy^2 + wz^2)) + u0*(((cos((wx^2 + wy^2 + wz^2)^(1/2)) - 1)*(wx^2 + wy^2))/(wx^2 + wy^2 + wz^2) + 1)))/(tz - X*((wy*sin((wx^2 + wy^2 + wz^2)^(1/2)))/(wx^2 + wy^2 + wz^2)^(1/2) + (wx*wz*(cos((wx^2 + wy^2 + wz^2)^(1/2)) - 1))/(wx^2 + wy^2 + wz^2)) + Y*((wx*sin((wx^2 + wy^2 + wz^2)^(1/2)))/(wx^2 + wy^2 + wz^2)^(1/2) - (wy*wz*(cos((wx^2 + wy^2 + wz^2)^(1/2)) - 1))/(wx^2 + wy^2 + wz^2)) + Z*(((cos((wx^2 + wy^2 + wz^2)^(1/2)) - 1)*(wx^2 + wy^2))/(wx^2 + wy^2 + wz^2) + 1))^2];
%     Jy = [ - (Y*(f*((2*wx*(cos((wx^2 + wy^2 + wz^2)^(1/2)) - 1)*(wx^2 + wz^2))/(wx^2 + wy^2 + wz^2)^2 - (2*wx*(cos((wx^2 + wy^2 + wz^2)^(1/2)) - 1))/(wx^2 + wy^2 + wz^2) + (wx*sin((wx^2 + wy^2 + wz^2)^(1/2))*(wx^2 + wz^2))/(wx^2 + wy^2 + wz^2)^(3/2)) - v0*(sin((wx^2 + wy^2 + wz^2)^(1/2))/(wx^2 + wy^2 + wz^2)^(1/2) + (wx^2*cos((wx^2 + wy^2 + wz^2)^(1/2)))/(wx^2 + wy^2 + wz^2) - (wx^2*sin((wx^2 + wy^2 + wz^2)^(1/2)))/(wx^2 + wy^2 + wz^2)^(3/2) + (2*wx*wy*wz*(cos((wx^2 + wy^2 + wz^2)^(1/2)) - 1))/(wx^2 + wy^2 + wz^2)^2 + (wx*wy*wz*sin((wx^2 + wy^2 + wz^2)^(1/2)))/(wx^2 + wy^2 + wz^2)^(3/2))) + Z*(v0*((2*wx*(cos((wx^2 + wy^2 + wz^2)^(1/2)) - 1)*(wx^2 + wy^2))/(wx^2 + wy^2 + wz^2)^2 - (2*wx*(cos((wx^2 + wy^2 + wz^2)^(1/2)) - 1))/(wx^2 + wy^2 + wz^2) + (wx*sin((wx^2 + wy^2 + wz^2)^(1/2))*(wx^2 + wy^2))/(wx^2 + wy^2 + wz^2)^(3/2)) - f*((wx^2*sin((wx^2 + wy^2 + wz^2)^(1/2)))/(wx^2 + wy^2 + wz^2)^(3/2) - (wx^2*cos((wx^2 + wy^2 + wz^2)^(1/2)))/(wx^2 + wy^2 + wz^2) - sin((wx^2 + wy^2 + wz^2)^(1/2))/(wx^2 + wy^2 + wz^2)^(1/2) + (2*wx*wy*wz*(cos((wx^2 + wy^2 + wz^2)^(1/2)) - 1))/(wx^2 + wy^2 + wz^2)^2 + (wx*wy*wz*sin((wx^2 + wy^2 + wz^2)^(1/2)))/(wx^2 + wy^2 + wz^2)^(3/2))) - X*(f*((wx*wz*cos((wx^2 + wy^2 + wz^2)^(1/2)))/(wx^2 + wy^2 + wz^2) - (wy*(cos((wx^2 + wy^2 + wz^2)^(1/2)) - 1))/(wx^2 + wy^2 + wz^2) - (wx*wz*sin((wx^2 + wy^2 + wz^2)^(1/2)))/(wx^2 + wy^2 + wz^2)^(3/2) + (wx^2*wy*sin((wx^2 + wy^2 + wz^2)^(1/2)))/(wx^2 + wy^2 + wz^2)^(3/2) + (2*wx^2*wy*(cos((wx^2 + wy^2 + wz^2)^(1/2)) - 1))/(wx^2 + wy^2 + wz^2)^2) + v0*((wx*wy*sin((wx^2 + wy^2 + wz^2)^(1/2)))/(wx^2 + wy^2 + wz^2)^(3/2) - (wx*wy*cos((wx^2 + wy^2 + wz^2)^(1/2)))/(wx^2 + wy^2 + wz^2) - (wz*(cos((wx^2 + wy^2 + wz^2)^(1/2)) - 1))/(wx^2 + wy^2 + wz^2) + (wx^2*wz*sin((wx^2 + wy^2 + wz^2)^(1/2)))/(wx^2 + wy^2 + wz^2)^(3/2) + (2*wx^2*wz*(cos((wx^2 + wy^2 + wz^2)^(1/2)) - 1))/(wx^2 + wy^2 + wz^2)^2)))/(tz - X*((wy*sin((wx^2 + wy^2 + wz^2)^(1/2)))/(wx^2 + wy^2 + wz^2)^(1/2) + (wx*wz*(cos((wx^2 + wy^2 + wz^2)^(1/2)) - 1))/(wx^2 + wy^2 + wz^2)) + Y*((wx*sin((wx^2 + wy^2 + wz^2)^(1/2)))/(wx^2 + wy^2 + wz^2)^(1/2) - (wy*wz*(cos((wx^2 + wy^2 + wz^2)^(1/2)) - 1))/(wx^2 + wy^2 + wz^2)) + Z*(((cos((wx^2 + wy^2 + wz^2)^(1/2)) - 1)*(wx^2 + wy^2))/(wx^2 + wy^2 + wz^2) + 1)) - ((X*((wx*wy*sin((wx^2 + wy^2 + wz^2)^(1/2)))/(wx^2 + wy^2 + wz^2)^(3/2) - (wx*wy*cos((wx^2 + wy^2 + wz^2)^(1/2)))/(wx^2 + wy^2 + wz^2) - (wz*(cos((wx^2 + wy^2 + wz^2)^(1/2)) - 1))/(wx^2 + wy^2 + wz^2) + (wx^2*wz*sin((wx^2 + wy^2 + wz^2)^(1/2)))/(wx^2 + wy^2 + wz^2)^(3/2) + (2*wx^2*wz*(cos((wx^2 + wy^2 + wz^2)^(1/2)) - 1))/(wx^2 + wy^2 + wz^2)^2) - Z*((2*wx*(cos((wx^2 + wy^2 + wz^2)^(1/2)) - 1)*(wx^2 + wy^2))/(wx^2 + wy^2 + wz^2)^2 - (2*wx*(cos((wx^2 + wy^2 + wz^2)^(1/2)) - 1))/(wx^2 + wy^2 + wz^2) + (wx*sin((wx^2 + wy^2 + wz^2)^(1/2))*(wx^2 + wy^2))/(wx^2 + wy^2 + wz^2)^(3/2)) + Y*(sin((wx^2 + wy^2 + wz^2)^(1/2))/(wx^2 + wy^2 + wz^2)^(1/2) + (wx^2*cos((wx^2 + wy^2 + wz^2)^(1/2)))/(wx^2 + wy^2 + wz^2) - (wx^2*sin((wx^2 + wy^2 + wz^2)^(1/2)))/(wx^2 + wy^2 + wz^2)^(3/2) + (2*wx*wy*wz*(cos((wx^2 + wy^2 + wz^2)^(1/2)) - 1))/(wx^2 + wy^2 + wz^2)^2 + (wx*wy*wz*sin((wx^2 + wy^2 + wz^2)^(1/2)))/(wx^2 + wy^2 + wz^2)^(3/2)))*(f*ty + X*(f*((wz*sin((wx^2 + wy^2 + wz^2)^(1/2)))/(wx^2 + wy^2 + wz^2)^(1/2) - (wx*wy*(cos((wx^2 + wy^2 + wz^2)^(1/2)) - 1))/(wx^2 + wy^2 + wz^2)) - v0*((wy*sin((wx^2 + wy^2 + wz^2)^(1/2)))/(wx^2 + wy^2 + wz^2)^(1/2) + (wx*wz*(cos((wx^2 + wy^2 + wz^2)^(1/2)) - 1))/(wx^2 + wy^2 + wz^2))) + tz*v0 + Y*(v0*((wx*sin((wx^2 + wy^2 + wz^2)^(1/2)))/(wx^2 + wy^2 + wz^2)^(1/2) - (wy*wz*(cos((wx^2 + wy^2 + wz^2)^(1/2)) - 1))/(wx^2 + wy^2 + wz^2)) + f*(((cos((wx^2 + wy^2 + wz^2)^(1/2)) - 1)*(wx^2 + wz^2))/(wx^2 + wy^2 + wz^2) + 1)) - Z*(f*((wx*sin((wx^2 + wy^2 + wz^2)^(1/2)))/(wx^2 + wy^2 + wz^2)^(1/2) + (wy*wz*(cos((wx^2 + wy^2 + wz^2)^(1/2)) - 1))/(wx^2 + wy^2 + wz^2)) - v0*(((cos((wx^2 + wy^2 + wz^2)^(1/2)) - 1)*(wx^2 + wy^2))/(wx^2 + wy^2 + wz^2) + 1))))/(tz - X*((wy*sin((wx^2 + wy^2 + wz^2)^(1/2)))/(wx^2 + wy^2 + wz^2)^(1/2) + (wx*wz*(cos((wx^2 + wy^2 + wz^2)^(1/2)) - 1))/(wx^2 + wy^2 + wz^2)) + Y*((wx*sin((wx^2 + wy^2 + wz^2)^(1/2)))/(wx^2 + wy^2 + wz^2)^(1/2) - (wy*wz*(cos((wx^2 + wy^2 + wz^2)^(1/2)) - 1))/(wx^2 + wy^2 + wz^2)) + Z*(((cos((wx^2 + wy^2 + wz^2)^(1/2)) - 1)*(wx^2 + wy^2))/(wx^2 + wy^2 + wz^2) + 1))^2, - (Y*(f*((2*wy*(cos((wx^2 + wy^2 + wz^2)^(1/2)) - 1)*(wx^2 + wz^2))/(wx^2 + wy^2 + wz^2)^2 + (wy*sin((wx^2 + wy^2 + wz^2)^(1/2))*(wx^2 + wz^2))/(wx^2 + wy^2 + wz^2)^(3/2)) - v0*((wx*wy*cos((wx^2 + wy^2 + wz^2)^(1/2)))/(wx^2 + wy^2 + wz^2) - (wz*(cos((wx^2 + wy^2 + wz^2)^(1/2)) - 1))/(wx^2 + wy^2 + wz^2) - (wx*wy*sin((wx^2 + wy^2 + wz^2)^(1/2)))/(wx^2 + wy^2 + wz^2)^(3/2) + (wy^2*wz*sin((wx^2 + wy^2 + wz^2)^(1/2)))/(wx^2 + wy^2 + wz^2)^(3/2) + (2*wy^2*wz*(cos((wx^2 + wy^2 + wz^2)^(1/2)) - 1))/(wx^2 + wy^2 + wz^2)^2)) - X*(f*((wy*wz*cos((wx^2 + wy^2 + wz^2)^(1/2)))/(wx^2 + wy^2 + wz^2) - (wx*(cos((wx^2 + wy^2 + wz^2)^(1/2)) - 1))/(wx^2 + wy^2 + wz^2) - (wy*wz*sin((wx^2 + wy^2 + wz^2)^(1/2)))/(wx^2 + wy^2 + wz^2)^(3/2) + (wx*wy^2*sin((wx^2 + wy^2 + wz^2)^(1/2)))/(wx^2 + wy^2 + wz^2)^(3/2) + (2*wx*wy^2*(cos((wx^2 + wy^2 + wz^2)^(1/2)) - 1))/(wx^2 + wy^2 + wz^2)^2) + v0*((wy^2*sin((wx^2 + wy^2 + wz^2)^(1/2)))/(wx^2 + wy^2 + wz^2)^(3/2) - (wy^2*cos((wx^2 + wy^2 + wz^2)^(1/2)))/(wx^2 + wy^2 + wz^2) - sin((wx^2 + wy^2 + wz^2)^(1/2))/(wx^2 + wy^2 + wz^2)^(1/2) + (2*wx*wy*wz*(cos((wx^2 + wy^2 + wz^2)^(1/2)) - 1))/(wx^2 + wy^2 + wz^2)^2 + (wx*wy*wz*sin((wx^2 + wy^2 + wz^2)^(1/2)))/(wx^2 + wy^2 + wz^2)^(3/2))) + Z*(v0*((2*wy*(cos((wx^2 + wy^2 + wz^2)^(1/2)) - 1)*(wx^2 + wy^2))/(wx^2 + wy^2 + wz^2)^2 - (2*wy*(cos((wx^2 + wy^2 + wz^2)^(1/2)) - 1))/(wx^2 + wy^2 + wz^2) + (wy*sin((wx^2 + wy^2 + wz^2)^(1/2))*(wx^2 + wy^2))/(wx^2 + wy^2 + wz^2)^(3/2)) - f*((wx*wy*sin((wx^2 + wy^2 + wz^2)^(1/2)))/(wx^2 + wy^2 + wz^2)^(3/2) - (wx*wy*cos((wx^2 + wy^2 + wz^2)^(1/2)))/(wx^2 + wy^2 + wz^2) - (wz*(cos((wx^2 + wy^2 + wz^2)^(1/2)) - 1))/(wx^2 + wy^2 + wz^2) + (wy^2*wz*sin((wx^2 + wy^2 + wz^2)^(1/2)))/(wx^2 + wy^2 + wz^2)^(3/2) + (2*wy^2*wz*(cos((wx^2 + wy^2 + wz^2)^(1/2)) - 1))/(wx^2 + wy^2 + wz^2)^2)))/(tz - X*((wy*sin((wx^2 + wy^2 + wz^2)^(1/2)))/(wx^2 + wy^2 + wz^2)^(1/2) + (wx*wz*(cos((wx^2 + wy^2 + wz^2)^(1/2)) - 1))/(wx^2 + wy^2 + wz^2)) + Y*((wx*sin((wx^2 + wy^2 + wz^2)^(1/2)))/(wx^2 + wy^2 + wz^2)^(1/2) - (wy*wz*(cos((wx^2 + wy^2 + wz^2)^(1/2)) - 1))/(wx^2 + wy^2 + wz^2)) + Z*(((cos((wx^2 + wy^2 + wz^2)^(1/2)) - 1)*(wx^2 + wy^2))/(wx^2 + wy^2 + wz^2) + 1)) - ((Y*((wx*wy*cos((wx^2 + wy^2 + wz^2)^(1/2)))/(wx^2 + wy^2 + wz^2) - (wz*(cos((wx^2 + wy^2 + wz^2)^(1/2)) - 1))/(wx^2 + wy^2 + wz^2) - (wx*wy*sin((wx^2 + wy^2 + wz^2)^(1/2)))/(wx^2 + wy^2 + wz^2)^(3/2) + (wy^2*wz*sin((wx^2 + wy^2 + wz^2)^(1/2)))/(wx^2 + wy^2 + wz^2)^(3/2) + (2*wy^2*wz*(cos((wx^2 + wy^2 + wz^2)^(1/2)) - 1))/(wx^2 + wy^2 + wz^2)^2) - Z*((2*wy*(cos((wx^2 + wy^2 + wz^2)^(1/2)) - 1)*(wx^2 + wy^2))/(wx^2 + wy^2 + wz^2)^2 - (2*wy*(cos((wx^2 + wy^2 + wz^2)^(1/2)) - 1))/(wx^2 + wy^2 + wz^2) + (wy*sin((wx^2 + wy^2 + wz^2)^(1/2))*(wx^2 + wy^2))/(wx^2 + wy^2 + wz^2)^(3/2)) + X*((wy^2*sin((wx^2 + wy^2 + wz^2)^(1/2)))/(wx^2 + wy^2 + wz^2)^(3/2) - (wy^2*cos((wx^2 + wy^2 + wz^2)^(1/2)))/(wx^2 + wy^2 + wz^2) - sin((wx^2 + wy^2 + wz^2)^(1/2))/(wx^2 + wy^2 + wz^2)^(1/2) + (2*wx*wy*wz*(cos((wx^2 + wy^2 + wz^2)^(1/2)) - 1))/(wx^2 + wy^2 + wz^2)^2 + (wx*wy*wz*sin((wx^2 + wy^2 + wz^2)^(1/2)))/(wx^2 + wy^2 + wz^2)^(3/2)))*(f*ty + X*(f*((wz*sin((wx^2 + wy^2 + wz^2)^(1/2)))/(wx^2 + wy^2 + wz^2)^(1/2) - (wx*wy*(cos((wx^2 + wy^2 + wz^2)^(1/2)) - 1))/(wx^2 + wy^2 + wz^2)) - v0*((wy*sin((wx^2 + wy^2 + wz^2)^(1/2)))/(wx^2 + wy^2 + wz^2)^(1/2) + (wx*wz*(cos((wx^2 + wy^2 + wz^2)^(1/2)) - 1))/(wx^2 + wy^2 + wz^2))) + tz*v0 + Y*(v0*((wx*sin((wx^2 + wy^2 + wz^2)^(1/2)))/(wx^2 + wy^2 + wz^2)^(1/2) - (wy*wz*(cos((wx^2 + wy^2 + wz^2)^(1/2)) - 1))/(wx^2 + wy^2 + wz^2)) + f*(((cos((wx^2 + wy^2 + wz^2)^(1/2)) - 1)*(wx^2 + wz^2))/(wx^2 + wy^2 + wz^2) + 1)) - Z*(f*((wx*sin((wx^2 + wy^2 + wz^2)^(1/2)))/(wx^2 + wy^2 + wz^2)^(1/2) + (wy*wz*(cos((wx^2 + wy^2 + wz^2)^(1/2)) - 1))/(wx^2 + wy^2 + wz^2)) - v0*(((cos((wx^2 + wy^2 + wz^2)^(1/2)) - 1)*(wx^2 + wy^2))/(wx^2 + wy^2 + wz^2) + 1))))/(tz - X*((wy*sin((wx^2 + wy^2 + wz^2)^(1/2)))/(wx^2 + wy^2 + wz^2)^(1/2) + (wx*wz*(cos((wx^2 + wy^2 + wz^2)^(1/2)) - 1))/(wx^2 + wy^2 + wz^2)) + Y*((wx*sin((wx^2 + wy^2 + wz^2)^(1/2)))/(wx^2 + wy^2 + wz^2)^(1/2) - (wy*wz*(cos((wx^2 + wy^2 + wz^2)^(1/2)) - 1))/(wx^2 + wy^2 + wz^2)) + Z*(((cos((wx^2 + wy^2 + wz^2)^(1/2)) - 1)*(wx^2 + wy^2))/(wx^2 + wy^2 + wz^2) + 1))^2, - (Z*(v0*((2*wz*(cos((wx^2 + wy^2 + wz^2)^(1/2)) - 1)*(wx^2 + wy^2))/(wx^2 + wy^2 + wz^2)^2 + (wz*sin((wx^2 + wy^2 + wz^2)^(1/2))*(wx^2 + wy^2))/(wx^2 + wy^2 + wz^2)^(3/2)) - f*((wx*wz*sin((wx^2 + wy^2 + wz^2)^(1/2)))/(wx^2 + wy^2 + wz^2)^(3/2) - (wx*wz*cos((wx^2 + wy^2 + wz^2)^(1/2)))/(wx^2 + wy^2 + wz^2) - (wy*(cos((wx^2 + wy^2 + wz^2)^(1/2)) - 1))/(wx^2 + wy^2 + wz^2) + (wy*wz^2*sin((wx^2 + wy^2 + wz^2)^(1/2)))/(wx^2 + wy^2 + wz^2)^(3/2) + (2*wy*wz^2*(cos((wx^2 + wy^2 + wz^2)^(1/2)) - 1))/(wx^2 + wy^2 + wz^2)^2)) - X*(v0*((wy*wz*sin((wx^2 + wy^2 + wz^2)^(1/2)))/(wx^2 + wy^2 + wz^2)^(3/2) - (wy*wz*cos((wx^2 + wy^2 + wz^2)^(1/2)))/(wx^2 + wy^2 + wz^2) - (wx*(cos((wx^2 + wy^2 + wz^2)^(1/2)) - 1))/(wx^2 + wy^2 + wz^2) + (wx*wz^2*sin((wx^2 + wy^2 + wz^2)^(1/2)))/(wx^2 + wy^2 + wz^2)^(3/2) + (2*wx*wz^2*(cos((wx^2 + wy^2 + wz^2)^(1/2)) - 1))/(wx^2 + wy^2 + wz^2)^2) + f*(sin((wx^2 + wy^2 + wz^2)^(1/2))/(wx^2 + wy^2 + wz^2)^(1/2) + (wz^2*cos((wx^2 + wy^2 + wz^2)^(1/2)))/(wx^2 + wy^2 + wz^2) - (wz^2*sin((wx^2 + wy^2 + wz^2)^(1/2)))/(wx^2 + wy^2 + wz^2)^(3/2) + (2*wx*wy*wz*(cos((wx^2 + wy^2 + wz^2)^(1/2)) - 1))/(wx^2 + wy^2 + wz^2)^2 + (wx*wy*wz*sin((wx^2 + wy^2 + wz^2)^(1/2)))/(wx^2 + wy^2 + wz^2)^(3/2))) + Y*(f*((2*wz*(cos((wx^2 + wy^2 + wz^2)^(1/2)) - 1)*(wx^2 + wz^2))/(wx^2 + wy^2 + wz^2)^2 - (2*wz*(cos((wx^2 + wy^2 + wz^2)^(1/2)) - 1))/(wx^2 + wy^2 + wz^2) + (wz*sin((wx^2 + wy^2 + wz^2)^(1/2))*(wx^2 + wz^2))/(wx^2 + wy^2 + wz^2)^(3/2)) - v0*((wx*wz*cos((wx^2 + wy^2 + wz^2)^(1/2)))/(wx^2 + wy^2 + wz^2) - (wy*(cos((wx^2 + wy^2 + wz^2)^(1/2)) - 1))/(wx^2 + wy^2 + wz^2) - (wx*wz*sin((wx^2 + wy^2 + wz^2)^(1/2)))/(wx^2 + wy^2 + wz^2)^(3/2) + (wy*wz^2*sin((wx^2 + wy^2 + wz^2)^(1/2)))/(wx^2 + wy^2 + wz^2)^(3/2) + (2*wy*wz^2*(cos((wx^2 + wy^2 + wz^2)^(1/2)) - 1))/(wx^2 + wy^2 + wz^2)^2)))/(tz - X*((wy*sin((wx^2 + wy^2 + wz^2)^(1/2)))/(wx^2 + wy^2 + wz^2)^(1/2) + (wx*wz*(cos((wx^2 + wy^2 + wz^2)^(1/2)) - 1))/(wx^2 + wy^2 + wz^2)) + Y*((wx*sin((wx^2 + wy^2 + wz^2)^(1/2)))/(wx^2 + wy^2 + wz^2)^(1/2) - (wy*wz*(cos((wx^2 + wy^2 + wz^2)^(1/2)) - 1))/(wx^2 + wy^2 + wz^2)) + Z*(((cos((wx^2 + wy^2 + wz^2)^(1/2)) - 1)*(wx^2 + wy^2))/(wx^2 + wy^2 + wz^2) + 1)) - ((X*((wy*wz*sin((wx^2 + wy^2 + wz^2)^(1/2)))/(wx^2 + wy^2 + wz^2)^(3/2) - (wy*wz*cos((wx^2 + wy^2 + wz^2)^(1/2)))/(wx^2 + wy^2 + wz^2) - (wx*(cos((wx^2 + wy^2 + wz^2)^(1/2)) - 1))/(wx^2 + wy^2 + wz^2) + (wx*wz^2*sin((wx^2 + wy^2 + wz^2)^(1/2)))/(wx^2 + wy^2 + wz^2)^(3/2) + (2*wx*wz^2*(cos((wx^2 + wy^2 + wz^2)^(1/2)) - 1))/(wx^2 + wy^2 + wz^2)^2) - Z*((2*wz*(cos((wx^2 + wy^2 + wz^2)^(1/2)) - 1)*(wx^2 + wy^2))/(wx^2 + wy^2 + wz^2)^2 + (wz*sin((wx^2 + wy^2 + wz^2)^(1/2))*(wx^2 + wy^2))/(wx^2 + wy^2 + wz^2)^(3/2)) + Y*((wx*wz*cos((wx^2 + wy^2 + wz^2)^(1/2)))/(wx^2 + wy^2 + wz^2) - (wy*(cos((wx^2 + wy^2 + wz^2)^(1/2)) - 1))/(wx^2 + wy^2 + wz^2) - (wx*wz*sin((wx^2 + wy^2 + wz^2)^(1/2)))/(wx^2 + wy^2 + wz^2)^(3/2) + (wy*wz^2*sin((wx^2 + wy^2 + wz^2)^(1/2)))/(wx^2 + wy^2 + wz^2)^(3/2) + (2*wy*wz^2*(cos((wx^2 + wy^2 + wz^2)^(1/2)) - 1))/(wx^2 + wy^2 + wz^2)^2))*(f*ty + X*(f*((wz*sin((wx^2 + wy^2 + wz^2)^(1/2)))/(wx^2 + wy^2 + wz^2)^(1/2) - (wx*wy*(cos((wx^2 + wy^2 + wz^2)^(1/2)) - 1))/(wx^2 + wy^2 + wz^2)) - v0*((wy*sin((wx^2 + wy^2 + wz^2)^(1/2)))/(wx^2 + wy^2 + wz^2)^(1/2) + (wx*wz*(cos((wx^2 + wy^2 + wz^2)^(1/2)) - 1))/(wx^2 + wy^2 + wz^2))) + tz*v0 + Y*(v0*((wx*sin((wx^2 + wy^2 + wz^2)^(1/2)))/(wx^2 + wy^2 + wz^2)^(1/2) - (wy*wz*(cos((wx^2 + wy^2 + wz^2)^(1/2)) - 1))/(wx^2 + wy^2 + wz^2)) + f*(((cos((wx^2 + wy^2 + wz^2)^(1/2)) - 1)*(wx^2 + wz^2))/(wx^2 + wy^2 + wz^2) + 1)) - Z*(f*((wx*sin((wx^2 + wy^2 + wz^2)^(1/2)))/(wx^2 + wy^2 + wz^2)^(1/2) + (wy*wz*(cos((wx^2 + wy^2 + wz^2)^(1/2)) - 1))/(wx^2 + wy^2 + wz^2)) - v0*(((cos((wx^2 + wy^2 + wz^2)^(1/2)) - 1)*(wx^2 + wy^2))/(wx^2 + wy^2 + wz^2) + 1))))/(tz - X*((wy*sin((wx^2 + wy^2 + wz^2)^(1/2)))/(wx^2 + wy^2 + wz^2)^(1/2) + (wx*wz*(cos((wx^2 + wy^2 + wz^2)^(1/2)) - 1))/(wx^2 + wy^2 + wz^2)) + Y*((wx*sin((wx^2 + wy^2 + wz^2)^(1/2)))/(wx^2 + wy^2 + wz^2)^(1/2) - (wy*wz*(cos((wx^2 + wy^2 + wz^2)^(1/2)) - 1))/(wx^2 + wy^2 + wz^2)) + Z*(((cos((wx^2 + wy^2 + wz^2)^(1/2)) - 1)*(wx^2 + wy^2))/(wx^2 + wy^2 + wz^2) + 1))^2, 0, f/(tz - X*((wy*sin((wx^2 + wy^2 + wz^2)^(1/2)))/(wx^2 + wy^2 + wz^2)^(1/2) + (wx*wz*(cos((wx^2 + wy^2 + wz^2)^(1/2)) - 1))/(wx^2 + wy^2 + wz^2)) + Y*((wx*sin((wx^2 + wy^2 + wz^2)^(1/2)))/(wx^2 + wy^2 + wz^2)^(1/2) - (wy*wz*(cos((wx^2 + wy^2 + wz^2)^(1/2)) - 1))/(wx^2 + wy^2 + wz^2)) + Z*(((cos((wx^2 + wy^2 + wz^2)^(1/2)) - 1)*(wx^2 + wy^2))/(wx^2 + wy^2 + wz^2) + 1)), v0/(tz - X*((wy*sin((wx^2 + wy^2 + wz^2)^(1/2)))/(wx^2 + wy^2 + wz^2)^(1/2) + (wx*wz*(cos((wx^2 + wy^2 + wz^2)^(1/2)) - 1))/(wx^2 + wy^2 + wz^2)) + Y*((wx*sin((wx^2 + wy^2 + wz^2)^(1/2)))/(wx^2 + wy^2 + wz^2)^(1/2) - (wy*wz*(cos((wx^2 + wy^2 + wz^2)^(1/2)) - 1))/(wx^2 + wy^2 + wz^2)) + Z*(((cos((wx^2 + wy^2 + wz^2)^(1/2)) - 1)*(wx^2 + wy^2))/(wx^2 + wy^2 + wz^2) + 1)) - (f*ty + X*(f*((wz*sin((wx^2 + wy^2 + wz^2)^(1/2)))/(wx^2 + wy^2 + wz^2)^(1/2) - (wx*wy*(cos((wx^2 + wy^2 + wz^2)^(1/2)) - 1))/(wx^2 + wy^2 + wz^2)) - v0*((wy*sin((wx^2 + wy^2 + wz^2)^(1/2)))/(wx^2 + wy^2 + wz^2)^(1/2) + (wx*wz*(cos((wx^2 + wy^2 + wz^2)^(1/2)) - 1))/(wx^2 + wy^2 + wz^2))) + tz*v0 + Y*(v0*((wx*sin((wx^2 + wy^2 + wz^2)^(1/2)))/(wx^2 + wy^2 + wz^2)^(1/2) - (wy*wz*(cos((wx^2 + wy^2 + wz^2)^(1/2)) - 1))/(wx^2 + wy^2 + wz^2)) + f*(((cos((wx^2 + wy^2 + wz^2)^(1/2)) - 1)*(wx^2 + wz^2))/(wx^2 + wy^2 + wz^2) + 1)) - Z*(f*((wx*sin((wx^2 + wy^2 + wz^2)^(1/2)))/(wx^2 + wy^2 + wz^2)^(1/2) + (wy*wz*(cos((wx^2 + wy^2 + wz^2)^(1/2)) - 1))/(wx^2 + wy^2 + wz^2)) - v0*(((cos((wx^2 + wy^2 + wz^2)^(1/2)) - 1)*(wx^2 + wy^2))/(wx^2 + wy^2 + wz^2) + 1)))/(tz - X*((wy*sin((wx^2 + wy^2 + wz^2)^(1/2)))/(wx^2 + wy^2 + wz^2)^(1/2) + (wx*wz*(cos((wx^2 + wy^2 + wz^2)^(1/2)) - 1))/(wx^2 + wy^2 + wz^2)) + Y*((wx*sin((wx^2 + wy^2 + wz^2)^(1/2)))/(wx^2 + wy^2 + wz^2)^(1/2) - (wy*wz*(cos((wx^2 + wy^2 + wz^2)^(1/2)) - 1))/(wx^2 + wy^2 + wz^2)) + Z*(((cos((wx^2 + wy^2 + wz^2)^(1/2)) - 1)*(wx^2 + wy^2))/(wx^2 + wy^2 + wz^2) + 1))^2];
%     omega=  [0 -wz wy;
%          wz 0 -wx;
%         -wy wx 0;]
%     e = 2 ;
%     v1 = (wx * omega +   cross(omega,(eye(3)-rotationMatrix)*e))/norm(omega)
%     v1 * rotationMatrix
% %


