function [ result D s] = image2world_ex2(image ,f,d, rotationMatrix,translationVector, cameraParams)
%

[r,face] = read_ply('teabox.ply');
face = face + 1 ;
% implement C = ?Q^?1 * q 

 P = cameraMatrix(cameraParams,rotationMatrix,translationVector);
 p = P';
 Q = p(:,1:3); q = p(:,4);
 Orig = -inv(Q) * q ;


% I = rgb2gray(image) ;
% 
% peak_thresh = 2;
% [f,d] = vl_sift(single(I),'PeakThresh', peak_thresh) ;
% extract the image features



% for each point use TriangleRayIntersection to check the intersection
% if there is the intersection, convert them to world coordinate system 
result=[];
s=[];
D = [] ;
for i=1:size(f,2)
    ray = inv( Q) * [f(1:2,i)' 1]';
    [intersect, t, u, v, xcoor] = TriangleRayIntersection (Orig', ray', r(face(:,1),:), r(face(:,2),:), r(face(:,3),:));
    as = find(intersect == 1);
    if ~isempty(as)
        s = [s f(1:2,i)];
        D = [D d(:,i)];
        project = Orig + ray * min(t(as)) ;
        result =[result ; project'];
    end    
        
end  
% figure(7)
% imshow(I);
% hold on ;
% 
% scatter(s(1,:),s(2,:));
% a = [1:10]'; b = num2str(a); c = cellstr(b);
% text(s(1,1:10),s(2,1:10), c);
% 
% see =[result ;Orig'];
% 
% figure(10)
% scatter3(result(:,1),result(:,2),result(:,3),1);
% a = [1:10]'; b = num2str(a); c = cellstr(b);
% text(result(1:10,1),result(1:10,2),result(1:10,3), c);
% 
% hold on ;
% trisurf(face, r(:,1),r(:,2),r(:,3),'FaceAlpha', 0.5);
% grid on ;
end

