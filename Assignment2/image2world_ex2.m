function [ result,d ] = image2world(r, face,image, rotationMatrix,translationVector, cameraParams)
%
% implement C = ?Q^?1 * q 
 P = cameraMatrix(cameraParams,rotationMatrix,translationVector);
 p = P';
 Q = p(:,1:3); q = p(:,4);
 Orig = -inv(Q) * q ;


I = rgb2gray(image) ;
% imshow(I);
% hold on ;

[f,d] = vl_sift(single(I)) ;
% extract the image features


% h1 = vl_plotframe(f) ;
% set(h1,'color','b','linewidth',3) ;

% for each point use TriangleRayIntersection to check the intersection
% if there is the intersection, convert them to world coordinate system 
result=[];
for i=1:size(f,2)
    ray =   inv( Q) * [f(1:2,i)' 1]';
    [intersect, t, u, v, xcoor] = TriangleRayIntersection (Orig', ray', r(face(:,1),:), r(face(:,2),:), r(face(:,3),:));
    as = find(intersect == 1);
    if ~isempty(as)
        project = Orig + ray * min(t(as)) ;
        result =[result ; project'];
    end    
        
end  
% see =[result ;Orig'];
% 
% figure(3)
% scatter3(see(:,1),see(:,2),see(:,3));
% hold on ;
% trisurf(face, r(:,1),r(:,2),r(:,3),'FaceAlpha', 0.5);
% grid on ;
end

