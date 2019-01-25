function e = projection_err(cameraParams,theta,M,m)
    rotationMatrix = rotationVectorToMatrix(theta(1:3));
    estimate_point = worldToImage(cameraParams,rotationMatrix,theta(4:6),M);
%     plot(estimate_point(:,1),estimate_point(:,2),'r*','Color','r');
%     hold on ;
    err = [estimate_point - m];
    
    
    e = [err(:,1);err(:,2)];
    
   
%     e = abs(err(:,1))+ abs(err(:,2));
end

