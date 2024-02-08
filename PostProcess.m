function [preImage, leftPreImage]= PostProcess(preImage,leftPreImage,stage)
   [~,sizeY] = size(preImage);
%    sizeY = 16;
   vectorLeft = leftPreImage(:,sizeY-stage);
   vectorPre = preImage(:,stage);
   tempMatrix = zeros(sizeY,stage*2);
   for i=1:sizeY
       linspaceVector = linspace(vectorLeft(i),vectorPre(i),stage*2+2);
       tempMatrix(i,:) = linspaceVector(2:stage*2+1);
   end
   leftPreImage(:,sizeY-stage+1:sizeY) = tempMatrix(:,1:stage);
   preImage(:,1:stage) = tempMatrix(:,stage+1:stage*2);
end