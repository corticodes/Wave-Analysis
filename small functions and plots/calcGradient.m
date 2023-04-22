function [grad_x,grad_y] = calcGradient(data)
%CALCGRADIENT calcs the spatial gradient of data
%   data is framHeihtXframeWidth
%   If padded=0 than grad_x,grad_y is (frameHeight-1)X(frameWidth-1)
%   else it is padded with zeros to be size(data)

grad_x=diff(data,1,2);
grad_y=diff(data,1,1);

% if padded
grad_x=padarray(grad_x,[0 1],0,'post');
grad_y=padarray(grad_y,[1 0],0,'post');
% end

end

