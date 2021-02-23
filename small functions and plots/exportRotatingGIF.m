function [] = exportRotatingGIF(f,fullExportPath)
%EXPORTROTATINGGIF creates a rotating gif from figure f and exports it to
%gif at fullExportPath (fullExportPath should end with '.gif')
%   taken from code i used in
%   '/media/sil2/Literature/Projects/corplex/progress reports/figs for PR
%   201021/scriptForPlots.m', don't rememeber where I got it from originally

initial_az=-30;
initial_el=20;
az = initial_az;
el = initial_el;
view([az,el])
degStep = 1;
detlaT = 0.1;
fCount = 71;
h = getframe(gcf);
[im,map] = rgb2ind(h.cdata,256,'nodither');
im(1,1,1,fCount) = 0;
k = 1;
% % % % % spin 45°
% % % % for i = 0:-degStep:-45
% % % %   az = i;
% % % %   view([az,el])
% % % %   h = getframe(f);
% % % %   im(:,:,1,k) = rgb2ind(h.cdata,map,'nodither');
% % % %   k = k + 1;
% % % % end
% % % % % tilt down
% % % % for i = 90:-degStep:15
% % % %   el = i;
% % % %   view([az,el])
% % % %   h = getframe(f);
% % % %   im(:,:,1,k) = rgb2ind(h.cdata,map,'nodither');
% % % %   k = k + 1;
% % % % end
% spin left
for i = az:-degStep:-90
  az = i;
  view([az,el])
  h = getframe(f);
  im(:,:,1,k) = rgb2ind(h.cdata,map,'nodither');
  k = k + 1;
end
% spin right
for i = az:degStep:0
  az = i;
  view([az,el])
  h = getframe(f);
  im(:,:,1,k) = rgb2ind(h.cdata,map,'nodither');
  k = k + 1;
end
%spin to initial
for i = az:-degStep:initial_az
  az = i;
  view([az,el])
  h = getframe(f);
  im(:,:,1,k) = rgb2ind(h.cdata,map,'nodither');
  k = k + 1;
end

% tilt up to original
% % % % for i = el:degStep:90
% % % %   el = i;
% % % %   view([az,el])
% % % %   h = getframe(f);
% % % %   im(:,:,1,k) = rgb2ind(h.cdata,map,'nodither');
% % % %   k = k + 1;
% % % % end
imwrite(im,map,fullExportPath,'DelayTime',detlaT,'LoopCount',inf)
imwrite(im,map,fullExportPath,'DelayTime',detlaT,'LoopCount',inf)



end

