%     dx_remesh = (c/fmax)/ppw_remesh;
%     tprune = round((ppw/ppw_remesh)/CFL);
%     indxs_t = 1:tprune:timesteps;
% 
%     [XX_out,YY_out] = uniformGrid(dx_remesh,xminmax,yminmax);
%     grid = zip2(XX_out,YY_out);
%     
%     p_srcs_remesh = zeros(size(p_srcs,1), length(XX_out(:)), length(indxs_t));
%     v_srcs_remesh = zeros(size(p_srcs,1), length(XX_out(:)), 2, length(indxs_t));
%     
%     for i=1:size(p_srcs,1)
%         [p,vx,vy] = remeshFields(X2D,Y2D,XX_out,YY_out,...
%             squeeze(p_srcs(i,:,:)),squeeze(vx_srcs(i,:,:)),squeeze(vy_srcs(i,:,:)));
%         p_srcs_remesh(i,:,:) = p(:,indxs_t);
%         v = zip2(vx,vy);
%         v_srcs_remesh(i,:,:,:) = v(:,:,indxs_t);
%     end    

function [XX_ip,YY_ip] = uniformGrid(dx,xminmax,yminmax)
    %t = linspace(0,tmax,ceil(tmax/dt));
    x = linspace(xminmax(1),xminmax(2),ceil((xminmax(2)-xminmax(1))/dx));
    y = linspace(xminmax(1),xminmax(2),ceil((yminmax(2)-yminmax(1))/dx));
    
    [XX_ip, YY_ip] = meshgrid(x,y);    
end

function [p,vx,vy] = remeshFields(X2D,Y2D,XX_out,YY_out,p,vx,vy)
    for i=1:size(p,2)
        p(:,i) = griddata(X2D(:),Y2D(:),p(:,i), XX_out(:),YY_out(:),'cubic');
        vx(:,i) = griddata(X2D(:),Y2D(:),vx(:,i),XX_out(:),YY_out(:),'cubic');
        vy(:,i) = griddata(X2D(:),Y2D(:),vy(:,i),XX_out(:),YY_out(:),'cubic');
    end
end