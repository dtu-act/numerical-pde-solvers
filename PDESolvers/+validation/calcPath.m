function [P1_path, P2_path, XY_path] = calcPath(P1,P2,X,Y,eval_on_path)
    assert(size(P1,1) == size(P2,1) && size(P1,2) == size(P2,2));
    
    NM = size(P1);
    
    switch eval_on_path
        case 'axial-vertical-center'
            indx = (NM(2)-1)/2+1;
            P1_path = P1(:,indx);
            P2_path = P2(:,indx);
            X_path = X(:,indx)';
            Y_path = Y(indx,:);
        case 'axial-horizontal-center'
            indx = (NM(1)-1)/2+1;
            P1_path = P1(indx,:);
            P2_path = P2(indx,:);
            X_path = X(indx,:)';
            Y_path = Y(:,indx);
        case 'diagonal'
            P1_path = diag(P1);
            P2_path = diag(P2);
            X_path = diag(X);
            Y_path = diag(Y);
        otherwise
            error('path not supported')
    end
    
    XY_path = [X_path, Y_path];
    
%     center_indx = (length(P1_diag)-1)/2+1;
%     
%     P1_diag = P1_diag(center_indx:-1:1);
%     P2_diag = P2_diag(center_indx:-1:1);
%     X_diag = X_diag(center_indx:-1:1);
%     Y_diag = Y_diag(center_indx:-1:1);
end