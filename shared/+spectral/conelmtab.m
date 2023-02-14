function EToV = conelmtab(noelms1,noelms2)

% total number of elements
M = 2*noelms1*noelms2;
EToV = zeros(M,3);

count = 1;
for j = 1 : noelms1
    % first elements in column
    % upper triangle
    au = count + noelms2+1;
    bu = count;
    cu = au+1;
    % lower triangle
    al = bu+1;
    bl = cu;
    cl  = bu;
    
    for i = 1 : 2: 2*noelms2-1
        % insert rows
        EToV( (j-1)*2*noelms2  +i,:)   = [au,bu,cu];
        EToV( (j-1)*2*noelms2  +i+1,:) = [al,bl,cl];
        
        % translate triangles
        % upper
        au = au+1;
        bu = bu+1;
        cu = cu+1;
        % derive lower
        al = bu+1;
        bl = cu;
        cl  = bu;
    end
    count = count + noelms2+1;
end