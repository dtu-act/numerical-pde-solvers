function [etovt] = createEToV(numElems1, numElems2)
    etov = zeros(3, 2 * numElems1 * numElems2);
    
    numNodes2 = numElems2 + 1;
    
    for j=0:numElems2 - 1
        for i=0:numElems1 - 1
            index = i * numElems2 + j;

            e2 = index * 2 + 2; % 1-indexed
            e1 = e2 - 1;
            
            % corners in quadrant 'index'
            n1 = (i + 1) * numNodes2 + j + 1; % 1-indexed
            n2 = i * numNodes2 + j     + 1;   % 1-indexed
            n3 = n1 + 1;
            n4 = n2 + 1;
            
            % first local triangle element
            etov(:,e1) = [n1,n2,n3];
            % second local triangle element
            etov(:,e2) = [n4,n3,n2];
        end
    end
    
    etovt = etov';
end