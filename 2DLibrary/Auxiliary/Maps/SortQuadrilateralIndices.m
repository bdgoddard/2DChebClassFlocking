function Ysort = SortQuadrilateralIndices(Y)
% Sort indices so that left, right, top and bottom make sense
% when mapped from (-1,-1), (1,-1), (1,1), (-1,1)

    Y1 = Y(:,1);
    Y2 = Y(:,2);

    [~,Ind1] = sort(Y1);
%     [~,Ind2] = sort(Y2);

    % split into 'left' and 'right'
    left = Ind1(1:2);
    right = Ind1(3:4);
%     bottom = Ind2(1:2);
%     top = Ind2(3:4);

    % find top and bottom of left and right pairs
    [~,topLeft] = max(Y2(left));
    [~,topRight] = max(Y2(right));
    bottomLeft = ~(topLeft - 1) + 1;
    bottomRight = ~(topRight - 1) + 1;
%     [~,bottomLeft] = min(Y2(left));
%     [~,bottomRight] = min(Y2(right));

    % find positions in original list
    topLeft = left(topLeft);
    bottomLeft = left(bottomLeft);
    topRight = right(topRight);
    bottomRight = right(bottomRight);
    
    Y1sort = Y1([bottomLeft bottomRight topRight topLeft]);
    Y2sort = Y2([bottomLeft bottomRight topRight topLeft]);
    
    Ysort = [Y1sort Y2sort];
end