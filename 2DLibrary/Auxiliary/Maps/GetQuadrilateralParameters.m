function A = GetQuadrilateralParameters(Y)
%function A = GetQuadrilateralParameters(Y)
%Input: Y(i,:) = (y_1,y_2), coordinates of i-th corner of quadrilateral
% in cyclic order mapped from (-1,-1), (1,-1), (1,1), (-1,1)
%Output: 
% parameters such that
% y = A * (1;x_1;x_2;x_1*x_2)

    B = [1,0,0,0;
         -1,1,0,0;
         -1,0,0,1;
         1,-1,1,-1];
         
    A = (B*Y)';
end