function I  = GetIndicesQuadrilateral(this)
%************************************************************************
%I  = GetIndicesBox(x1,x2)
%INPUT:
% x1 - grid points in computational space in 1st variable, x1 in [-1 1]
% x2 - grid points in computational space in 2nd variable, x2 in [-1 1]
%OUTPUT: structure with ...
% - indices of 'left', 'right', 'bottom', 'top' boundaries
% - bound = includes left, right, top and bottom boundaries
% - normal vectors: normalLeft,normalRight,normalTop, normalBottom
% - normal = normal vectors for full boundary
% - corners 
% - indices of 'finite', 'infinite' boundaries
% - normal vectors: normalFinite, normalInfinite
%************************************************************************

        x1_kv = this.Pts.x1_kv;
        x2_kv = this.Pts.x2_kv;
        
        y1_kv = this.Pts.y1_kv;
        y2_kv = this.Pts.y2_kv;
        
        % left, right, top, bottom come from ordering of computational
        % domain since corners are sorted into this order in construction
        % of Quadrilateral
        left    = (x1_kv == min(x1_kv));
        right   = (x1_kv == max(x1_kv)); 
        bottom  = (x2_kv == min(x2_kv));
        top     = (x2_kv == max(x2_kv));
                
        bound   = (left | right | bottom | top);
        
        leftFinite   = any(isfinite(y1_kv(left)));
        rightFinite  = any(isfinite(y1_kv(right)));
        topFinite    = any(isfinite(y2_kv(top)));
        bottomFinite = any(isfinite(y2_kv(bottom)));                
        
        finite1 = ( (leftFinite & left) | (rightFinite & right) );
        finite2 = ( (topFinite & top) | (bottomFinite & bottom) );        
        finite  = (finite1 | finite2);
               
        infinite = ( (~leftFinite & left) | (~rightFinite & right) ...
                   | (~topFinite & top) | (~bottomFinite & bottom) );
                              
        corners = ((right & top) | (top & left) | (left & bottom) | (bottom & right));
        
        Z             = zeros(length(x1_kv));

        % normals in y1 direction
        nx1Left       = Z;
        nx1Right      = Z;
        nx1Top        = Z;
        nx1Bottom     = Z;
       
        % normals in y2 direction
        nx2Left       = Z;
        nx2Right      = Z;
        nx2Top        = Z;
        nx2Bottom     = Z;
        
        tl = (top & left);
        tr = (top & right);
        bl = (bottom & left);
        br = (bottom & right);
        
        % use that (y2,-y1) . (y1,y2) = (-y2,y1) . (y1,y2) = 0
        % and pick sign such that the component we know the direction of is
        % correct
        
        n1Temp = y2_kv(tl) - y2_kv(bl); % dy2
        n2Temp = y1_kv(tl) - y1_kv(bl); % dy1
        normN = sqrt(n1Temp^2+n2Temp^2);
        
        n1 = n1Temp/normN;
        n2 = n2Temp/normN;

        nx1Left(left,left)     = -sign(n1)*n1*speye(sum(left));
        nx2Left(left,left)     = sign(n1)*n2*speye(sum(left));
        
        % same for right edge
        n1Temp = y2_kv(tr) - y2_kv(br);
        n2Temp = y1_kv(tr) - y1_kv(br);
        normN = sqrt(n1Temp^2+n2Temp^2);
        
        n1 = n1Temp/normN;
        n2 = n2Temp/normN;

        nx1Right(right,right)  = sign(n1)*n1*speye(sum(right));
        nx2Right(right,right)  = -sign(n1)*n2*speye(sum(right));

        % and for top
        n1Temp = y2_kv(tr) - y2_kv(tl);
        n2Temp = y1_kv(tr) - y1_kv(tl);
        normN = sqrt(n1Temp^2+n2Temp^2);
        
        n1 = n1Temp/normN;
        n2 = n2Temp/normN;

        nx1Top(top,top)  = -sign(n2)*n1*speye(sum(top));
        nx2Top(top,top)  = sign(n2)*n2*speye(sum(top));
        
        % and for bottom
        n1Temp = y2_kv(br) - y2_kv(bl);
        n2Temp = y1_kv(br) - y1_kv(bl);
        normN = sqrt(n1Temp^2+n2Temp^2);
        
        n1 = n1Temp/normN;
        n2 = n2Temp/normN;

        nx1Bottom(bottom,bottom)  = sign(n2)*n1*speye(sum(bottom));
        nx2Bottom(bottom,bottom)  = -sign(n2)*n2*speye(sum(bottom));
        
         
        nf1 = leftFinite*nx1Left(finite,:) + rightFinite*nx1Right(finite,:) + ...
              topFinite*nx1Top(finite,:)   + bottomFinite*nx1Bottom(finite,:);
        nf2 = leftFinite*nx2Left(finite,:) + rightFinite*nx2Right(finite,:) + ...
              topFinite*nx2Top(finite,:)   + bottomFinite*nx2Bottom(finite,:);
                  
        normalFinite = sparse([nf1 nf2] );
                    
        normalFinite1 = sparse( nf1 );
        normalFinite2 = sparse( nf2 );

        normalInfinite = sparse([Z Z]);
        
        nx1 = nx1Left(bound,:) + nx1Right(bound,:) + nx1Top(bound,:) + nx1Bottom(bound,:) ;
        nx2 = nx2Left(bound,:) + nx2Right(bound,:) + nx2Top(bound,:) + nx2Bottom(bound,:) ;
        
        % normalise at corners
        norm = sqrt(nx1.^2 + nx2.^2);
        mask = norm>0;
        nx1(mask) = nx1(mask)./norm(mask);
        nx2(mask) = nx2(mask)./norm(mask);
        
        
        I = struct('left',left,'right',right,'bottom',bottom,'top',top,...
            'bound',bound,...
            'normalLeft',  sparse( [nx1Left(left,:) nx2Left(left,:) ] ),...
            'normalRight', sparse( [nx1Right(right,:) nx2Right(right,:)] ),...
            'normalTop',   sparse( [nx1Top(top,:)  nx2Top(top,:)] ),...
            'normalBottom',sparse( [nx1Bottom(bottom,:)  nx2Bottom(bottom,:)] ),...
            'normalCorners',sparse( [nx1Bottom(corners,:)  nx2Bottom(corners,:)] ),...
            'normal',sparse( [ nx1 nx2 ] ),...
            'corners',corners, ...
            'finite',finite,'infinite',infinite, ...
            'normalFinite',normalFinite,'normalInfinite',normalInfinite,...
            'finite1',finite1,'finite2',finite2,...
            'normalFinite1',normalFinite1,'normalFinite2',normalFinite2);
        
    end