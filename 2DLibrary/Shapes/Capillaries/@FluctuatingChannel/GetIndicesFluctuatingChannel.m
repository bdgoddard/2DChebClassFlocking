function I  = GetIndicesFluctuatingChannel(this)
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

        nx1Left       = Z;
        nx1Right      = Z;
        
        nx1Top       = Z;
        nx1Bottom      = Z;
        nx2Top        = Z;
        nx2Bottom     = Z;
                
        ft      = this.ft;
        fb      = this.fb;
        ftParams = this.ftParams;
        fbParams = this.fbParams;
        
        x1 = this.Pts.x1;
        
        [~, DtopDx1] = ft(x1,ftParams); 
        [~, DbottomDx1] = fb(x1,fbParams);
        
        [~,dy1dx1] =  LinearMap(x1,0,this.L1);
        
        % normals should be in physical space
        DtopDy1 = DtopDx1./dy1dx1;
        DbottomDy1 = DbottomDx1./dy1dx1;
        
        n1TempT = -DtopDy1;
        n2TempT = 1;
        normNT = sqrt(n1TempT.^2+n2TempT.^2);
        n1T = n1TempT./normNT;
        n2T = n2TempT./normNT;
        
        n1TempB = -DbottomDy1;
        n2TempB = 1;
        normNB = sqrt(n1TempB.^2+n2TempB.^2);
        n1B = n1TempB./normNB;
        n2B = n2TempB./normNB;
        
        nx1Left(left,left)     = -speye(sum(left));
        nx1Right(right,right)  = speye(sum(right));
        
        nx1Top(top,top)           = diag(n1T);
        nx1Bottom(bottom,bottom)  = -diag(n1B);
        nx2Top(top,top)           = diag(n2T);
        nx2Bottom(bottom,bottom)  = -diag(n2B);
        
        nf1 = leftFinite*nx1Left(finite,:) + rightFinite*nx1Right(finite,:) + ...
              topFinite*nx1Top(finite,:) + bottomFinite*nx1Bottom(finite,:);
        nf2 = topFinite*nx2Top(finite,:)   + bottomFinite*nx2Bottom(finite,:);
        
        normalFinite = sparse([nf1 nf2] );
                    
        normalFinite1 = sparse( nf1 );
        normalFinite2 = sparse( nf2 );
                    
        normalInfinite = sparse( [((~leftFinite)*nx1Left(infinite,:) + (~rightFinite)*nx1Right(infinite,:)) ...
                        ((~topFinite)*nx2Top(infinite,:) + (~bottomFinite)*nx2Bottom(infinite,:))] );
        
        I = struct('left',left,'right',right,'bottom',bottom,'top',top,...
            'bound',bound,...
            'normalLeft',  sparse( [nx1Left(left,:) Z(left,:) ] ),...
            'normalRight', sparse( [nx1Right(right,:) Z(right,:)] ),...
            'normalTop',   sparse( [nx1Top(top,:)  nx2Top(top,:)] ),...
            'normalBottom',sparse( [nx1Bottom(bottom,:)  nx2Bottom(bottom,:)] ),...
            'normal',sparse( [(nx1Left(bound,:)+nx1Right(bound,:)+nx1Top(bound,:)+nx1Bottom(bound,:)) (nx2Top(bound,:)+nx2Bottom(bound,:))] ),...
            'corners',corners, ...
            'finite',finite,'infinite',infinite, ...
            'normalFinite',normalFinite,'normalInfinite',normalInfinite,...
            'finite1',finite1,'finite2',finite2,...
            'normalFinite1',normalFinite1,'normalFinite2',normalFinite2);
        
    end