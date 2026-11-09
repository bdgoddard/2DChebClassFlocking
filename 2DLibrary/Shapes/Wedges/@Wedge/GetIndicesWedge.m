function I  = GetIndicesWedge(this)
%************************************************************************
%I  = GetIndicesWedge(x1,x2)
%INPUT:
% x1 - grid points in computational space in 1st variable, x1 in [-1 1]
% x2 - grid points in computational space in 2nd variable, x2 in [-1 1]
%OUTPUT: structure with ...
% - indices of 'radial1', 'radial2', 'inner', 'outer' boundaries
% - bound = includes radial_1, radial_2, outer and inner boundaries
% - normal vectors: normalRadial_1, normalRadial_2, normalouter, normalinner
% - normal = normal vectors for full boundary
% - corners
% - indices of 'finite', 'infinite' boundaries
% - normal vectors: normalFinite, normalInfinite
%************************************************************************
%         N1  = this.Pts.N1;
%         N2  = this.Pts.N2;

x1_kv = this.Pts.x1_kv; % r in comp space
x2_kv = this.Pts.x2_kv; % theta in comp space
y1_kv = this.Pts.y1_kv; % r in phys space
y2_kv = this.Pts.y2_kv; % theta in phys space

% make logicals first
outer   = (x1_kv == max(x1_kv)); % out radius = const
inner   = (x1_kv == min(x1_kv)); % in radius = const
radial1 = (x2_kv == min(x2_kv)); % theta1 = const
radial2 = (x2_kv == max(x2_kv)); % theta2 = const

bound   = (radial1 | outer | radial2 | inner);

radial1Finite  = any(isfinite(y2_kv(radial1)));
radial2Finite  = any(isfinite(y2_kv(radial2)));
outerFinite    = any(isfinite(y1_kv(outer)));
innerFinite    = any(isfinite(y2_kv(inner)));

finite1 = ( (radial1Finite & radial1) | (radial2Finite & radial2) );
finite2 = ( (outerFinite & outer) | (innerFinite & inner) );
finite  = (finite1 | finite2);

infinite = ( (~radial1Finite & radial1) | (~radial2Finite & radial2) ...
             | (~outerFinite & outer) | (~innerFinite & inner) );

corners = ((radial1 & outer) | (outer & radial2) | (radial2 & inner) | (inner & radial1));

% named corners
outerradial2 = (outer & radial2);
outerradial1 = (outer & radial1);
innerradial2 = (inner & radial2);
innerradial1 = (inner & radial1);
%------------------------------------------------------------------------
% construct normals for plotting
%------------------------------------------------------------------------

%--------------------------------------------------------------------------
% set up placeholders
%--------------------------------------------------------------------------

Z = zeros(length(x1_kv));

% normals in r (y1) direction
nx1radial2 = Z;
nx1radial1 = Z;
nx1outer   = Z;
nx1inner   = Z;

% normals in theta (y2) direction
nx2radial2 = Z;
nx2radial1 = Z;
nx2outer   = Z;
nx2inner   = Z;

%--------------------------------------------------------------------------
% construct normals for each edge
%--------------------------------------------------------------------------

% Outer normal

% outer normal is in theta direction
n2Temp = y2_kv(outer); % thetas
n1Temp = ones(size(n2Temp)); % rs
nx1outer(outer,outer) = diag(n1Temp);
nx2outer(outer,outer) = diag(n2Temp);

% convert to cartesian
normalOutRPol  = sparse( [nx1outer(outer,:)  nx2outer(outer,:)] );
normalOutRCart = this.GetCartVec(normalOutRPol);

% Inner normal

% inner normal is in theta + pi direction (to be outtward pointing)
n2Temp = y2_kv(inner) + pi; % thetas
n1Temp = ones(size(n2Temp)); % rs
nx1inner(inner,inner) = diag(n1Temp);
nx2inner(inner,inner) = diag(n2Temp);

% convert to cartesian
normalInRPol = sparse( [nx1inner(inner,:)  nx2inner(inner,:)] );
normalInRCart = this.GetCartVec(normalInRPol);

% Radial1 normal

% radial normal is in theta + 3*pi/2 direction (e.g., at theta = 0 the
% normal should point directly downwards)
n2Temp = y2_kv(radial1) + 3*pi/2; % thetas
n1Temp = ones(size(n2Temp)); % rs
nx1radial1(radial1,radial1) = diag(n1Temp);
nx2radial1(radial1,radial1) = diag(n2Temp);

% convert to cartesian
normalRadial_1Pol = sparse( [nx1radial1(radial1,:) nx2radial1(radial1,:)] );
normalRadial_1Cart = this.GetCartVec(normalRadial_1Pol);


% Radial2 normal

% radial normal is in theta + pi/2 direction (e.g., at theta = pi/2 the
% normal should point directly left)
n2Temp = y2_kv(radial2) + pi/2; % thetas
n1Temp = ones(size(n2Temp)); %rs
nx1radial2(radial2,radial2) = diag(n1Temp);
nx2radial2(radial2,radial2) = diag(n2Temp);

% convert to cartesian
normalRadial_2Pol = sparse( [nx1radial2(radial2,:) nx2radial2(radial2,:)] );
normalRadial_2Cart = this.GetCartVec(normalRadial_2Pol);


%--------------------------------------------------------------------------
% construct the full set of theta's and r's
%--------------------------------------------------------------------------

% we can safely add these as they live in different columns, except the
% corners, and we will overwrite those
nThetas = nx2outer(bound,:) + nx2inner(bound,:) + nx2radial2(bound,:) + nx2radial1(bound,:);

% we can or these as they are 0 or 1 and then we don't  need to worry about
% fixing the corners
nRadii  = nx1outer(bound,:) | nx1inner(bound,:) | nx1radial2(bound,:) | nx1radial1(bound,:);

%--------------------------------------------------------------------------
% deal with the corners
%--------------------------------------------------------------------------

psiOuterRadial1 = cornerAngle(nx2radial1(outerradial1,outerradial1), nx2outer(outerradial1,outerradial1));
psiOuterRadial2 = cornerAngle(nx2radial2(outerradial2,outerradial2), nx2outer(outerradial2,outerradial2));
psiInnerRadial2 = cornerAngle(nx2radial2(innerradial2,innerradial2), nx2inner(innerradial2,innerradial2));
psiInnerRadial1 = cornerAngle(nx2radial1(innerradial1,innerradial1), nx2inner(innerradial1,innerradial1));

[~, outerradial1Corner] = ismember(find(outerradial1),find(bound));
[~, outerradial2Corner] = ismember(find(outerradial2),find(bound));
[~, innerradial2Corner] = ismember(find(innerradial2),find(bound));
[~, innerradial1Corner] = ismember(find(innerradial1),find(bound));

nThetas(outerradial1Corner,:) = psiOuterRadial1;
nThetas(outerradial2Corner,:) = psiOuterRadial2;
nThetas(innerradial2Corner,:) = psiInnerRadial2;
nThetas(innerradial1Corner,:) = psiInnerRadial1;

%--------------------------------------------------------------------------
% construct the full normals
%--------------------------------------------------------------------------

normalPol = sparse([nRadii, nThetas]);
normalCart = this.GetCartVec(normalPol);

%------------------------------------------------------------------------
% construct normals for polar computations 
%------------------------------------------------------------------------

nxThetaRadial2 = Z;
nxThetaRadial1 = Z;
nxrOuter = Z;
nxrInner = Z;

% set the theta components for the flat edges - the normal is in the theta
% direction
nxThetaRadial2(radial2,radial2) = speye(nnz(radial2));
nxThetaRadial1(radial1,radial1) = -speye(nnz(radial1));

% construct the radial parts for the curved edges - the normal is in the
% radial direction
nxrOuter(outer,outer) = speye(nnz(outer));
nxrInner(inner,inner) = -speye(nnz(inner));

%------------------------------------------------------------------------
% construct normals finite and infinite
%------------------------------------------------------------------------

nf1 = radial2Finite*nx1radial2(finite,:) + radial1Finite*nx1radial1(finite,:) + ...
      outerFinite*nx1outer(finite,:)   + innerFinite*nx1inner(finite,:);
nf2 = radial2Finite*nx2radial2(finite,:) + radial1Finite*nx2radial1(finite,:) + ...
      outerFinite*nx2outer(finite,:)   + innerFinite*nx2inner(finite,:);

normalFinite = sparse([nf1 nf2]);

normalFinite1 = sparse( nf1 );
normalFinite2 = sparse( nf2 );

normalInfinite = sparse([Z Z]);

%------------------------------------------------------------------------
% Put everything into a structure
%------------------------------------------------------------------------

I = struct('radial_2',radial2,'radial_1',radial1,'inR',inner,'outR',outer,...
        'bound',bound,...
        'normalRadial_2',sparse([Z(radial2,:)  nxThetaRadial2(radial2,:)]),...
        'normalRadial_1',sparse([Z(radial1,:)  nxThetaRadial1(radial1,:)]),...
        'normalOutR',sparse([nxrOuter(outer,:)  Z(outer,:)]),...
        'normalInR',sparse([nxrInner(inner,:)  Z(inner,:)]),...
        'normal',sparse([(nxrOuter(bound,:)+nxrInner(bound,:)) (nxThetaRadial1(bound,:) + nxThetaRadial2(bound,:))]),...
        'normalRadial_2Cart', normalRadial_2Cart,...
        'normalRadial_1Cart', normalRadial_1Cart,...
        'normalOutRCart', normalOutRCart,...
        'normalInRCart', normalInRCart,...
        'normalCart', normalCart,...
        'corners',corners, ...
        'finite',finite,'infinite',infinite, ...
        'normalFinite',normalFinite,'normalInfinite',normalInfinite,...
        'finite1',finite1,'finite2',finite2,...
        'normalFinite1',normalFinite1,'normalFinite2',normalFinite2);


    function psi = cornerAngle(theta1,theta2)
        
        phi = theta2 - theta1;
        if(abs(phi)<=pi)
            psi = theta1 + phi/2; % simple case
        elseif pi<phi && phi<3/2*pi
            psi = theta1 + (2*pi-phi)/2; % flipped case 3rd quadrant
        else
            psi = theta1-(2*pi-phi)/2;  % flipped case 4th quadrant
        end
        psi = mod(psi,2*pi); % probably not needed
    end



end