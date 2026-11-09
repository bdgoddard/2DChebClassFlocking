function HI = ComputeHI2D(rho,Dmu,IntMatrHI)

% note that Dmu is v in the inertial case, but the computation of the term
% needed to give the HI takes exactly the same form.

%     nSpecies=size(rho,2);

    rhoDmu = rho.*Dmu;

    
%     HI_11=zeros(size(rho));
%     HI_12=zeros(size(rho));
    
%     for iSpecies=1:nSpecies
%         for jSpecies=1:nSpecies
%             
%             HI_11(:,iSpecies) =HI_11(:,iSpecies) ...
%               + IntMatrHI(iSpecies,jSpecies).HIInt11*rhoMat;
%             HI_12(:,iSpecies) =HI_12(:,iSpecies) ...
%               + IntMatrHI(iSpecies,jSpecies).HIInt12*rhoDmu(:,jSpecies);
%         end
%     end    
    
    HI_12 = IntMatrHI.HIInt12*rhoDmu;

    M =  IntMatrHI.HIInt11;
    
    M11 = M(1:end/2,1:end/2);
    M21 = M(end/2+1:end,1:end/2);
    M12 = M(1:end/2,end/2 + 1:end);
    M22 = M(end/2 + 1:end,end/2 + 1:end);
    
    rho1 = rho(1:end/2);
    
    M11rho = M11*rho1;
    M12rho = M12*rho1;
    M21rho = M21*rho1;
    M22rho = M22*rho1;
    
    v1 = Dmu(1:end/2);
    v2 = Dmu(end/2+1 : end);
    
    HI_11 = [ (M11rho.*v1 + M12rho.*v2) ; (M21rho.*v1 + M22rho.*v2) ];
    
    HI = HI_11 + HI_12;
    
end
