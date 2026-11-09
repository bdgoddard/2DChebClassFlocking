clear all
AddPaths;

% set up the Shape
shapeName = 'Quadrilateral';
geom = eval([shapeName '.DocumentationPlotDefaultGeometry']);
geom.N = 4*geom.N; % default is normally quite small for plotting clarity
shapeFn = str2func(shapeName);
if(nargin(shapeFn) == 0) % specified by argument list
    aShape = shapeFn(geom{:});
else % specified by structure
    aShape = shapeFn(geom);
end

% compute operators
[Pts,Diff,Int,Ind] = aShape.ComputeAll;

% compute interpolation for plotting and checking accuracy
Interp = aShape.ComputeInterpolationMatrix((-1:0.02:1)',(-1:0.02:1)',true,true);

% shorthands
y1 = Pts.y1_kv;
y2 = Pts.y2_kv;
Origin = aShape.Origin;

% determine exact reference
[fExact,DiffExact] = ExactDifferentiation(y1,y2,Origin(1),Origin(2));

% shorthands
Dy1 = Diff.Dy1;
Dy2 = Diff.Dy2;
DDy1 = Diff.DDy1;
DDy2 = Diff.DDy2;

% approximate derivatives in y1, y2
D1Approx = Dy1 * fExact;
D2Approx = Dy2 * fExact;
DD1Approx = DDy1 * fExact;
DD2Approx = DDy2 * fExact;

disp(['Dy1 pointwise error = ' num2str(max( abs(DiffExact.D1 - D1Approx) ))]);
disp(['Dy2 pointwise error = ' num2str(max( abs(DiffExact.D2 - D2Approx) ))]);
disp(['DDy1 pointwise error = ' num2str(max( abs(DiffExact.DD1 - DD1Approx) ))]);
disp(['DDy2 pointwise error = ' num2str(max( abs(DiffExact.DD2 - DD2Approx) ))]);

% function [f,D1f,D2f] = ExactDifferentiation(y1,y2,y10,y20)
%     f = ones(size(y1));
%     D1f = zeros(size(f));
%     D2f = zeros(size(f));
% end



% function [f,D1f,D2f] = ExactDifferentiation(y1,y2,y10,y20)
%     f = y1;
%     D1f = ones(size(f));
%     D2f = zeros(size(f));
% end


function [f,Diff] = ExactDifferentiation(y1,y2,y10,y20)
    f = exp( - ( (y1-y10).^2 + (y2-y20).^2 ) );
    Diff.D1 = - 2 * (y1-y10) .* f;
    Diff.D2 = - 2 * (y2-y20) .* f;
    Diff.DD1 = - 2 .* f - 2 * (y1-y10) .* Diff.D1;
    Diff.DD2 = - 2 .* f - 2 * (y2-y20) .* Diff.D2;
end

    function y = ExactSolution(y1,y2,t)
        y10 = 5;
        y20  = 3;
        y = 1/(4*pi*t) * exp(-((y1-y10).^2 + (y2-y20).^2)/4/t);
       
    end
