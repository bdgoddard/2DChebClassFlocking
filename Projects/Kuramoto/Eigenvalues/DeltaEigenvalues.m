% Housekeeping
clear all;
AddPaths;
saveDirBase = ['Kuramoto' filesep 'Stability'];

opts = KuramotoDefaults();

% set up the periodic interval on which to solve the eigenvalue problem
geom.N = 50;  % number of points to use for the pseudospectral method
geom.yMin = 0;    % boundaries of periodic domain
geom.yMax = 2*pi;

geomPlot = geom;   % number of points to use for the plotting
geomPlot.N = 1000;

aLine = FourierLine(geom);
aLine.ComputeAll(geomPlot);


% extract the relevant operators
x = aLine.Pts.y;
Dx = aLine.Diff.Dy;
DDx = aLine.Diff.DDy;
Int = aLine.Int;


% ------------------------ Linear operators ------------------------------

% set the value of n
n=1;

% compute the eigenvalues

delta_values = 0.01:0.01:0.25;
num_deltas = length(delta_values);
first_eigenvalue = zeros(1, num_deltas);

% Iterate over delta values
for i = 1:num_deltas
    delta = delta_values(i);
    eta = delta^2;

    % Define the linear operator

    % Linear operator for W(x) = -cos(nx)
    L = DDx - diag(n^2*delta*cos(n*x) + n^2*delta^2*sin(n*x).^2);

    
    % Calculate eigenvalues
    e = sort(real(eig(L)));
    
    % Store the largest eigenvalue
    first_eigenvalue(i) = e(end-1);
end

% Plot the second eigenvalue against delta
plot(delta_values, first_eigenvalue);
xlabel('$\delta$', 'Interpreter', 'latex');
ylabel('First Eigenvalue');
Wstring = ['$W(x) = -\cos(' num2str(n) 'x)$'];
title(['First Eigenvalue vs. $\delta$, ' Wstring], 'Interpreter', 'latex', 'FontSize', 14);

hold on

% fplot(@(x) -1 + (18/5)*x^2, [0.01, 0.25])
% fplot(@(x) -1 + n^2/2*( n^2/(n^2-4) - 1)*x^2, [0.01, 0.25])
if (n==2)
    analyticExpression = -(n^2)/2*(1/2 + delta_values);
else
    analyticExpression = -1 + (n^2)/2*( (n^2)/(n^2-4) - 1)*delta_values.^2;
end
    
plot(delta_values, analyticExpression);

legend('Numerical eigenvalues', 'Analytical expression','Interpreter', 'latex', 'FontSize', 12)


