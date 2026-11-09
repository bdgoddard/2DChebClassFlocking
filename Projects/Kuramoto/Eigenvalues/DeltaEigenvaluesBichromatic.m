% Housekeeping
clear all;
AddPaths;
saveDirBase = ['Kuramoto' filesep 'Stability'];

opts = KuramotoDefaults();

% set up the periodic interval on which to solve the eigenvalue problem
geom.N = 100;  % number of points to use for the pseudospectral method
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

% compute the eigenvalues

delta_values = 0.001:0.01:0.25;
eta_values = sqrt(delta_values);

num_deltas = length(delta_values);
first_eigenvalue = zeros(1, num_deltas);

evalsList = zeros(3,num_deltas);

% Iterate over delta values
for i = 1:num_deltas
    delta = delta_values(i);

    % Define the linear operator
    %L = DDx - diag(n^2*delta*cos(n*x) + n^2*delta^2*sin(n*x).^2);
    L = DDx - diag(sqrt(3/2)*delta*sin(x).^2 + (3/2)^(1/4)*sqrt(delta)*cos(x) ...
                   + 2*sqrt(3/2)*delta*cos(2*x));

    % Calculate eigenvalues
    e = sort(real(eig(L)));
    
    % Store the largest eigenvalue
    first_eigenvalue(i) = e(end-2);
    evalsList(:,i) = e(end-2:end);
end

figure
% Plot the second eigenvalue against delta
%plot(eta_values, first_eigenvalue);
%hold on
plot(eta_values, evalsList);
xlabel('$\eta$', 'Interpreter', 'latex');
ylabel('Eigenvalues');
Wstring = ['$W(x) = -\cos(x) - 1/2 \cos(2x)$'];
title(['First Eigenvalue vs. $\eta$, ' Wstring], 'Interpreter', 'latex', 'FontSize', 14);

hold on

analyticExpression = -1 - delta_values*( (36^(1/4))/12 + 3/2*sqrt(3/2) );

plot(eta_values, analyticExpression);

legend('e(end-2)', 'e(end-1)', 'e(end)','Analytical expression','Interpreter', 'latex', 'FontSize', 12)


