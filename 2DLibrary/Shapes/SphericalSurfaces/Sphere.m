classdef Sphere < SpectralSpectral
    %Computes the integral over the surface of a sphere, where it is
    %assumed that the x-y-plane is a axis of symmetry for the function to
    %evaluate. (otherwise inherit from SpectralFourier, not
    %SpectralEvenFourier)
    %
    %here, we use spherical coordinates and fix the radial variable to a
    %cosntant. We then consider the projection onto the following plane:
    %
    %here, x = r sin(theta)*cos(phi)
    %      y = r cos(theta)
    %
    % we consider invariance in direction z = r sin(theta)*sin(phi)
    % where phi in [0,pi] and theta in [0,2pi]

    properties
        theta1 = 0 % angle range
        theta2 = pi
        theta0 = 0 % angle offset
        R % radius
        PtsCart % Cartesian points
        volume = true
    end

    methods

        %------------------------------------------------------------------
        % Constructor
        %------------------------------------------------------------------

        function this = Sphere(Geometry)

            this@SpectralSpectral(Geometry.N(1),Geometry.N(2));

            % set parameters
            if(isfield(Geometry,'theta1'))
                this.theta1 = Geometry.theta1;
            end

            if(isfield(Geometry,'theta2'))
                this.theta2 = Geometry.theta2;
            end

            this.R      = Geometry.R;

            if(isfield(Geometry,'Origin'))
                this.Origin = Geometry.Origin;
            end

            if(isfield(Geometry,'theta0'))
                this.theta0 = Geometry.theta0;
            end

            if(isfield(Geometry,'volume'))
                this.volume = Geometry.volume;
            end

            % set up points using method from Shape (and PhysSpace below)
            InitializationPts(this);

            this.polar = 'sphSurf';
        end
    end

    methods (Access = public)

        %------------------------------------------------------------------
        % Overridden initialisations
        %------------------------------------------------------------------

        function [int,area] = ComputeIntegrationVector(this)

            R   = this.R; %#ok<*PROP>
            y1  = this.Pts.y1_kv;
            y2  = this.Pts.y2_kv;
            th1 = this.theta1;
            th2 = this.theta2;

            int = ComputeIntegrationVector@SpectralSpectral(this);
            int = int.*(R^2*sin(y2)'.*sin(y1').^2); %weight from spherical coordinates

            % adjust integration weights if we're using a volume
            % and compute exact area
            if(this.volume)
                int  = int.*(2*R*sin(y1).*sin(y2))';
                area = pi*R^3/6*( (cos(th2)*(cos(2*th2)-5)) - ...
                    (cos(th1)*(cos(2*th1)-5)));
            else
                area = R^2/2*(2*(th2-th1)-(sin(2*th2)-sin(2*th1)));
            end

            % update the Int
            this.Int = int;

            % display the error in the area if requested
            if(nargout == 2)
                if(area == 0)
                    disp(['Ball: Area = 0, Error in computation (absolute):',num2str(area-sum(int))]);
                else
                    disp(['Ball: Error in computation of area (ratio):',num2str(1-sum(int)/area)]);
                end
            end

        end

        % shape only designed for integration
        function ComputeConvolutionMatrix(this) %#ok
            exc = MException('Sphere: ComputeConvolutionMatrix','case not implemented');
            throw(exc);
        end

        function Ind    = ComputeIndices(this) %#ok
            exc = MException('Sphere: ComputeIndices','case not implemented');
            throw(exc);
        end
        function Diff   = ComputeDifferentiationMatrix(this) %#ok
            exc = MException('Sphere: ComputeDifferentiationMatrix','case not implemented');
            throw(exc);
        end

        %------------------------------------------------------------------
        % Maps
        %------------------------------------------------------------------

        function [y,dy,dx,ddx,dddx,ddddx] = PhysSpace1(this,x)
            % linear map from theta1 to theta 2
            % note that theta0 is implemented in Shape when we get
            % Cartesian coordinates
            [y,dy,dx,ddx,dddx,ddddx] = LinearMap(x,this.theta1,this.theta2);
        end

        function xf = CompSpace1(this,y)
            % inverse linear map
            xf  = InvLinearMap(y,this.theta1,this.theta2);
        end

        function [y,dy,dx,ddx,dddx,ddddx] = PhysSpace2(~,x)
            % linear map
            [y,dy,dx,ddx,dddx,ddddx] = LinearMap(x,0,pi);
        end

        function xf = CompSpace2(~,y)
            % inverse linear map
            xf = InvLinearMap(y,0,pi);
        end
        
    end

end


