classdef Ball < SpectralEvenFourier
    %Computes the integral over the surface of a sphere, where it is
    %assumed that the x-y-plane is a axis of symmetry for the function to
    %evaluate. (otherwise inherit from SpectralFourier, not
    %SpectralEvenFourier)
    %
    %Here, we use spherical coordinates and fix the radial variable to a
    %cosntant. We then consider the projection onto the following plane:
    %
    %here, x = r sin(theta)*cos(phi)
    %      y = r cos(theta)
    %
    %We consider invariance in direction z = r sin(theta)*sin(phi)
    %where phi in [0,pi] and theta is in a given range

    properties
        theta1, theta2 % angle range
        R % radius
        theta0 = 0 % angle offset, used by Shape in GetCartPts
        PtsCart % Cartesian points
    end

    methods

        %------------------------------------------------------------------
        % Constructor
        %------------------------------------------------------------------

        function this = Ball(Geometry)

            this@SpectralEvenFourier(Geometry.N(1),Geometry.N(2));

            % set parameters
            this.theta1 = Geometry.theta1;
            this.theta2 = Geometry.theta2;
            this.R      = Geometry.R;

            if(isfield(Geometry,'Origin'))
                this.Origin = Geometry.Origin;
            end

            if(isfield(Geometry,'theta0'))
                this.theta0 = Geometry.theta0;
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

            % weight the integration for spherical coordinates
            int = ComputeIntegrationVector@SpectralEvenFourier(this);
            int = int.*(2*this.R^2*sin(this.Pts.y1_kv)'); %weight from spherical coordinates

            % update the Int
            this.Int = int;

            % exact area
            area = 2*this.R^2*pi*(cos(this.theta1)-cos(this.theta2));

            % display the error in the area if requested
            if(nargout == 2)
                if(area == 0)
                    disp(['Ball: Area = 0, Error in computation (absolute):',num2str(area-sum(int))]);
                else
                    disp(['Ball: Error in computation of area (%):',num2str(1-sum(int)/area)]);
                end
            end

        end

        % shape only designed for integration
        function ComputeConvolutionMatrix(this) %#ok
            exc = MException('Ball: ComputeConvolutionMatrix','case not implemented');
            throw(exc);
        end
        function Ind    = ComputeIndices(this) %#ok
            exc = MException('Ball: ComputeIndices','case not implemented');
            throw(exc);
        end
        function Diff   = ComputeDifferentiationMatrix(this) %#ok
            exc = MException('Ball: ComputeDifferentiationMatrix','case not implemented');
            throw(exc);
        end

        %------------------------------------------------------------------
        % Maps
        %------------------------------------------------------------------

        function [y1,dy1,dx,ddx,dddx,ddddx] = PhysSpace1(this,x1)
            % x1 is Chebyshev so in [-1,1]
            % linear map from theta1 to theta 2
            % note that theta0 is implemented in Shape when we get
            % Cartesian coordinates
            [y1,dy1,dx,ddx,dddx,ddddx] = LinearMap(x1,this.theta1,this.theta2);
        end

        function xf = CompSpace1(this,y1)
            % inverse linear map
            xf  = InvLinearMap(y1,this.theta1,this.theta2);
        end

        function [th,dth,dx,ddx,dddx,ddddx] = PhysSpace2(~,xT)
            % x2 is Fourier so in [0,1]
            % linear map from 0 to 2\pi
            [th,dth,dx,ddx,dddx,ddddx] = LinearMap01(xT,0,2*pi);
        end

        function xf = CompSpace2(~,th)
            % inverse linear map (on [0,1])
            xf = InvLinearMap01(th,0,2*pi);
        end

    end

end


