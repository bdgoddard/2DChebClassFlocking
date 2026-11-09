classdef Segment2 < Ball
    %Computes the integral over the surface of a sphere, where it is
    %assumed that the x-y-plane is a axis of symmetry for the function to
    %evaluate. (otherwise inherit from SpectralFourier, not
    %SpectralEvenFourier)
    %See Ball

    methods

        %------------------------------------------------------------------
        % Constructor
        %------------------------------------------------------------------

        function this = Segment2(Geometry)
            this@Ball(Geometry);
        end
    end

    methods (Access = public)

        %------------------------------------------------------------------
        % Overridden initialisations
        %------------------------------------------------------------------


        function [int,area] = ComputeIntegrationVector(this)

            % reweight the integration from ball
            int = ComputeIntegrationVector@Ball(this);
            int = int.*(sin(this.Pts.y1_kv)'.*sin(this.Pts.y2_kv)'/2);

            % exact area
            area = this.R^2*(this.theta2-this.theta1) + ...
                this.R^2/2*(sin(2*this.theta1) - sin(2*this.theta2));

            % display the error in the area if requested
            if(nargout == 2)

                if(area == 0)
                    disp(['Segment2: Area = 0, Error in computation (absolute):',num2str(area-sum(int))]);
                else
                    disp(['Segment2: Error in computation of area (%):',num2str(1-sum(int)/area)]);
                end

            end

        end
    end

end


