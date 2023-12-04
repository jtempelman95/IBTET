classdef excitation

    properties
        omega
        periods
        amp
        FLOC
        tau
        t0
        steadystate
        band
        BandLoc
    end
    
    methods 

        function self = excitation(omega,periods,amp,FLOC,tau,t0,steadystate)
            self.omega = omega;
            self.periods = periods;
            self.amp = amp;
            self.FLOC = FLOC;

            if nargin < 7
                self.steadystate = 0;
            else
                self.steadystate = steadystate;
            end
            if nargin < 6
                self.t0 = [];
            else
                self.t0 = t0;
            end
            if nargin < 5
                self.tau = [];
            else
                self.tau = tau;
            end
        end

        function Fvec = Ffunction(self,tf)
            fin = self.omega/2/pi;
            Ncy = self.periods;
            amp = self.amp;
            
            Fvecf = @(t,Ncy,fin,amp) amp*0.5*(heaviside(t)-...
            heaviside(t-2*pi/(2*pi*fin)*Ncy)).*(1-cos(2*pi*fin/Ncy*t)).*sin(2*pi*fin*t);
        
            Fvec =Fvecf(tf,Ncy,fin,amp);
        end

        function TF = tforce(self)
            periods = self.periods;
            freq = self.omega;
            TF =2*pi*periods/freq;
        end
    end

end