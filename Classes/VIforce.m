classdef VIforce

    properties
        TimeVI
        rescoeff
        E
        nu 
        kc
        radius
        sided
        dnominal
        VIlocs
        VIsites;
        config;
        dVI    
        d
    end
    
    methods 
        function self = VIforce(TimeVI, E, nu, radius ,rescoef , dnominal,...
                sided, VIlocs, dVI)
        self.TimeVI     = TimeVI;
        self.E          = E;
        self.nu         = nu;
        self.rescoeff   = rescoef;
        self.kc         = 2*E*sqrt(radius)/(3*(1-nu^2)); 
        self.dnominal   = dnominal;
        self.d          = dnominal;
        self.VIlocs     = VIlocs;
        self.VIsites    = VIlocs;
        self.dVI        = dVI;
        self.d(VIlocs)  = dVI;
        self.radius     = radius;
        self.sided      = sided;
        end


%         function self = VIforceVarKC(TimeVI, E, nu, radius ,rescoef , dnominal,...
%         sided, VIlocs, dVI)
%         self.TimeVI     = TimeVI;
%         self.E          = E;
%         self.nu         = nu;
%         self.rescoeff   = rescoef;
%         self.kc         = 2*E*sqrt(radius)/(3*(1-nu)); 
%         self.dnominal   = dnominal;
%         self.d          = dnominal;
%         self.VIlocs     = VIlocs;
%         self.VIsites    = VIlocs;
%         self.dVI        = dVI;
%         self.d(VIlocs)  = dVI;
%         self.radius     = radius;
%         self.sided      = sided;
%         end

    end
end