classdef VI_LatticeSys  < matlab.mixin.Copyable 
    % Class defintion to setup matrices and simulation functions for VI
    % lattice system
    % ------------------------------------------
    % Recover stiffness and mass matrices 
    % ------------------------------------------
    %   k1      - intercell force
    %   k2      - intracell force
    %   m1(m2)  - mass 1(2)
    %   Kg      - grounding stifness
    %   grounded- if grounded or not
    %   config  - lattice configuration
    %   ndof    - size of lattice (degrees of freedom)
    %   VIsties - locations of vibro impacts (if different nominal
    %   stiffness desired)
    %   KVI     - stiffness at VI sites.

    properties 
        config
        ndof
        m
        k
        c % c(1) - stiffenss prop, c(2) - mass prop
        kg
        kNL_g
        grounded
        VIsites
        KVI
        K
        M
        VIforce
        excitation
        simdata
        statematrix
        Minv
       % wavelet_params
        Dispersion
        SimParams
    end

   % Data set for vibro impact lattice 
   methods 

%%
    % ------------------------------------------
    function self = VI_LatticeSys(ndof,Config,m,k,VIsites,KVI,kg,grounded)
    % ------------------------------------------
    % Initialize Class Object
    % ------------------------------------------
            if nargin <  8
                self.grounded = 0;
            else
                self.grounded = grounded;
            end
            if nargin < 7
                self.kg = k2;
            end
            if nargin < 6
                KVI = k2;
            else
                self.kg       = kg;
            end
            if nargin < 5
                VIloc = 0;
            end
            self.config = Config;
            self.ndof   = ndof;
            self.m      = m;
            self.k      = k;
            self.VIsites = VIsites;
            self.KVI    = KVI;
    end
%%
    % ------------------------------------------
    function [K,M,C] = system_mats(self,Ccoeff,Cnorm)
    % ------------------------------------------
    % Function: Matrix Assembly;
    % ------------------------------------------
    if nargin < 3;
        Cnorm = .1;
    end
    if nargin < 2;
        Ccoeff = 0;
    end
        k1 = self.k(1);
        k2 = self.k(2);
        if self.config ==3; k3 = self.k(3);end
        if self.config ==4; k3 = self.k(3);end
        if self.config ==4; k4 = self.k(4);end
        m1 = self.m(1);
        m2 = self.m(2);
        if self.config == 3; m3 = self.m(3);end
        if self.config == 4; m3 = self.m(3); m4 =self.m(4);end
        Kg = self.kg;
        ndof = self.ndof;
        K = zeros(ndof);
        m = zeros(ndof,1);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
            if self.config == 1
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
    %       CONFIGURATION 1          %
    % --------------------------------
        for i = 1:2:ndof
             K(i,i)     = k1+k2;
             K(i,i+1)   = -k2;
             if i > 2
             K(i,i-1)   = -k1;
             end
             m(i)       = m1;
        end
        for i = 2:2:ndof
            K(i,i)      = k1+k2;
            if i < ndof
             K(i,i+1)   = -k1;
            end
             K(i,i-1)   = -k2;
             m(i)       = m2;
        end
        % -------------------------------
        % Define BCs
        % -------------------------------
        %koff(ndof-1)    = -k1;              
        K(1,1)          = k1+k1;  %Fixed BC           
        K(end)          = k1+k2;  %Fixed BC     
        m(ndof-1)       = m1;                
        m(ndof)         = m2;
        % -------------------------------
            elseif self.config == 2
        % -------------------------------
        % CONFIGURATION 2                *
        % -------------------------------
        for i = 1:2:ndof-3  % (a sites)
            K(i,i)      = k1+k2;
            koff(i)     = -k1;
            K(i,i+1)    = -k2;
            K(i,i+2)    = -k1;       
            if i >2
            K(i,i)      = k2 + 2*k1;
            K(i,i-2)    = -k1;
            end
            m(i)        = m1;
        end
        for i = 2:2:ndof-2  %n(b sites)
            K(i,i)      = k2;
            K(i,i-1)    = -k2;
            m(i)        = m2;
        end
        % -------------------------------       
        % Define BCs                            
        % -------------------------------       
        K(1)            =   2*k1+k2     ; % Fixed BC
        K(ndof,ndof)    =   k2          ; 
        K(ndof-1,ndof-1)= 2*k1 + k2    ;   % Fixed BC
        K(ndof-1,ndof-3)=  -k1;
        K(end,ndof-1 )  =   -k2;
        K(end)          =   k2;
        K(end-1,end)    = -k2 ;
        m(ndof-1)       =   m1;
        m(ndof)         =   m2;
        
        
        % Assign springs of the VI sites
        for RR = [self.VIsites]+1
            K(RR,RR)        =  self.KVI          ;
            K(RR,RR-1)      = -self.KVI         ;
            K(RR-1,RR-1)    =  2*k1+self.KVI     ;
            K(RR-1,RR)      = -self.KVI          ;
        end
        if self.KVI~=self.k(2);disp('warning: Nominal VI stiffness different than rest of lattice');end


        % -------------------------------
            elseif self.config == 3
        % -------------------------------
        % CONFIGURATION 3                *
        % -------------------------------
            for i = 1:3:ndof-3 % (a sites)
                K(i,i)      = k1+k2;
                koff(i)     = -k1;
                K(i,i+1)    = -k2;
                K(i,i+3)    = -k1;
                if i >2
                    K(i,i)      = k2 + 2*k1;
                    K(i,i-3)    = -k1;
                end
                m(i)        = m1;
            end
            for i = 2:3:ndof-3  %(b sites)
                    K(i,i)      = k3+k2;
                    K(i,i-1)    = -k2;
                    K(i,i+1)    = -k3;
                    m(i)        = m2;
            end
              for i = 3:3:ndof-3  %(c sites)
                    K(i,i)      = k3;
                    K(i,i-1)    = -k3;
                    m(i)        = m3;
              end

        m(ndof-2)       =   m3;
        m(ndof-1)       =   m1;
        m(ndof)         =   m2;
        % -------------------------------       
        % Define BCs                            
        % -------------------------------       
        % fixed bc
        K(1,1) = 2*k1+k2;
        % Third to last row
        K(end-2,end-2) = 2*k1+k2;
        K(end-2,end-5) = -k1;
        K(end-2,end-1) = -k2;
        
        % secind to last row
        K(end-1,end-1) = k2+k3;
        K(end-1,end-2) = -k2;
        K(end-1,end)   = -k3;
        % last row
        K(end,end) = k3;
        K(end,end-1) = -k3;

        % Assign springs of the VI sites
        for RR = [self.VIsites]+1
            K(RR,RR)        =  self.KVI          ;
            K(RR,RR-1)      = -self.KVI         ;
            K(RR-1,RR-1)    =  2*k1+self.KVI     ;
            K(RR-1,RR)      = -self.KVI          ;
        end

        % -------------------------------
            elseif self.config == 4
         % -------------------------------
        % CONFIGURATION 4               *
        % -------------------------------
            for i = 1:4:ndof-4 % (a sites)
                K(i,i)          = k1+k2;
                koff(i)         = -k1;
                K(i,i+1)        = -k2;
                if i >2
                    K(i,i)      =  k2 + k1;
                    K(i,i-1)    = -k1;
                end
                m(i)        = m1;
            end
            for i = 2:4:ndof-4  %(b sites)
                    K(i,i)      = k3+2*k2;
                    K(i,i-1)    = -k2;
                    K(i,i+1)    = -k3;
                    K(i,i+2)    = -k2;
                    m(i)        =  m2;
            end
              for i = 3:4:ndof-4  %(c sites)
                    K(i,i)      = k3;
                    K(i,i-1)    = -k3;
                    m(i)        = m3;
              end
              for i = 4:4:ndof-4  %(c sites)
                    K(i,i)      = k1+k2;
                    K(i,i-2)    = -k2;
                    K(i,i+1)    = -k1;
                    m(i)        = m4;
              end

              
        m(ndof-3)       =   m1;
        m(ndof-2)       =   m3;
        m(ndof-1)       =   m1;
        m(ndof)         =   m2;
        % -------------------------------       
        % Define BCs                            
        % -------------------------------       
        % fixed bc
        K(1,1) = k1+k2;


        % Fourth to last row
        K(end-3,end-3) = k1+k2;
        K(end-3,end-2) = -k2;
        K(end-3,end-4) = -k1;

        % Third to last row
        K(end-2,end-2) = 2*k2+k3;
        K(end-2,end-3) = -k2;
       K(end-2,end-1) = -k3;
        K(end-2,end)   = -k2;
        
        % secind to last row
        K(end-1,end-1) =  k3;
        K(end-1,end-2) = -k3;

        % last row
        K(end,end)  =   k1 + k2;
        K(end,end-2) = -k2;

        % Assign springs of the VI sites
       if self.KVI~=self.k(3)
                    
           disp('warning: Nominal VI stiffness different than rest of lattice');
        for RR = [self.VIsites]+1
            K(RR,RR)        =  self.KVI          ;
            K(RR,RR-1)      = -self.KVI         ;
            K(RR-1,RR-1)    =  k2+k4+self.KVI     ;
            K(RR-1,RR)      = -self.KVI          ;
        end
        end
        % --------------------------------
        end
        % --------------------------------
        
        if self.grounded ==1
            if self.config<3 
            K = K + eye(size(K))*Kg;
            elseif self.config == 4
                gvec = zeros(self.ndof,1);
                asites = 1:4:self.ndof;
                gvec(asites) = Kg;
                K = K + diag(gvec);

            end
        end
        % -------------------------------       
        % MASS MATRIX ASSEMBLY
        % -------------------------------       
        M = diag(m);


        % -------------------------------       
        % Damping Matrix
        % -------------------------------    
        if nargin == 2
        C = Ccoeff(1)*K + Ccoeff(2)*M;
        else
            C = [];
        end
        if nargin == 3
            C = C/norm(C)*Cnorm;
        end
        if ~isempty(self.c)
            C = K*self.c(1) + M*self.c(2);
        end
      end

%%
   % ------------------------------------------
   function [As,Minv] = simulation_mats(self)
   % ------------------------------------------
   % Get the system matrices for simulation
   % ------------------------------------------
   if sum(self.K(:)) == 0
   [K,M,C] =            self.system_mats();
   else
            K = self.K;
            M = self.M;
            C = self.C;
   end
   Minv = inv(M);
       if isempty(C)
        As = [zeros(self.ndof) eye(self.ndof);-Minv*K zeros(size(M))] ;
       else
        As = [zeros(self.ndof) eye(self.ndof);-Minv*K -Minv*C]  ;
       end
   end
%%
   % ------------------------------------------
   function vg = group_vel(self,omega,band,bandloc)
   % ------------------------------------------
   % Get the group velocity of a traveling wave
   % ------------------------------------------
  
   
   % Get the wave vector
   % - this requires the dispersion function to have been run
   kappa  = self.Dispersion.wavevec;
   
   % Find what band the frequency is on
   DISP = self.Dispersion.DISP;

   for k = 1:size(DISP,2)
       if omega >= min(DISP(:,k)) & omega <= max(DISP(:,k))
           band = k;
       end
   end

   Band = DISP(:,band);

   if nargin>3
        idx = find(kappa > bandloc);
        idx = idx(1);
   else
        [val,idx]=min(abs(omega-Band));
       
   end
   
   dwdk = diff(Band)/mean(diff(kappa));

   vg = abs(dwdk(idx));
   
   end


   %%
   % ------------------------------------------
   function k = K_max_group_vel(self,band)
   % ------------------------------------------
   % Get the group velocity of a traveling wave
   % ------------------------------------------
  
   
   % Get the wave vector
   % - this requires the dispersion function to have been run
   kappa  = self.Dispersion.wavevec;
   kappaidx  = find(0<kappa & kappa<pi);
   kappa = kappa(kappaidx);

   % Find what band the frequency is on
   DISP = self.Dispersion.DISP;

   Band = DISP(kappaidx,band);
   
   dwdk = abs(diff(Band)/mean(diff(kappa)));

   idx = find(dwdk == max(dwdk));

   k = kappa(idx);
   
   end
   

   %%
   % ------------------------------------------
   function [Phi,omega] = modal(self)
   % ------------------------------------------
   % Get the mass orthonormalized modal matrix
   % ------------------------------------------
   % K - stiffness matrix of mechanical sysmet
   % M - mass matrix of mechanical system\
   [K,M] = self.system_mats();
   [Phi,Eval] = eig(K,M);
    
   for  k = 1:length(Phi)
        Phi(:,k) = Phi(:,k)/sqrt(Phi(:,k)'*M*Phi(:,k));
   end
        omega = sqrt(diag(Eval));
   end

     % ------------------------------------------
     function omega = freqs(self,norm)
   % ------------------------------------------
   % Get the mass orthonormalized modal matrix
   % ------------------------------------------
   % K - stiffness matrix of mechanical sysmet
   % M - mass matrix of mechanical system\
   [K,M]     = self.system_mats();
   [Phi,Eval] = eig(K,M);
   end


%%
   % ------------------------------------------
   function [DISP,wavevec] = dispersion(self,dkappa)
        % ------------------------------------------
        % Get the infinite dispersion of the system
        % ------------------------------------------
        % K - stiffness matrix of mechanical sysmet
        % M - mass matrix of mechanical system
        if nargin < 2
            dkappa = pi/100;
        end

        m1 = self.m(1); m2 = self.m(2);
        k1 = self.k(1); k2 = self.k(2);
        if self.config >2; k3 =self.k(3); m3 = self.m(3);end
        if self.config >3; k4 =self.k(4); m3 = self.m(3);end
        grounded = self.grounded;
        Kg = self.kg;
        R = 1;
        wavevec = -2*pi:dkappa:2*pi;
        Acoustic = zeros(size(wavevec));
        Optical  = zeros(size(wavevec));
        % --------------------------------------------
        % Loop through and compute eigenspectrum
        % --------------------------------------------
        if self.config <3
        M = diag([m1, m2]);
        elseif self.config == 3
        M = diag([m1,m2,self.m(3)]);
        elseif self.config == 4
        M = diag([m1,m2,self.m(3),self.m(4)]);
        end
        if isempty(self.grounded);grounded =0;end
        if isempty(self.kg);Kg =0; end;
        for mu = wavevec;
            if self.config ==1
                H = [Kg*grounded+k1+ k2 -k1-k2*exp(-sqrt(-1)*mu);
                -k1-k2*exp(sqrt(-1)*mu) Kg*grounded+k1+k2];
            elseif self.config == 2
                H = [Kg*grounded+2*k1+k2-k1*exp(-sqrt(-1)*mu)-k1*exp(sqrt(-1)*mu)     -k2;
                -k2                  Kg*grounded+k2];
            elseif self.config ==3
                 H = [Kg*grounded+2*k1+k2-k1*exp(-sqrt(-1)*mu)-k1*exp(sqrt(-1)*mu)     -k2 0 ;
                -k2                  k2+k3     -k3           ;
                  0                        -k3              k3                    ];
            elseif self.config ==4
                     H  = [   Kg*grounded+k1+k2 -k2 0 -k1*exp(-sqrt(-1)*mu);
                        -k2    k2+k2+k3  -k3  -k2;
                        0       -k3       k3    0;
                        -k1*exp(sqrt(-1)*mu)   -k2 0 k1+k2];
            end
            [~,eval] = eig(H,M);
            ev = sqrt(sort(diag(eval)));
            Acoustic(R) = ev(1);
            Optical(R)  = ev(2);
            DISP(R,:) = sort(abs(ev));
            kappa(R) = mu;
            R = R + 1;
        end
   end
   
%%
       % ------------------------------------------
       function simsave(self, tsim, IC, opts)
       % ------------------------------------------
       % Simulate the system
       % ------------------------------------------
       % K - stiffness matrix of mechanical sysmet
       % M - mass matrix of mechanical system
       % C - Damping matrix of mechanical system
       % opts - options of ode integrator
       
       if nargin < 4
           opts = self.SimParams.opts;
       end 
       if nargin < 3   
            IC = self.SimParams.IC;
       end
       if nargin < 2
           tsim = self.SimParams.tsim;
       end
       if ~isempty(self.excitation), Fin = self.excitation; else, disp('Need to define forcing profle'),end
       if ~isempty(self.VIforce),    FVI = self.VIforce;   else, disp('Need to define nonlinear contacts'),end
       
       if isempty(self.Minv); [As,Minv] = self.simulation_mats();
       else As = self.statematrix; Minv = self.Minv;
       end
        
       if self.excitation.steadystate ~= 1
       opts.Events = @myEvent;
       end
    
       % Getting the 'primer data' for the right boundary condition
       [tset, yset] = ode45(@(t,y) sys2_global(t,y,Minv,As,Fin,FVI,self),tsim(1:1000),IC,odeset('RelTol',opts.RelTol))     ;
       global maxy;
       maxy = max(yset(:));

       [tsim,ysim]  = ode78(@(t,y) sys2_global(t,y,Minv,As,Fin,FVI,self),tsim,IC,opts)     ;
        
       
           self.simdata.time         = tsim                              ;
           self.simdata.velocity     = ysim(:,(1+self.ndof):self.ndof*2) ;
           self.simdata.displacement = ysim(:,1:self.ndof)               ;
           
       end


%%
      % ------------------------------------------
       function [tsim,ysim] = simulate(self, tsim, IC, opts)
       % ------------------------------------------
       % Call a simulation code
       % ------------------------------------------
       if nargin < 4
           opts = self.SimParams.opts;
       end 
       if nargin < 3   
            IC = self.SimParams.IC;
       end
       if nargin < 2
            tsim = self.SimParams.tsim;
       end
            [tsim,ysim] = simulate_pers(self, tsim, IC, opts);
       end

       % ------------------------------------------
       function [tsim,ysim] = simulate_global(self, tsim, IC, opts)
       % ------------------------------------------
       % Simulate the system (global variable code)
       % ------------------------------------------
       % K - stiffness matrix of mechanical sysmet
       % M - mass matrix of mechanical system
       % C - Damping matrix of mechanical system
       % opts - options of ode integrator
         if nargin < 4
           opts = self.SimParams.opts;
       end 
       if nargin < 3   
            IC = self.SimParams.IC;
       end
       if nargin < 2
           tsim = self.SimParams.tsim;
       end
       if ~isempty(self.excitation), Fin = self.excitation; else, disp('Need to define forcing profle'),end
       if ~isempty(self.VIforce),    FVI = self.VIforce;   else, disp('Need to define nonlinear contacts'),end
       
       if isempty(self.Minv); [As,Minv] = self.simulation_mats();
       else As = self.statematrix; Minv = self.Minv;
       end
       
       if self.excitation.steadystate ~= 1
       opts.Events = @myEvent;
       end
    
       global Index UDOTMINUS %#ok<*NUSED> 
       clear Index UDOTMINUS
       % Getting the 'primer data' for the right boundary condition
       [~, yset] = ode45(@(t,y) sys2_global(t,y,Minv,As,Fin,FVI,self),tsim(1:500),IC,odeset('RelTol',opts.RelTol))     ;
       global maxy;
       maxy = max(yset(:));
       
       n = numel(self.SimParams.Integrator);
       if n == 5
            if self.SimParams.Integrator == 'ode78'
            [tsim,ysim]  = ode78(@(t,y) sys2_global(t,y,Minv,As,Fin,FVI,self),tsim,IC,opts)     ;
            elseif self.SimParams.Integrator == 'ode89'
            [tsim,ysim]  = ode89(@(t,y) sys2_global(t,y,Minv,As,Fin,FVI,self),tsim,IC,opts)     ;
            elseif self.SimParams.Integrator == 'ode45'   
            [tsim,ysim]  = ode78(@(t,y) sys2_global(t,y,Minv,As,Fin,FVI,self),tsim,IC,opts)     ;
            end
       elseif n == 6
            if self.SimParams.Integrator == 'ode113'
            [tsim,ysim]  = ode113(@(t,y) sys2_global(t,y,Minv,As,Fin,FVI,self),tsim,IC,opts)     ;
            elseif self.SimParams.Integrator == 'ode15s'
            [tsim,ysim]  = ode15s(@(t,y) sys2_global(t,y,Minv,As,Fin,FVI,self),tsim,IC,opts)     ;
            end
       end
       end

       % ------------------------------------------
       function [tsim,ysim] = simulate_pers(self, tsim, IC, opts)
       % ------------------------------------------
       % Simulate the system  (Persistitnnt Variable code) || DEFUALT
       % ------------------------------------------
       % K - stiffness matrix of mechanical sysmet
       % M - mass matrix of mechanical system
       % C - Damping matrix of mechanical system
       % opts - options of ode integrator
         if nargin < 4
           opts = self.SimParams.opts;
       end 
       if nargin < 3   
            IC = self.SimParams.IC;
       end
       if nargin < 2
           tsim = self.SimParams.tsim;
       end
       if ~isempty(self.excitation), Fin = self.excitation; else, disp('Need to define forcing profle'),end
       if ~isempty(self.VIforce),    FVI = self.VIforce;   else, disp('Need to define nonlinear contacts'),end
       
       if isempty(self.Minv); [As,Minv] = self.simulation_mats();
       else As = self.statematrix; Minv = self.Minv;
       end
        
       if self.excitation.steadystate ~= 1
       opts.Events = @myEvent;
       end
    

       % Getting the 'primer data' for the right boundary condition
       [~, yset] = ode45(@(t,y) sys2_persistent(t,y,Minv,As,Fin,FVI,self),tsim(1:500),IC,odeset('RelTol',opts.RelTol))     ;
       global maxy;
       maxy = max(yset(:));
       
       n = numel(self.SimParams.Integrator);
       if n == 5
            if self.SimParams.Integrator == 'ode78'
            [tsim,ysim]  = ode78(@(t,y) sys2_persistent(t,y,Minv,As,Fin,FVI,self),tsim,IC,opts)     ;
            elseif self.SimParams.Integrator == 'ode89'
            [tsim,ysim]  = ode89(@(t,y) sys2_persistent(t,y,Minv,As,Fin,FVI,self),tsim,IC,opts)     ;
            elseif self.SimParams.Integrator == 'ode45'   
            [tsim,ysim]  = ode45(@(t,y) sys2_persistent(t,y,Minv,As,Fin,FVI,self),tsim,IC,opts)     ;
            end
       elseif n == 6
            if self.SimParams.Integrator == 'ode113'
            [tsim,ysim]  = ode113(@(t,y) sys2_persistent(t,y,Minv,As,Fin,FVI,self),tsim,IC,opts)     ;
            elseif self.SimParams.Integrator == 'ode15s'
            [tsim,ysim]  = ode15s(@(t,y) sys2_persistent(t,y,Minv,As,Fin,FVI,self),tsim,IC,opts)     ;
            end
       end
       end



       % ------------------------------------------
       function [tsim,ysim] = simulate_bilinear(self, tsim, IC, opts)
       % ------------------------------------------
       % Simulate the system  (Persistitnnt Variable code) || DEFUALT
       % ------------------------------------------
       % K - stiffness matrix of mechanical sysmet
       % M - mass matrix of mechanical system
       % C - Damping matrix of mechanical system
       % opts - options of ode integrator
         if nargin < 4
           opts = self.SimParams.opts;
       end 
       if nargin < 3   
            IC = self.SimParams.IC;
       end
       if nargin < 2
           tsim = self.SimParams.tsim;
       end
       if ~isempty(self.excitation), Fin = self.excitation; else, disp('Need to define forcing profle'),end
       if ~isempty(self.VIforce),    FVI = self.VIforce;   else, disp('Need to define nonlinear contacts'),end
       
       if isempty(self.Minv); [As,Minv] = self.simulation_mats();
       else As = self.statematrix; Minv = self.Minv;
       end
        
       if self.excitation.steadystate ~= 1
       opts.Events = @myEvent;
       end
    

       % Getting the 'primer data' for the right boundary condition
       [~, yset] = ode45(@(t,y) sys2_persistent_bilinear(t,y,Minv,As,Fin,FVI,self),tsim(1:500),IC,odeset('RelTol',opts.RelTol))     ;
       global maxy;
       maxy = max(yset(:));
       
       n = numel(self.SimParams.Integrator);
       if n == 5
            if self.SimParams.Integrator == 'ode78'
            [tsim,ysim]  = ode78(@(t,y) sys2_persistent_bilinear(t,y,Minv,As,Fin,FVI,self),tsim,IC,opts)     ;
            elseif self.SimParams.Integrator == 'ode89'
            [tsim,ysim]  = ode89(@(t,y) sys2_persistent_bilinear(t,y,Minv,As,Fin,FVI,self),tsim,IC,opts)     ;
            elseif self.SimParams.Integrator == 'ode45'   
            [tsim,ysim]  = ode45(@(t,y) sys2_persistent_bilinear(t,y,Minv,As,Fin,FVI,self),tsim,IC,opts)     ;
            end
       elseif n == 6
            if self.SimParams.Integrator == 'ode113'
            [tsim,ysim]  = ode113(@(t,y) sys2_persistent_bilinear(t,y,Minv,As,Fin,FVI,self),tsim,IC,opts)     ;
            elseif self.SimParams.Integrator == 'ode15s'
            [tsim,ysim]  = ode15s(@(t,y) sys2_persistent_bilinear(t,y,Minv,As,Fin,FVI,self),tsim,IC,opts)     ;
            end
       end
       end


       % ------------------------------------------
       function [tsim,ysim] = simulate_noglobal(self, tsim, IC, opts)
       % ------------------------------------------
       % Simulate the system  (No udot updating)
       % ------------------------------------------
       % K - stiffness matrix of mechanical sysmet
       % M - mass matrix of mechanical system
       % C - Damping matrix of mechanical system
       % opts - options of ode integrator
         if nargin < 4
           opts = self.SimParams.opts;
       end 
       if nargin < 3   
            IC = self.SimParams.IC;
       end
       if nargin < 2
           tsim = self.SimParams.tsim;
       end
       if ~isempty(self.excitation), Fin = self.excitation; else, disp('Need to define forcing profle'),end
       if ~isempty(self.VIforce),    FVI = self.VIforce;   else,  disp('Need to define nonlinear contacts'),end
       
       if isempty(self.Minv); [As,Minv] = self.simulation_mats();
       else As = self.statematrix; Minv = self.Minv;
       end
        
       if self.excitation.steadystate ~= 1
       opts.Events = @myEvent;
       end
       
       % Getting the 'primer data' for the right boundary condition
       [~, yset] = ode45(@(t,y) sys2_global(t,y,Minv,As,Fin,FVI,self),tsim(1:500),IC,odeset('RelTol',opts.RelTol))     ;
       global maxy;
       maxy = max(yset(:));
       
       n = numel(self.SimParams.Integrator);
       if n == 5
            if self.SimParams.Integrator == 'ode78'
            [tsim,ysim]  = ode78(@(t,y) sys2_noglobal(t,y,Minv,As,Fin,FVI,self),tsim,IC,opts)     ;
            elseif self.SimParams.Integrator == 'ode89'
            [tsim,ysim]  = ode89(@(t,y) sys2_noglobal(t,y,Minv,As,Fin,FVI,self),tsim,IC,opts)     ;
            elseif self.SimParams.Integrator == 'ode45'   
            [tsim,ysim]  = ode78(@(t,y) sys2_noglobal(t,y,Minv,As,Fin,FVI,self),tsim,IC,opts)     ;
            end
       elseif n == 6
            if self.SimParams.Integrator == 'ode113'
            [tsim,ysim]  = ode113(@(t,y) sys2_noglobal(t,y,Minv,As,Fin,FVI,self),tsim,IC,opts)     ;
            elseif self.SimParams.Integrator == 'ode15s'
            [tsim,ysim]  = ode15s(@(t,y) sys2_noglobal(t,y,Minv,As,Fin,FVI,self),tsim,IC,opts)     ;
            end
       end
       end

   end
end