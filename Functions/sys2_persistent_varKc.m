function dx=sys2_persistent(t,x,Minv,As,excitation, VIforce,lattice)
% ---------------------------
% Simulaiton of system 2
% ---------------------------
persistent Index_pers UDOTMINUS_pers
%global udt2

            % global FLOC eigenfreqs tau t0 steadystate TimeVI sided\
            omega       = excitation.omega;
            periods     = excitation.periods;
            amp         = excitation.amp;
            steadystate = excitation.steadystate;
            FLOC        = excitation.FLOC;
            tau         = excitation.tau;
            t0          = excitation.t0;
            TimeVI      = VIforce.TimeVI;
            rescoeff    = VIforce.rescoeff;
            d           = VIforce.d;
            kc          = VIforce.kc;
            VIsites     = VIforce.VIsites;
            sided       = VIforce.sided;
            config      = lattice.config;

              wext = omega;
              S = length(Minv)                                  ;
              F =0                                              ;
              Force = zeros(size(As,1),1)                       ;
            
              % Sinusoid forcing
              if t <periods*2*pi/wext
               Force(end/2+FLOC) = amp*sin(wext*t)              ;
              end

              % Gaussian tone burst
              St = amp*exp(-sqrt(-1)*omega*t-(t-t0).^2/tau^2)         ;
              St2 = amp*exp(-sqrt(-1)*omega*.75*t-(t-t0).^2/tau^2)    ;
              %%% COMMENT TO USE SINUSOID %%%
              if steadystate == 1
             %   Force(end/2+FLOC)  = amp*sin(wext*t)                                    ;
              else
             %   Force(end/2+FLOC) = real(St)+0*real(St2)  ;%+ 20*(rand-.5)*abs(real(St));
              end       
              fin = omega/2/pi;
              Ncy  = periods;
              f2 =amp*0.5*(heaviside(t)-...
              heaviside(t-2*pi/(2*pi*fin)*Ncy)).*(1-cos(2*pi*fin/Ncy*t)).*sin(2*pi*fin*t);
              Force(end/2+FLOC) = f2;
              
              % +++++++++++++++++++++++++
              % LINEAR SYSTEM
              % +++++++++++++++++++++++++
              Force(end/2+1:end) = Minv* Force(end/2+1:end)     ;
              dx = As*x + Force                                 ;
             
            %++++++++++++++++++++++++++++++++++++++++++
            %             HERTZIAN FORCES                                       
            % +++++++++++++++++++++++++++++++++++++++++
            % global kc d rescoeff Index_pers VIsites UDOTMINUS_pers
            
            r = rescoeff;
            % Loop over all impact indicies
            %if t == 0; Index_pers = ones(1,S);   end
            
            if t < eps;
            Index_pers(VIsites) = 1;
            end
            
            objcounter.jj = 0;
            if t <= TimeVI
                VI_lst = zeros(numel(VIsites),1e6);
                objcounter.jj = 0;
            else
                currentcount = objcounter.jj ;
                objcounter.jj = currentcount +1;
            end
            if t > TimeVI
            %ii =1;
            for ii = 1:numel(VIsites)
                if numel(kc)>1, KC = kc(ii);else,KC=kc;end
                    SITE = VIsites(ii);
                    % Check difference between impact locations and`
                    % velocities (The odd number must be higher than the
                    % even number for this to make sence. Check VIsties if
                    % not sure - currenlty set up so VIsites are even
                    % numbers
                        
                    delta = x(SITE+1) - x(SITE)                 ;
                    ddelta = dx(SITE+1) - dx(SITE)              ;
                    %if config ==3; delta = -delta; ddelta = -ddelta;end;
                        
                    % ----------------------------------------------------
                    % Check if contact is initiated. If no contact, set to
                    % 1, If contact, start adding to string and only use
                    % the first term
                    % ----------------------------------------------------
                     if  abs(delta) < d(SITE)           % e.g., if no contact
                         Index_pers(SITE) = 1    ;
                     elseif  -delta > - d(SITE)         % e.g., if no contact
                         Index_pers(SITE) = 1    ;
                     end

                   if ddelta <= 0.00001; ddelta = .00001; end 
                   UDOTMINUS_pers(ii,Index_pers(SITE)) = ddelta;  

                       if sided ==1;
                        fnl_hertz =             -KC*( max(-delta-d(SITE),0) .^(3/2) ) ; 
                       elseif sided == -1
                           fnl_hertz =          -KC*( -max(delta-d(SITE),0) .^(3/2) ) ; 
                       elseif sided == 2
                        fnl_hertz =             -KC*( max(-delta-d(SITE),0) .^(3/2) - max(delta-d(SITE),0) .^(3/2) ) ;
                       end
               
                    % if config ==3; fnl_hertz = -fnl_hertz;end
                    fnl_hertz =         real(fnl_hertz)             ; 
                    fnl_hertz =         fnl_hertz*(1 + 3*(1-r).*ddelta./2./UDOTMINUS_pers(ii,1)); % Dissipation term
                    Fvec2 =             zeros(size(As,1),1)         ;
                    
                    % Force applied to host mass
                    Fvec2(SITE+S) =     fnl_hertz                   ;
                    
                    % Force applied to resonator mass
                    Fvec2(SITE+S+1) =   -fnl_hertz                  ;   
                    
                    % Hertzian forcing vector 
                    Fvec2(end/2+1:end) = Minv* Fvec2(end/2+1:end)   ;  
                    Index_pers(SITE) = Index_pers(SITE) + 1                   ;
                 %   ii = ii + 1;

            dx =  dx + Fvec2  ;
            end
        end
end
