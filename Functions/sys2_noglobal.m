function dx=sys2_noglobal(t,x,Minv,As,excitation, VIforce)
% ---------------------------
% Simulaiton of system 2
% ---------------------------
global udt3
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

             % -------------------------------------------------; 
             %              Linear system                       ;
             % -------------------------------------------------; 
              wext = omega                                      ;
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
              if lattice.excitation.steadystate == 1
              Force(end/2+FLOC)  = amp*sin(wext*t)                                    ;
              else
              fin = omega/2/pi;
              Ncy  = periods;
              f2 =amp*0.5*(heaviside(t)-...
              heaviside(t-2*pi/(2*pi*fin)*Ncy)).*(1-cos(2*pi*fin/Ncy*t)).*sin(2*pi*fin*t);
              Force(end/2+FLOC) = f2;
              end 
              
              % +++++++++++++++++++++++++
              % LINEAR SYSTEM
              % +++++++++++++++++++++++++
              Force(end/2+1:end) = Minv* Force(end/2+1:end)     ;
              dx = As*x + Force        ;     
             
            %++++++++++++++++++++++++++++++++++++++++++
            %             HERTZIAN FORCES                                       
            % +++++++++++++++++++++++++++++++++++++++++
            % global kc d rescoeff Index VIsites UDOTMINUS
            % global Index UDOTMINUS
            r = rescoeff;
            % Loop over all impact indicies
            %if t == 0; Index = ones(1,S);   end
            
            Index(VIsites) = 1;
            
            objcounter.jj = 0;
            if t <= TimeVI
                VI_lst = zeros(numel(VIsites),1e6);
     %           objcounter.jj = 0;
            else
                currentcount = objcounter.jj ;
      %          objcounter.jj = currentcount +1;
      %          impact_counter.j = impact_counter.j + 1;
            end
            if t > TimeVI
            ii =1;
       %     impact_counter.j
            for SITE = VIsites
                    % Check difference between impact locations and`
                    % velocities
                    delta = x(SITE+1) - x(SITE)                 ;
                    ddelta = dx(SITE+1) - dx(SITE)              ;
                    
                    % ----------------------------------------------------
                    % Check if contact is initiated. If no contact, set to
                    % 1, If contact, start adding to string and only use
                    % the first term
                    % ----------------------------------------------------
                     if  abs(delta) < d(SITE)           % e.g., if no contact
                         Index(SITE) = 1    ;
                         iscontact = 1      ;
                     elseif  -delta > - d(SITE)         % e.g., if no contact
                         Index(SITE) = 1    ;
                         iscontact = 1      ;
                     else 
                         iscontact = 0      ;
                     end
                    
                   %  VI_lst(ii,j) = iscontact;
                    % if iscontact(i) ==1 & iscontact(i-)
                    
                   if ddelta <= 0.00001; ddelta = .00001; end 
                   UDOTMINUS(ii,Index(SITE)) = ddelta;  
                    
                   if rem(SITE,2) ~= 0 % Check if even number (host site)
                       if sided ==1;
                        fnl_hertz =             -kc*( max(-delta-d(SITE),0) .^(3/2) ) ; 
                       elseif sided == -1
                           fnl_hertz =          -kc*( -max(delta-d(SITE),0) .^(3/2) ) ; 
                       elseif sided == 2
                        fnl_hertz =             -kc*( max(-delta-d(SITE),0) .^(3/2) - max(delta-d(SITE),0) .^(3/2) ) ;
                       end
                   else
                       %fnl_hertz =         -kc*( - max(delta-d(SITE),0) .^(3/2) ) ;
                   end
                    fnl_hertz =         real(fnl_hertz)             ; 
                    fnl_hertz =         fnl_hertz*(1 + 3*(1-r).*ddelta./2./ddelta); % Dissipation term
                    Fvec2 =             zeros(size(As,1),1)         ;
                    
                    % Force applied to host mass
                    Fvec2(SITE+S) =     fnl_hertz                   ;
                    
                    % Force applied to resonator mass
                    Fvec2(SITE+S+1) =   -fnl_hertz                  ;   
                    
                    % Hertzian forcing vector 
                    Fvec2(end/2+1:end) = Minv* Fvec2(end/2+1:end)   ;  
                    dx =  dx + Fvec2                                ;
                    Index(SITE) = Index(SITE) + 1                   ;
                    ii = ii + 1;
            end
            end
end
