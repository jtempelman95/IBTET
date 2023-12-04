classdef SimParams

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
       Integrator
       opts
       IC
       tsim
    end
end
