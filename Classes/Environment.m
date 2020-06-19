classdef Environment
    %   This class initializes the environmental parameters 
    properties
        g = 9.80665;                            %   m/s^2 - Gravitational acceleration
        rho = 1.225;                            %   kg/m^3 - Air density
        z_0 = 0;                                %   m - Initial altitude above sea level
        z_r = 36*.0254;                         %   in -> m - Rail length
        phi_0 = 0*pi/180;                       %   deg -> rad - Launch angle off vertical
        t_f = 15;                               %   s - Simulation end time
        R = 8.314;                              %   N.m/mol.K - Universal gas constant
        M = 28.97;                              %   kg/mol - Molar mass of air
        T_0 = 288.15;                           %   K = Base temperature
        v_w1 = 5;                               %   m/s - Wind speed at ground level 
        v_w2 = 5;                               %   m/s - Wind speed at altitude z_w
        z_w = 5000/3.2808;                      %   ft -> m - Wind speed altitude reference
    end
    methods 
        function obj = Environment()            %   Constructor
        end
    end
end

