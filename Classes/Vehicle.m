classdef Vehicle
    %   This class initializes the vehicle parameters 
    properties
        m = 2.74*0.453592;                      %   lb -> kg - Vehicle mass
        L = 44*0.0254;                          %   in -> m - Vehicle length
        D = 4*0.0254;                           %   in -> m - Airframe diameter 
        Cd = .55;                               %   Drag coefficient 
        CG_1 = 21.70*.0254;                     %   in -> m - Center of gravity on pad
        CG_2 = 17.60*.0254;                     %   in -> m - Center of gravity after burnout
        CP = 19.78*.0254;                       %   in -> m - Center of Pressure
    end
    properties (Dependent)
        A_ref                                   %   m^2 - Reference area
        L_CG                                    %   m - Change in CG position 
    end
    methods 
        function obj = Vehicle()                %   Constructor
            
        end
        function val = get.A_ref(obj)           %   Get reference area 
            val = pi*obj.D^2/4;
        end
        function val = get.L_CG(obj)            %   Get change in CG 
            val = obj.CG_1-obj.CG_2;
        end
    end
end

