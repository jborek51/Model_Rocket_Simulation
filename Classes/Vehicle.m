classdef Vehicle
    %   This class initializes the vehicle parameters 
    properties
        m = 2.74*0.453592;                      %   lb -> kg - Vehicle mass
        L = 44*0.0254;                          %   in -> m - Vehicle length
        D = 4*0.0254;                           %   in -> m - Airframe diameter 
        A_ref                                   %   m^2 - Reference area
        Cd = .55;                               %   Drag coefficient 
        CG_1 = 21.70*.0254;                     %   in -> m - Center of gravity on pad
        CG_2 = 17.60*.0254;                     %   in -> m - Center of gravity after burnout
        L_CG                                    %   m - Change in CG position 
        CP = 19.78*.0254;                       %   in -> m - Center of Pressure
    end
    methods 
        function obj = Vehicle()                %   Constructor
            obj.A_ref = pi*obj.D^2/4;           %   Set reference area 
            obj.L_CG = obj.CG_1-obj.CG_2;       %   Set CG variation length 
        end
    end
end

