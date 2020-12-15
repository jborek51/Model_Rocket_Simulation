classdef wind3D
    %   This class initializes the environmental parameters 
    properties (SetAccess = private)
        grav
        density
        gndAlt
        railLength
        railElevation 
        railAzimuth 
        simDuration
        R
        M
        temp
        wind0
        windZ
    end
    properties (Dependent)
        windVecX
        windVecY
        windVecZ
    end
    methods 
        function obj = wind3D()
            obj.grav            = SIM.parameter('Value',9.8066,'Unit','m/s^2','Description','Acceleration due to gravity');
            obj.density         = SIM.parameter('Value',1.225,'Unit','kg/m^3','Description','Density at sea level');
            obj.gndAlt          = SIM.parameter('Value',0,'Unit','m','Description','Launch field altitude above sea level');
            obj.railLength      = SIM.parameter('Unit','m','Description','Launch rail length');
            obj.railElevation   = SIM.parameter('Unit','deg','Description','Launch rail elevation angle');
            obj.railAzimuth     = SIM.parameter('Unit','deg','Description','Launch rail azimuth angle');
            obj.simDuration     = SIM.parameter('Value',50,'Unit','s','Description','Launch rail azimuth angle');
            obj.R               = SIM.parameter('Value',8.3145,'Unit','J/(mol*K)','Description','Universal gas constant');
            obj.M               = SIM.parameter('Value',0.0289644,'Unit','','Description','Molar mass of air');
            obj.temp            = SIM.parameter('Value',288.15,'Unit','K','Description','Ground temperature');
            obj.wind0           = SIM.parameter('Value',[0,0,0],'Unit','m/s','Description','Wind vector at ground level');
            obj.windZ           = SIM.parameter('Value',[0,0,0],'Unit','m/s','Description','Wind vector at 5000 ft');
        end
        function val = get.windVecX(obj)
            V = [obj.wind0.Value(1) obj.windZ.Value(1)];
            val = SIM.parameter('Value',V,'Unit','(m/s)/m','Description','Wind X component look-up vector');
        end
        function val = get.windVecY(obj)
            V = [obj.wind0.Value(2) obj.windZ.Value(2)];
            val = SIM.parameter('Value',V,'Unit','(m/s)/m','Description','Wind Y component look-up vector');
        end
        function val = get.windVecZ(obj)
            V = [obj.wind0.Value(3) obj.windZ.Value(3)];
            val = SIM.parameter('Value',V,'Unit','(m/s)/m','Description','Wind Z component look-up vector');
        end
    end
end

