classdef wind3D
    %   This class initializes the environmental parameters 
    properties
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
            obj.M               = SIM.parameter('Value',28.97,'Unit','','Description','Molar mass of air');
            obj.temp            = SIM.parameter('Value',288.15,'Unit','K','Description','Ground temperature');
            obj.wind0           = SIM.parameter('Value',[0,0,0],'Unit','m/s','Description','Wind speed at ground level');
            obj.windZ           = SIM.parameter('Value',[0,0,0],'Unit','m/s','Description','Wind speed at 5000 ft');
        end
    end
end

