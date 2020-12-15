classdef motor
    %   This class initializes the motor parameters for the Estes F15
    properties
        motorName
        efficiency
        motorMass
        fuelMass
        thrustVec
        burnTimeVec
        position
    end
    properties (Dependent)
        burnTime
    end
    methods 
        function obj = motor()
            obj.motorName   = SIM.parameter('Description','Filename for the motor selection');
            obj.efficiency  = SIM.parameter('Value',1,'Unit','','Description','Motor efficiency scalar gain');
            obj.position    = SIM.parameter('Value',[0;0;0],'Unit','m','Description','Motor attachment location vector from nose');
        end
                
        function setMotorName(obj,val,units)
            if ~endsWith(val,'.mat')
                val = [val '.mat'] ;
            end
            obj.motorName.setValue(val,units);
        end
        
        function val = get.burnTime(obj)
            val = SIM.parameter('Value',obj.burnTimeVec.Value(end-2),'Unit','s','Description','Total motor burn time');
        end
        
        function plotThrustCurve(obj)
            figure; hold on; grid on;
            plot(obj.burnTimeVec.Value,obj.thrustVec.Value,'b-')
            xlabel('Time [s]'); ylabel('Thrust [N]');
        end
    end
end

