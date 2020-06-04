classdef F15_Motor
    %   This class initializes the motor parameters for the Estes F15
    properties
        m = .102;                               %   kg - Total motor mass
        m_f = .06;                              %   kg - Mass of fuel
        F_th = [0,7.6380,12.2530,16.3910,...    %   N - Motor thrust curve 
            20.2100,22.7560,25.2600,23.0740,...
            20.8450,19.0930,17.5000,16.2250,...
            15.4270,14.9480,14.6270,15.7410,...
            14.7850,14.6230,14.3030,14.1410,...
            13.8190,13.3380,13.3340,13.0130,...
            9.3520,4.8950,0,0]';
        T_b = [0,.1480,.2280,.2940,.3530,...    %   s - Thrust curve time reference 
            .3820,.4190,.4770,.5200,.5930,...
            .6880,.8550,1.0370,1.2050,1.4230,...
            1.4520,1.5030,1.7360,1.9550,...
            2.2100,2.4940,2.7630,3.1200,...
            3.3820,3.4040,3.4180,3.4500,...
            3.5000]';
        eta = 1;                                %   Motor efficiency
    end
    properties (Dependent)
        t_b                                     %   s - Total motor burn time
    end
    methods 
        function obj = F15_Motor()              %   Constructor
            
        end
        function val = get.t_b(obj)             %   Get total burn time 
            val = obj.T_b(end);
        end
    end
end

