%%  Save Output Variable Structure
if computer == 1
    results_path = 'C:\Users\John Jr\Dropbox (UNC Charlotte)\Projects\Rocket\Model_Rocket_Simulation\Data Files\Results';
else
    results_path = 'C:\Users\jbore\Dropbox (UNC Charlotte)\Projects\Rocket\Model_Rocket_Simulation\Data Files\Results';
end
cond = 1;   test = 0;
while cond == 1
    filename = strcat('1-DOF_Output_',num2str(test),'.mat');
    if exist(filename,'file')
        test = test + 1;
    else
        cond = 0;
    end
end
save(strcat(results_path,'\1-DOF_Output_',num2str(test),'.mat'),'Sim')