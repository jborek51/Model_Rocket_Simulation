%%  Load data from excel 
A = xlsread('Motor_Data.xlsx',2,'F3:H29');
F_t = A(:,2);
T_b = A(:,1);
m_m = A(3,3);
m_f = A(4,3);
%%  Save mat file
fpath = fullfile(fileparts(which('Motor_Data.xlsx')));
save([fpath '\E30.mat'],'m_m','m_f','F_t','T_b')
