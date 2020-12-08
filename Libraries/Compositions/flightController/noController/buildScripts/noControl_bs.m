CONTROLLER = "noControl";

fltCtrl = CTR.FPID('deg','deg');

%% Save
saveBuildFile('fltCtrl',mfilename,'variant',"CONTROLLER");
