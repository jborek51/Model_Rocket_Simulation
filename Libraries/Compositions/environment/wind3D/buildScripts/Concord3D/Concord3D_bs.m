%% Set up environment

Env = ENV.wind3D;

Env.gndAlt.setValue(215,'m');
Env.railLength.setValue(1,'m');
Env.railElevation.setValue(90,'deg');
Env.railAzimuth.setValue(0,'deg');

ENVIRONMENT = 'wind3D';

saveBuildFile('Env',mfilename,'variant',["ENVIRONMENT"]);
