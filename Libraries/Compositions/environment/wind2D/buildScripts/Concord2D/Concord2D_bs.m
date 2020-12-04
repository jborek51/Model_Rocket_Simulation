%% Set up environment

Env = ENV.wind2D;

Env.gndAlt.setValue(215,'m');
Env.railLength.setValue(1,'m');
Env.railElevation.setValue(90,'deg');

ENVIRONMENT = 'wind2D';

saveBuildFile('Env',mfilename,'variant',["ENVIRONMENT"]);
