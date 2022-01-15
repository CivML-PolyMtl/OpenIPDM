% Test Input:
y_perfect = [100,NaN,100,100,NaN,100,NaN,100,NaN,100,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN];
y_avg = [75,NaN,77,77,NaN,68,NaN,70,NaN,70,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN];
y_poor = [40,NaN,38,39,NaN,32,NaN,35,NaN,33,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN];
A = [1,1,0.5;0,1,1;0,0,1];
C = [1,0,0];
Q = [5.0625e-06,1.0125e-05,1.0125e-05;1.0125e-05,2.025e-05,2.025e-05;1.0125e-05,2.025e-05,2.025e-05];
R = 0;
Re = [5.86,0,5.94,5.79,0,5.79,0,5.79,0,3.03,0,0,0,0,0,0,0,0];
init_x = [0;0;0];
init_V = [0,0,0;0,0,0;0,0,0];
InpecBiase = [-0.28,0,-0.38,-0.91,0,-0.91,0,-0.91,0,-0.91,0,0,0,0,0,0,0,0];
OptmInsp = [0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0];
RU = 1;
InspBU = 0;
ObsYears = [1;0;1;1;0;1;0;1;0;1;0;0;0;0;0;0;0;0];
yearly = 2011:2028;
InspectorLabel = [5033851,NaN,5038915,5638405,NaN,5638405,NaN,5638405,NaN,6241477,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN];
StructureInd = 2;
ElementInd = 1;
TableEstimatedParameters = [0.0045,1.02,0.05,0.0245,0,0,0.1499,4];
RegressionModel = [];
AllAtt = [5,65,48.643155,100];
app.BatchMode = 0;
% Synthetic data config.
TableOfParameters_syn={0.005,1,0.1,0.1,NaN,NaN,0,4;0.0045,1.1,0.05,0.05,0,0,0.149,4};
InpecBiaseTrue=[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0];
ReTrue = [5.86,0,5.94,5.79,0,5.79,0,5.79,0,3.03,0,0,0,0,0,0,0,0,0];
QTrue = [6.25e-06,1.25e-05,1.25e-05;1.25e-05,2.5e-05,2.5e-05;1.25e-05,2.5e-05,2.5e-05];
SynDatabaseState = {[92.9400561665594,92.6342076771046,92.3106400730089,91.9640348261534,91.5924095422602,91.1952533480615,90.7711946263530,90.3207098022512,89.8372531805094,89.3179245173635,88.7664458391173,88.1805562193170,87.5540803909949,86.8831292471732,86.1673385455020,85.4049967181738,84.5903624490245,83.7216495730393,82.8051685603412,81.8430870111265,80.8347851895649,79.7825270833101,78.6872246312203,77.5479804906387,76.3667909446664,75.1481008834698,73.8956054396116,72.6063924593577,71.2783038187650,69.9101914641703,68.5021425432142,67.0524581473626,65.5572553843431,64.0150749676374,62.4225366151891,60.7807697948801,59.0870796916254,57.3400641537716,55.5413632995327,53.6931260886879,51.7915026842101,49.8360678068577,47.8275502531829,45.7710007194603,43.6696207198008,41.5222411154200,39.3210506358579,37.0692200350985,34.7705510130335,32.4207133560287,30.0191574124942,27.5666401649323,25.0633330565865,22.5110044626773,19.9061656204284,17.2467610056797,14.5290907713314,11.7476829239281,8.90278692302259,6.00291183266361;
    -0.298749994620669,-0.313782424225328,-0.334157683816114,-0.359106200579393,-0.384415434859324,-0.410888582429457,-0.436599293189451,-0.465966973836230,-0.502224711666259,-0.535064681579914,-0.567885568977270,-0.605751003706595,-0.648718450311186,-0.693209594761225,-0.738465121103913,-0.787379019007983,-0.842049264640543,-0.894070896082165,-0.938592972083874,-0.985479094973456,-1.03061912527685,-1.07393225128345,-1.11752990804088,-1.16091752876022,-1.20066886175927,-1.23550340393208,-1.27025976896665,-1.30876786608133,-1.34700380635671,-1.38859240782125,-1.42753070249222,-1.47198094848363,-1.51829990161598,-1.56751841507093,-1.61697842153579,-1.66667936013740,-1.72096897576996,-1.77349923610526,-1.82298753027666,-1.87435739025210,-1.92854759358407,-1.98274756198175,-2.03341349381506,-2.07893723544391,-2.12343772211704,-2.17318009034552,-2.22780391959033,-2.27478188610657,-2.32356230875016,-2.37575860322845,-2.42725138894772,-2.47793696909682,-2.52797173420695,-2.57781174490003,-2.63149749884136,-2.68822794003628,-2.74844165739344,-2.81486733149197,-2.87309440383173,-2.92636378078238;
    -0.0104604720540189,-0.0169785625095703,-0.0244135626236428,-0.0244727597814508,-0.0248699308359065,-0.0266246382716833,-0.0223012542061485,-0.0320157066661744,-0.0388073107071254,-0.0318276535793420,-0.0347407191480089,-0.0423124982578283,-0.0473222195024573,-0.0458464682362531,-0.0451779577862427,-0.0547932586315570,-0.0570569968032651,-0.0458851851375325,-0.0442569837729446,-0.0464809075177096,-0.0466283644611168,-0.0421276427177773,-0.0460272794859279,-0.0411155494938087,-0.0377691501347274,-0.0329526226427196,-0.0357800701606031,-0.0389956609550714,-0.0392691696654911,-0.0382987063456348,-0.0423199476384211,-0.0458681342508215,-0.0466833820276393,-0.0529347629140625,-0.0495102447040411,-0.0513261822951509,-0.0521756495139781,-0.0532248304596326,-0.0494183719395441,-0.0537530279511786,-0.0537335520957352,-0.0538982712803006,-0.0466232758378507,-0.0431578782896012,-0.0412974359701974,-0.0550195220649902,-0.0504927005804953,-0.0461800177718383,-0.0531473635991036,-0.0505687639468319,-0.0493812579531231,-0.0496072231176881,-0.0474417582599916,-0.0561876016957217,-0.0544001153460669,-0.0590654332632583,-0.0637487664400033,-0.0654872191575301,-0.0541703625417607,-0.0477094773587394]};
GraphMatrix = [false,false;false,false;false,false;false,false;false,false];
yearly = 1965:1982;
FigureID = 0;

%% Test #1: Average-health-state-real
% Test 1 verifies the performance of the deterioration model using data for a
% structural element in an average deterioration state
[loglik,~,~,~,~,~,~,~,~,~,~,~,~,x,s_Xsmooth,~]=KFsKF(app,y_avg, A, C, Q, R, Re,...
 init_x, init_V,InpecBiase,[],OptmInsp,RU,InspBU,ObsYears,yearly,...
 InspectorLabel,StructureInd,ElementInd,TableEstimatedParameters(1,end),...
 [],[],[],RegressionModel,AllAtt,TableEstimatedParameters,1,1);
assert(ceil(loglik) == -14);
assert(ceil(s_Xsmooth(1,1)) == 2);
assert(ceil(s_Xsmooth(2,1)) == 1);
assert(ceil(s_Xsmooth(3,1)) == 1);
assert(ceil(s_Xsmooth(1,end)) == 3);
assert(ceil(s_Xsmooth(2,end)) == 1);
assert(ceil(s_Xsmooth(3,end)) == 1);
assert(floor(x(1,1)) == 76);
assert(floor(x(2,1)) == -1);
assert(floor(x(3,1)) == -1);
assert(floor(x(1,end)) == 67);
assert(floor(x(2,end)) == -1);
assert(floor(x(3,end)) == -1);

%% Test #2: Minimum-health-state-real
% Test 2 verifies the performance of the deterioration model using data for a
% structural element in an minimum deterioration state
[loglik,~,~,~,~,~,~,~,~,~,~,~,~,x,s_Xsmooth,Vvalues]=KFsKF(app,y_poor, A, C, Q, R, Re,...
 init_x, init_V,InpecBiase,[],OptmInsp,RU,InspBU,ObsYears,yearly,...
 InspectorLabel,StructureInd,ElementInd,TableEstimatedParameters(1,end),...
 [],[],[],RegressionModel,AllAtt,TableEstimatedParameters,1,1);
assert(ceil(loglik) == -11);
assert(ceil(s_Xsmooth(1,1)) == 2);
assert(ceil(s_Xsmooth(2,1)) == 1);
assert(ceil(s_Xsmooth(3,1)) == 1);
assert(ceil(s_Xsmooth(1,end)) == 3);
assert(ceil(s_Xsmooth(2,end)) == 1);
assert(ceil(s_Xsmooth(3,end)) == 1);
assert(floor(x(1,1)) == 38);
assert(floor(x(2,1)) == -1);
assert(floor(x(3,1)) == -1);
assert(floor(x(1,end)) == 28);
assert(floor(x(2,end)) == -1);
assert(floor(x(3,end)) == -1);

%% Test #3: Perfect-health-state-real
% Test 3 verifies the performance of the deterioration model using data for a
% structural element in a perfect health state
[loglik,~,~,~,~,~,~,~,~,~,~,~,~,x,s_Xsmooth,Vvalues]=KFsKF(app,y_perfect, A, C, Q, R, Re,...
 init_x, init_V,InpecBiase,[],OptmInsp,RU,InspBU,ObsYears,yearly,...
 InspectorLabel,StructureInd,ElementInd,TableEstimatedParameters(1,end),...
 [],[],[],RegressionModel,AllAtt,TableEstimatedParameters,1,1);
assert(ceil(loglik) == -10);
assert(ceil(s_Xsmooth(1,1)) == 2);
assert(ceil(s_Xsmooth(2,1)) == 1);
assert(ceil(s_Xsmooth(3,1)) == 1);
assert(ceil(s_Xsmooth(1,end)) == 2);
assert(ceil(s_Xsmooth(2,end)) == 1);
assert(ceil(s_Xsmooth(3,end)) == 1);
assert(floor(x(1,1)) == 106);
assert(floor(x(2,1)) == -1);
assert(floor(x(3,1)) == -1);
assert(floor(x(1,end)) == 100);
assert(floor(x(2,end)) == -1);
assert(floor(x(3,end)) == -1);


%% Test #4: Average-health-state-Synthetic
% Test 4 verifies the performance of the deterioration model using data for a
% structural element in an average deterioration state
[loglik,MAcc,AccVal,AccST,PWithinCI,PWithinCITrue,Ebar,Ebarbar,...
BoundViolation,XtrueValue,XValue,XbtrueTr,XbbtrueTr,XestimateTr,~,~,...
Erbias, Erbarbias, AccValBias]=KFsKF([],y_avg, A, C, Q, R, Re,init_x, init_V,InpecBiase,InpecBiaseTrue,OptmInsp,RU,...
InspBU,ObsYears,yearly,InspectorLabel,StructureInd,ElementInd,...
TableOfParameters_syn{2,8},ReTrue,QTrue,SynDatabaseState,RegressionModel,AllAtt,...
TableOfParameters_syn,GraphMatrix,[],FigureID);
assert(ceil(loglik) == -15);
assert(floor(XestimateTr(1,1)) == 75);
assert(floor(XestimateTr(2,1)) == -1);
assert(floor(XestimateTr(3,1)) == -1);
assert(floor(XestimateTr(1,end)) == 63);
assert(floor(XestimateTr(2,end)) == -2);
assert(floor(XestimateTr(3,end)) == -1);
assert(floor(MAcc) == 17);
assert(Erbias(1,1)<=0.025);
assert(ceil(XbtrueTr(end))==-1);
%% Test #5: Minimum-health-state-Synthetic
% Test 5 verifies the performance of the deterioration model using data for a
% structural element in an minimum deterioration state
[loglik,MAcc,AccVal,AccST,PWithinCI,PWithinCITrue,Ebar,Ebarbar,...
BoundViolation,XtrueValue,XValue,XbtrueTr,XbbtrueTr,XestimateTr,~,~,...
Erbias, Erbarbias, AccValBias]=KFsKF([],y_poor, A, C, Q, R, Re,init_x, init_V,InpecBiase,InpecBiaseTrue,OptmInsp,RU,...
InspBU,ObsYears,yearly,InspectorLabel,StructureInd,ElementInd,...
TableOfParameters_syn{2,8},ReTrue,QTrue,SynDatabaseState,RegressionModel,AllAtt,...
TableOfParameters_syn,GraphMatrix,[],FigureID);
assert(ceil(loglik) == -13);
assert(floor(XestimateTr(1,1)) == 39);
assert(floor(XestimateTr(2,1)) == -1);
assert(floor(XestimateTr(3,1)) == -1);
assert(floor(XestimateTr(1,end)) == 27);
assert(floor(XestimateTr(2,end)) == -1);
assert(floor(XestimateTr(3,end)) == -1);
assert(floor(MAcc) == 54);
assert(Erbias(1,1)<=-0.0997);
assert(ceil(XbtrueTr(end))==-1);
%% Test #6: Perfect-health-state-Synthetic
% Test 6 verifies the performance of the deterioration model using data for a
% structural element in a perfect health state
[loglik,MAcc,AccVal,AccST,PWithinCI,PWithinCITrue,Ebar,Ebarbar,...
BoundViolation,XtrueValue,XValue,XbtrueTr,XbbtrueTr,XestimateTr,~,~,...
Erbias, Erbarbias, AccValBias]=KFsKF([],y_perfect, A, C, Q, R, Re,init_x, init_V,InpecBiase,InpecBiaseTrue,OptmInsp,RU,...
InspBU,ObsYears,yearly,InspectorLabel,StructureInd,ElementInd,...
TableOfParameters_syn{2,8},ReTrue,QTrue,SynDatabaseState,RegressionModel,AllAtt,...
TableOfParameters_syn,GraphMatrix,[],FigureID);
assert(ceil(loglik) == -12);
assert(floor(XestimateTr(1,1)) == 101);
assert(floor(XestimateTr(2,1)) == -1);
assert(floor(XestimateTr(3,1)) == -1);
assert(floor(XestimateTr(1,end)) == 94);
assert(floor(XestimateTr(2,end)) == -1);
assert(floor(XestimateTr(3,end)) == -1);
assert(floor(MAcc) == 10);
assert(Erbias(1,1)<=0.181);
assert(ceil(XbtrueTr(end))==-1);