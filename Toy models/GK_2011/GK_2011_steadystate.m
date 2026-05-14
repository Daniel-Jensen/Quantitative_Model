function [ys, params, check] = GK_2011_steadystate(ys_init, exo_ss, M_, options_)

%% DO NOT CHANGE THIS PART.
%%
%% Here we load the values of the deep parameters in a loop.
NumberOfParameters = M_.param_nbr;
for ii = 1:NumberOfParameters
    paramname = M_.param_names{ii};
    eval([ paramname ' = M_.params(' int2str(ii) ');']);
    eval([ 'params.' paramname ' = M_.params(' int2str(ii) ');']);
end
% initialize indicator
check = 0;


%% THIS BLOCK IS MODEL SPECIFIC.
%%
%% Here the user has to define the steady state.
%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
starting.K0 = K0;
starting.L0 = L0;
starting.chi0 =  chi0;
starting.lambda0  =  lambda0;
starting.omega0 = omega0;

%Setting target moments
moments.L_mom    =   L_mom;
moments.RkmR_mom =   RkmR_mom;
moments.phi_mom  =   phi_mom;
switch_print    =   'no';
switches.switch_print    =  switch_print;  %'test' for too many details, 'details', 'yes', 'no'
%Calibrating the parameters to hit the moments
switch switch_print
    case 'details'
        options     =   optimset('Display','iter');
    case 'test'
        options     =   optimset('Display','iter');
    otherwise
        options     =   optimset('Display','off');
end;

[XX,diff,exitf]     =   fsolve(@f_mom,[omega0 lambda0 chi0],options,params,starting,moments,switches);


omega  =   XX(1);
lambda =   XX(2);
chi    =   XX(3);

%params   =   params;
params.omega    =   omega;
params.lambda   =   lambda;
params.chi      =   chi;

[XX_ss,fval,exitf]  =   fsolve(@f_KL,[K0 L0],options,params);

K    =   log(XX_ss(1));
L    =   log(XX_ss(2));


%Calculating the equilibrium K and N values given the parameters
switch switches.switch_print
    case 'test'
        options     =   optimset('Display','iter');
    otherwise
        options     =   optimset('Display','off');
end;

%prices value from author simulations
%K    =  log(5.661570773251615);
Kf   =  K;
%L    =  log(0.333333333333694);
Lf   =  L;
Pm   =  log((epsilon-1)/epsilon);
X    =  log(1/((epsilon-1)/epsilon));
D    =  log(1);
Q    =  log(1);
Qf   =  log(1);
delta= log(deltai);
deltaf= delta;
%Macrovariables
Y   =  log(exp(K)^alfa*exp(L)^(1-alfa));
Y_ss= Y;
Yf  =  Y;
Ym  =  Y;
Ymf =  Y;
I   =  log(deltai*exp(K));
I_ss= exp(I);
If  =  I;
In  =  0;
In_f=  0;
G   =  log(G_over_Y*exp(Y));
G_ss= exp(G);
Gf  =  G;
C   =  log(exp(Y)-exp(I)-exp(G));
Cf  =  C;
varrho  =   log((1-betta*hh)*((1-hh)*exp(C))^(-sig));
varrhof = varrho;
Lambda  =  log(1);
Lambdaf =  Lambda;
R       =  log(1/betta);
Rf      =  R;
Rk      =  log(exp(Pm)*alfa*exp(Y)/exp(K)+1-deltai);
Rkf     =  Rk;
RkmR    =  log(exp(Rk)-exp(R));
RkmRf   =  RkmR;
w       =  log(exp(Pm)*(1-alfa)*exp(Y)/exp(L));
wf      =  w;
VMPK    =  log(exp(Pm)*alfa*exp(Y)/exp(K));
VMPKf   =  VMPK;
U       =  log(1);
Uf      =  log(1);
if sig==1
    Welf    =   (log((1-hh)*exp(C))-chi*exp(L)^(1+varphi)/(1+varphi))/(1-betta);
    Welff=Welf;
else
    Welf    =   (((1-hh)*exp(C))^(1-sig)/(1-sig)-chi*exp(L)^(1+varphi)/(1+varphi))/(1-betta);
    Welff=Welf;
end
%Variable capacity utilization parameters
b       =    exp(Pm)*alfa*exp(Y)/exp(K);
delta_c=exp(delta)-b/(1+zetta);
%Intermediaries
aa     =   lambda*betta*theta*exp(RkmR);
bb     =   -(1-theta)*(lambda-betta*exp(RkmR));
cc     =   (1-theta);
phi    =   log((-bb-sqrt(bb^2-4*aa*cc))/(2*aa));
phif   =   phi;
z      =   log(exp(RkmR)*exp(phi)+exp(R));
zf     =   z;
x      =   z;
xf     =   x;
nu     =   log(((1-theta)*betta*exp(RkmR))/(1-betta*theta*exp(x)));
nuf    =   nu;
eta    =   log((1-theta)/(1-betta*theta*exp(z)));
etaf   =   eta;
N      =   log(omega*exp(K)/(1-theta*(exp(RkmR)*exp(phi)+exp(R))));
Nf     =   N;
Ne     =   log(theta*exp(z)*exp(N));
Nef    =   Ne;
Nn     =   log(omega*exp(K));
Nnf    =   Nn;
Keff   =   K;
Kefff  =   Keff;   %Effective capital
prem   =   log(exp(Rk)/exp(R));
premf  =   prem;     %Premium

%Pricing variables
F   =   log(exp(Y)*exp(Pm)/(1-betta*gam));
Z   =   log(exp(Y)/(1-betta*gam));
infl  =   0;
inflstar =0;
i   =   R;
%MMB varaibles
interest= 400*(i-log(1/betta));
inflation  = 0;
inflationq  = 0;
output     = 0;
outputgap  = 0;


%shocks

a=0;
ksi=0;
g=0;

e_a     =   0;
e_ksi   =   0;
e_g     =   0;
e_Ne    =   0;
e_i     =   0;

%%
%% END OF THE MODEL SPECIFIC BLOCK.


%% DO NOT CHANGE THIS PART.
%%
%% Here we define the steady state values of the endogenous variables of
%% the model.
%%
%% DO NOT CHANGE THIS PART.
%%
%% Here we define the steady state values of the endogenous variables
%% and parameters of the model.
%% end own model equations

params=NaN(NumberOfParameters,1);
for iter = 1:length(M_.params) %update parameters set in the file
    eval([ 'params(' num2str(iter) ',1) = ' M_.param_names{iter} ';' ])
end

NumberOfEndogenousVariables = M_.orig_endo_nbr; %auxiliary variables are set automatically
for ii = 1:NumberOfEndogenousVariables
    varname = M_.endo_names{ii};
    eval(['ys(' int2str(ii) ',1) = ' varname ';']);
end                                                           % End of the loop.
%%
%% END OF THE SECOND MODEL INDEPENDENT BLOCK.