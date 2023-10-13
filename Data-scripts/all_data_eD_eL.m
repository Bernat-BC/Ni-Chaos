% Select the folder GAIO-master. This only needs to be done once every
% session and can be commented after. 
% Instead, you can also add the GAIO-master folder to your path.

addpath(genpath(uigetdir(matlab.desktop.editor.getActiveFilename)));

clear;
F = findall(0,'type','figure','tag','TMWWaitbar');
delete(F);

E1 = linspace(0.000251, 0.00157, 25);
E2 = linspace(0.00157, 0.004, 25);
E = cat(2, E1, E2);                                                        % This choice is made for plotting purposes, E can be chosen as a linspace from any two e_min, e_max.

LYAPS_m = zeros(length(E),1);
LYAPS_w = zeros(length(E), 1);

Lyap_exp = zeros(length(E), 1);
Lyap_con = zeros(length(E), 1);
Escape_exp = zeros(length(E), 1);
Escape_con = zeros(length(E), 1);
Lambdas_exp = zeros(length(E), 1);
Lambdas_con = zeros(length(E), 1);

%GAIO parameters:
depth = 14;                                                                 % We use 2^depth boxes in total
n = 10;                                                                    % Number of particles in a box
noise = 500;                                                               % Noise realisations considered                                  

a = 3.83;                                                                  % Logistic parameter

f = @(x,omega)(a*x(:,1).*(1-x(:,1)) + omega(:));                           % Logistic function
L = @(x) log(abs(a*(1-2*x)));                                              % Lyapunov Exponent fn


% We define the indicator function for the regions Mexp and Mcont.
Inexp = @IN_exp_two_three_sure;                                            %for all w, it must expand in three iterates AND in two it must also expand. This ensures far away enough from attractor
Incon = @IN_con_compl_two;


% Other definitions for these regions include: (see bottom of script)
% Inexp = @IN_exp_sure;               
% Incon = @IN_con_compl_sure;                                              %exists w st third iterate is contracting.

% Inexp = @IN_exp_compl_sure;
% Incon = @IN_con_sure;

% Incon = @IN_con_compl_two;
% Inexp = @IN_exp_two_sure_three_compl;

% Inexp = @IN_exp_compl_two_or_three;
% Incon = @IN_con_two_or_three;

h = waitbar(0,' Computing...');
for t = 1:length(E)
    eps = E(t);
    max_val = f(1/2, eps);                                                 % Left boundary of invariant region
    min_val = f(f(1/2, eps), -eps);                                        % Right boundary of invariant region

    x = linspace(-1,1,n);                                                  % To create the boxes and...
    c = (max_val + min_val)/2;                                             % .
    r = (max_val - min_val)/2;                                             % .
    B = c + r.*x(:);                                                       % .  
    tree = Tree(c, r);                                                     % ...generate the box collection
    sd = 8;
    for i = 1:depth
        tree.set_flags('all', sd);                                         % Flag all leaves
        tree.subdivide;                                                    % Subdivide flagged leaves
    end

    Pt = tpmatrix1d(tree, f, B, depth, 0, noise, eps);                     % Trans prob matrix
    [xh,lambda] = eigs(Pt, 1,'LR');                                        % max eigenvalue and eigenvector


    sm = abs(xh(:,1))/norm(xh(:,1),1);                                     % Stationary measure
    Xs = linspace(min_val, max_val, length(sm));                           % xs where SM is evaluated
    Lyap = L(Xs)*sm;                                                       % (global) Lyapunov exponent

    mask = sparse(Inexp(Xs', eps));
    Mask = sparse(mask*mask');
    Pt_cond = Mask.*Pt;

    mask2 = sparse(Incon(Xs', eps));
    Mask2 = sparse(mask2*mask2');
    Pt_cond2 = Mask2.*Pt;
    
    [xh_cond,lambda_cond] = eigs(sparse(Pt_cond), 1,'LR');                 % Max eval and evec
    [xh_cond2,lambda_cond2] = eigs(sparse(Pt_cond2), 1,'LR'); 

    qsm = abs(xh_cond(:,1))/norm(xh_cond(:,1),1);                          % Quasistat measure IN
    qsm2 = abs(xh_cond2(:,1))/norm(xh_cond2(:,1),1);                       % Quasistat measure OUT

    [xh_cond_right,lambda_cond_right] = eigs(sparse(...
                                            transpose(Pt_cond)),1,'LR');   % Right ev for em
    ev_cond_right = abs(xh_cond_right(:,1))/norm(xh_cond_right(:,1),1);    % norm 
    [xh_cond_right2,lambda_cond_right2] = eigs(sparse(...
                                           transpose(Pt_cond2)),1,'LR');
    ev_cond_right2 = abs(xh_cond_right2(:,1))/norm(xh_cond_right2(:,1),1);

    qem = ev_cond_right.* qsm/(norm(ev_cond_right.* qsm, 1));              % Qerg measure
    qem2 = ev_cond_right2.* qsm2/(norm(ev_cond_right2.* qsm2, 1));               

    escape_t = 1/(1-lambda_cond);                                          % Expected Escape time from region "In"
    escape_t2 = 1/(1-lambda_cond2);                                        % Expected Escape time from region "OUT"
    
    Lyap_cond = L(Xs)*qem;                                                 % Conditioned Lyapunov exponent
    Lyap_cond2 = L(Xs)*qem2;

    L_weighted = (1/(escape_t + escape_t2))*(Lyap_cond * escape_t...
                                            + Lyap_cond2 * escape_t2);
                                                                           
    LYAPS_m(t) = Lyap;
    LYAPS_w(t) = L_weighted;
    
    Lyap_exp(t) = Lyap_cond;
    Lyap_con(t) = Lyap_cond2;
    
    Lambdas_exp(t) = lambda_cond;
    Lambdas_con(t) = lambda_cond2;
    
    Escape_exp(t) = escape_t;
    Escape_con(t) = escape_t2;
    
    waitbar(t/length(E), h);
end
close(h);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% For plotting purposes (see Plots/plots_eD_eL.m):
E_bis = cat(2,E1(1:2:length(E1)), E2);
Lyap_con_bis = cat(1,Lyap_con(1:2:length(E1)),...
                     Lyap_con(length(E1)+1:length(Lyap_con)));
Lyap_exp_bis = cat(1,Lyap_exp(1:2:length(E1)),...
                     Lyap_exp(length(E1)+1:length(Lyap_exp)));
Lyap_w_bis = cat(1,LYAPS_w(1:2:length(E1)),...
                   LYAPS_w(length(E1)+1:length(LYAPS_w)));


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Here we define several regions one can choose from when defining Inexp
% and Incon

%Surely Expanding on 3rd iterate
function IN_exp_sure = IN_exp_sure(x, eps)                                 % Definition of the region IN (expanding = True)
a = 3.83;
f = @(x,omega)(a*x(:,1).*(1-x(:,1)) + omega);
% If one of them is contracting on 3rd, then it's not expanding anymore.
aux = (abs(a*(1-2.*x).*a.*(1-2.*f(x, -eps)).*a.*(1-2.*f(f(x, -eps), +eps))) < 1)...
   +(abs(a*(1-2.*x).*a.*(1-2.*f(x, +eps)).*a.*(1-2.*f(f(x, +eps), -eps))) < 1)...
   +(abs(a*(1-2.*x).*a.*(1-2.*f(x, -eps)).*a.*(1-2.*f(f(x, -eps), -eps))) < 1)...
   +(abs(a*(1-2.*x).*a.*(1-2.*f(x, +eps)).*a.*(1-2.*f(f(x, +eps), +eps))) < 1);
IN_exp_sure = (aux < 1);
end

% Contracting (non-surely): Complement to Surely Expanding on 3rd iterate
function IN_con_compl_sure = IN_con_compl_sure(x, eps)                
a = 3.83;
f = @(x,omega)(a*x(:,1).*(1-x(:,1)) + omega);
% If one of them is contracting, then it's contracting.
aux = (abs(a*(1-2.*x).*a.*(1-2.*f(x, -eps)).*a.*(1-2.*f(f(x, -eps), +eps))) < 1)...
    +(abs(a*(1-2.*x).*a.*(1-2.*f(x, +eps)).*a.*(1-2.*f(f(x, +eps), -eps))) < 1)...
    +(abs(a*(1-2.*x).*a.*(1-2.*f(x, -eps)).*a.*(1-2.*f(f(x, -eps), -eps))) < 1)...
    +(abs(a*(1-2.*x).*a.*(1-2.*f(x, +eps)).*a.*(1-2.*f(f(x, +eps), +eps))) < 1);
IN_con_compl_sure = 1 - (aux < 1);
end

%Surely contracting on third iterate
function IN_con_sure = IN_con_sure(x, eps)               
a = 3.83;
f = @(x,omega)(a*x(:,1).*(1-x(:,1)) + omega);
% If one of them is expanding, then it's not contracting anymore.
aux = (abs(a*(1-2.*x).*a.*(1-2.*f(x, -eps)).*a.*(1-2.*f(f(x, -eps), +eps))) > 1)...
    +(abs(a*(1-2.*x).*a.*(1-2.*f(x, +eps)).*a.*(1-2.*f(f(x, +eps), -eps))) > 1)...
    +(abs(a*(1-2.*x).*a.*(1-2.*f(x, -eps)).*a.*(1-2.*f(f(x, -eps), -eps))) > 1)...
    +(abs(a*(1-2.*x).*a.*(1-2.*f(x, +eps)).*a.*(1-2.*f(f(x, +eps), +eps))) > 1);
IN_con_sure = (aux < 1);
end

%Expanding (non-surely): complement to surely expanding on third iterate
function IN_exp_compl_sure = IN_exp_compl_sure(x, eps)               
a = 3.83;
f = @(x,omega)(a*x(:,1).*(1-x(:,1)) + omega);
% If one of them is expanding, then it's not contracting anymore.
aux = (abs(a*(1-2.*x).*a.*(1-2.*f(x, -eps)).*a.*(1-2.*f(f(x, -eps), +eps))) > 1)...
    +(abs(a*(1-2.*x).*a.*(1-2.*f(x, +eps)).*a.*(1-2.*f(f(x, +eps), -eps))) > 1)...
    +(abs(a*(1-2.*x).*a.*(1-2.*f(x, -eps)).*a.*(1-2.*f(f(x, -eps), -eps))) > 1)...
    +(abs(a*(1-2.*x).*a.*(1-2.*f(x, +eps)).*a.*(1-2.*f(f(x, +eps), +eps))) > 1);
IN_exp_compl_sure = 1-(aux < 1);
end

%Expanding (non-surely): but surely on 2
function IN_exp_two_sure_three_compl = IN_exp_two_sure_three_compl(x, eps)               
a = 3.83;
f = @(x,omega)(a*x(:,1).*(1-x(:,1)) + omega);
% If one of them is expanding, then it's not contracting anymore.
aux3 = (abs(a*(1-2.*x).*a.*(1-2.*f(x, -eps)).*a.*(1-2.*f(f(x, -eps), +eps))) > 1)...
    +(abs(a*(1-2.*x).*a.*(1-2.*f(x, +eps)).*a.*(1-2.*f(f(x, +eps), -eps))) > 1)...
    +(abs(a*(1-2.*x).*a.*(1-2.*f(x, -eps)).*a.*(1-2.*f(f(x, -eps), -eps))) > 1)...
    +(abs(a*(1-2.*x).*a.*(1-2.*f(x, +eps)).*a.*(1-2.*f(f(x, +eps), +eps))) > 1);
aux2 = ((abs(a*(1-2.*x).*a.*(1-2.*f(x, +eps))) < 1)...
   +(abs(a*(1-2.*x).*a.*(1-2.*f(x, -eps)))<1));

aux = (aux2<1).*(1-(aux3 < 1));

IN_exp_two_sure_three_compl = aux ;
end

%Surely Expanding on third or second iterates
function IN_exp_two_three_sure = IN_exp_two_three_sure(x, eps)                
a = 3.83;
f = @(x,omega)(a*x(:,1).*(1-x(:,1)) + omega);
% If one of them is contracting in the 2nd or 3rd iterate, then it's not Expanding anymore.
aux3 = ((abs(a*(1-2.*x).*a.*(1-2.*f(x, -eps)).*a.*(1-2.*f(f(x, -eps), +eps))) < 1)...
   +(abs(a*(1-2.*x).*a.*(1-2.*f(x, +eps)).*a.*(1-2.*f(f(x, +eps), -eps))) < 1)...
   +(abs(a*(1-2.*x).*a.*(1-2.*f(x, -eps)).*a.*(1-2.*f(f(x, -eps), -eps))) < 1)...
   +(abs(a*(1-2.*x).*a.*(1-2.*f(x, +eps)).*a.*(1-2.*f(f(x, +eps), +eps))) < 1));
aux2 = ((abs(a*(1-2.*x).*a.*(1-2.*f(x, +eps))) < 1)...
   +(abs(a*(1-2.*x).*a.*(1-2.*f(x, -eps)))<1));
% aux1 = (abs(a*(1-2.*x)) < 1);
aux = aux3 + aux2;
IN_exp_two_three_sure = (aux < 1);
end

% Contracting: Complement to Surely Expanding on third and second iterates
function IN_con_compl_two = IN_con_compl_two(x, eps)
a = 3.83;
f = @(x,omega)(a*x(:,1).*(1-x(:,1)) + omega);
aux3 = ((abs(a*(1-2.*x).*a.*(1-2.*f(x, -eps)).*a.*(1-2.*f(f(x, -eps), +eps))) < 1)...
   +(abs(a*(1-2.*x).*a.*(1-2.*f(x, +eps)).*a.*(1-2.*f(f(x, +eps), -eps))) < 1)...
   +(abs(a*(1-2.*x).*a.*(1-2.*f(x, -eps)).*a.*(1-2.*f(f(x, -eps), -eps))) < 1)...
   +(abs(a*(1-2.*x).*a.*(1-2.*f(x, +eps)).*a.*(1-2.*f(f(x, +eps), +eps))) < 1));
aux2 = ((abs(a*(1-2.*x).*a.*(1-2.*f(x, +eps))) < 1)...
   +(abs(a*(1-2.*x).*a.*(1-2.*f(x, -eps)))<1));
aux = aux3 + aux2;
IN_con_compl_two = 1 - (aux < 1);
end

%Surely contracting on third or second iterate
function IN_con_two_or_three = IN_con_two_or_three(x, eps)               
a = 3.83;
f = @(x,omega)(a*x(:,1).*(1-x(:,1)) + omega);
% If one of them is expanding, then it's not contracting anymore.
aux3 = (abs(a*(1-2.*x).*a.*(1-2.*f(x, -eps)).*a.*(1-2.*f(f(x, -eps), +eps))) < 1)...
    +(abs(a*(1-2.*x).*a.*(1-2.*f(x, +eps)).*a.*(1-2.*f(f(x, +eps), -eps))) < 1)...
    +(abs(a*(1-2.*x).*a.*(1-2.*f(x, -eps)).*a.*(1-2.*f(f(x, -eps), -eps))) < 1)...
    +(abs(a*(1-2.*x).*a.*(1-2.*f(x, +eps)).*a.*(1-2.*f(f(x, +eps), +eps))) < 1);
aux2 = ((abs(a*(1-2.*x).*a.*(1-2.*f(x, +eps))) < 1)...
   +(abs(a*(1-2.*x).*a.*(1-2.*f(x, -eps)))<1));
aux = ((aux3 > 0) + (aux2>0))>0;
IN_con_two_or_three = aux;
end

function IN_exp_compl_two_or_three = IN_exp_compl_two_or_three(x, eps)               
a = 3.83;
f = @(x,omega)(a*x(:,1).*(1-x(:,1)) + omega);
% If one of them is expanding, then it's not contracting anymore.
aux3 = (abs(a*(1-2.*x).*a.*(1-2.*f(x, -eps)).*a.*(1-2.*f(f(x, -eps), +eps))) < 1)...
      +(abs(a*(1-2.*x).*a.*(1-2.*f(x, +eps)).*a.*(1-2.*f(f(x, +eps), -eps))) < 1)...
      +(abs(a*(1-2.*x).*a.*(1-2.*f(x, -eps)).*a.*(1-2.*f(f(x, -eps), -eps))) < 1)...
      +(abs(a*(1-2.*x).*a.*(1-2.*f(x, +eps)).*a.*(1-2.*f(f(x, +eps), +eps))) < 1);
aux2 = ((abs(a*(1-2.*x).*a.*(1-2.*f(x, +eps))) < 1)...
   +(abs(a*(1-2.*x).*a.*(1-2.*f(x, -eps)))<1));
aux =  1 - (((aux3 > 0) + (aux2>0))>0);
IN_exp_compl_two_or_three = aux;
end


