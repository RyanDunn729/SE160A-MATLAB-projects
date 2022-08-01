%% Preface %%

clear all; close all; clc;
format long;
Name = {'Ryan Dunn'};
PID = {'A15600858'};

%% Input Processing %%

inFile  = 'SE160A_5_Wing_Bending_Input.xlsx';
outFile = 'SE160A_5_Wing_Bending_Output.xlsx';
insheet = readcell(inFile,'Sheet',1);
title = insheet(7,4);
unitsout = insheet{11,5};
    
% Stringers
Yst = cell2mat(insheet(15,5:8));
Zst = cell2mat(insheet(16,5:8));
Ast = cell2mat(insheet(17,5:8));
Iyyst = cell2mat(insheet(18,5:8));
Izzst = cell2mat(insheet(19,5:8));
Iyzst = cell2mat(insheet(20,5:8));
Est = cell2mat(insheet(21,5:8))*10^6;
Ytst = cell2mat(insheet(22,5:8));
Utst = cell2mat(insheet(23,5:8));
Ycst = cell2mat(insheet(24,5:8));
Ucst = cell2mat(insheet(25,5:8));

% Skins
tsk = cell2mat(insheet(29,5:8));
Gsk = cell2mat(insheet(30,5:8));
Tysk = cell2mat(insheet(31,5:8));
Tusk = cell2mat(insheet(32,5:8));

Winggeo = cell2mat(insheet(36:41,5));
L = cell2mat(insheet(36,5));
c = cell2mat(insheet(37,5));
Wwing = cell2mat(insheet(38,5));
nload = cell2mat(insheet(39,5));
SFy = cell2mat(insheet(40,5));
SFu = cell2mat(insheet(41,5));

Wingaero = cell2mat(insheet(45:51,5));
z0 = cell2mat(insheet(45,5));
z2 = cell2mat(insheet(46,5));
z4 = cell2mat(insheet(47,5));
y0 = cell2mat(insheet(48,5));
ynth = cell2mat(insheet(49,5));
nth = cell2mat(insheet(50,5));
m0 = cell2mat(insheet(51,5));

%% Wing Properties

EYc = sum([Est.*Yst.*Ast])/sum([Est.*Ast]);
EZc = sum([Est.*Zst.*Ast])/sum([Est.*Ast]);

% Stringers w/ Shift
Iyystleading = Iyyst+Ast.*((Zst-EZc).^2);
Izzstleading = Izzst+Ast.*((Yst-EYc).^2);
Iyzstleading = Iyzst+Ast.*((Yst-EYc).*(Zst-EZc));
EAst = Est.*Ast;
EA = sum(EAst);
EIyy = sum(Est.*Iyystleading);
EIzz = sum(Est.*Izzstleading);
EIyz = sum(Est.*Iyzstleading);

% Skin (straight segments)
if Yst(2)>=c/2
    xpts = [Yst(2) Yst(3) Yst(4)];
    ypts = [Zst(2) Zst(3) Zst(4)];
    for n=1:2
        x1=xpts(n);
        x2=xpts(n+1);
        y1=ypts(n);
        y2=ypts(n+1);
        Lseg(n) = ((x2-x1)^2+(y2-y1)^2)^(1/2);
    end
    L23 = Lseg(1);
    L34 = Lseg(2);
    Asq = 0;
    Atr = 0.5*(Zst(4)-Zst(2))*(c-Yst(2));
elseif Yst(2)<(c/2)
    xpts = [Yst(2) c/2 Yst(3) c/2 Yst(4)];
    ypts = [Zst(2) Zst(2) Zst(3) Zst(4) Zst(4)];
    for n=1:4
        x1=xpts(n);
        x2=xpts(n+1);
        y1=ypts(n);
        y2=ypts(n+1);
        Lseg(n) = ((x2-x1)^2+(y2-y1)^2)^(1/2);
    end
    L23 = Lseg(1)+Lseg(2);
    L34 = Lseg(3)+Lseg(4);
    Asq = (Zst(4)-Zst(3))*(c/2-Yst(2));
    Atr = 0.5*(Zst(4)-Zst(2))*(c-Yst(2));
end

% Skin (ellipse)
a = Yst(2);
b = Zst(4);
Lell = pi*(3*(a+b)-((3*a+b)*(a+3*b))^(1/2));
Aell = pi*a*b;
Area = Asq+Atr+Aell/2;
L12 = Lell/4;
L41 = Lell/4;
skin = [L12 L23 L34 L41];
perim = sum(skin);

GJ = 4*(Area^2)/sum([skin./(Gsk.*tsk)])*10^6;

%% Stress Resultants @ Root

Lift = @(x) z0+z2*(x/L).^2+z4*(x/L).^4-nload*Wwing;
Drag = @(x) y0+ynth*(x/L).^(nth);

MxEqn = @(x) (c/4-EYc)*(z0+z2*(x/L).^2+z4*(x/L).^4)+m0-(c/2-EYc)*nload*Wwing;
MyEqn = @(x) -x.*(z0+z2*(x/L).^2+z4*(x/L).^4-nload*Wwing);
MzEqn = @(x) x.*(y0+ynth*(x/L).^(nth));

Px = 0;
Py = integral(Drag,0,L);
Pz = integral(Lift,0,L);
Mx = integral(MxEqn,0,L);
My = integral(MyEqn,0,L);
Mz = integral(MzEqn,0,L);

%% Stringer Stress Analysis

Sxx = (Px+(My*Est.*(Zst-EZc)./EIyy)-(Mz*Est.*(Yst-EYc)./EIzz))*10^-3;
Ststar = [min([Utst(1)/SFu,Ytst(1)/SFy]) min([Utst(2)/SFu,Ytst(2)/SFy])...
          min([Utst(2)/SFu,Ytst(3)/SFy]) min([Utst(4)/SFu,Ytst(4)/SFy])];
Scstar = [max([Ucst(1)/SFu,Ycst(1)/SFy]) max([Ucst(2)/SFu,Ycst(2)/SFy])...
          max([Ucst(3)/SFu,Ycst(3)/SFy]) max([Ucst(4)/SFu,Ycst(4)/SFy])];
for n=1:4
    if Sxx(n)<0
        MSxx(n) = Scstar(n)/Sxx(n)-1;
    elseif Sxx(n)>0
        MSxx(n) = Ststar(n)/Sxx(n)-1;
    else
        MSxx(n) = Inf;
    end
end

%% Wing Skin Stress Analysis

A12 = Aell/4;
A14 = Aell/4;

% if Yst(2)<(c/2)
%     A23 = Zst(4)*(c/2-Yst(2)) + 0.5*Zst(4)*c/2;
%     A34 = A23;
% elseif Yst(2)>=(c/2)
    A23 = 0.5*Zst(4)*(c-Yst(2));
    A34 = A23;
% end
A1234 = [A12 A23 A34 A14];
Mxshift = Pz*(EYc-Yst(2))+Mx;

Mat = [1 -1 0 0;
       0 1 -1 0;
       0 0 1 -1;
       2*A12 2*A23 2*A34 2*A14];
Vymat = (Py/EIzz*[EAst(2)*(Yst(2)-EYc) EAst(3)*(Yst(3)-EYc) EAst(4)*(Yst(4)-EYc) 0])';
Vzmat = (Pz/EIyy*[EAst(2)*Zst(2) EAst(3)*Zst(3) EAst(4)*Zst(4) 0])';
Mmat = [0 0 0 Mxshift]';
qflow = (inv(Mat)*(Vymat+Vzmat+Mmat))';
Tsk = qflow./(tsk)*10^-3;
Tstar = [min([Tusk(1)/SFu,Tysk(1)/SFy]), min([Tusk(2)/SFu,Tysk(2)/SFy]),...
         min([Tusk(3)/SFu,Tysk(3)/SFy]) min([Tusk(4)/SFu,Tysk(4)/SFy])];
MSxy = Tstar./abs(Tsk)-1;

%% Shear centers

% Define unit shear flows
Vec = (1/EIyy*[EAst(2)*(Zst(2)-EZc) EAst(3)*(Zst(3)-EZc) EAst(4)*(Zst(4)-EZc) 0])';
q1 = (inv(Mat)*Vec)';
Vec = (1/EIzz*[EAst(2)*(Yst(2)-EYc) EAst(3)*(Yst(3)-EYc) EAst(4)*(Yst(4)-EYc) 0])';
q2 = (inv(Mat)*Vec)';
Vec = [0 0 0 1]';
q3 = (inv(Mat)*Vec)';
ymotion = sum(q1.*skin./(Gsk.*tsk));
zmotion = sum(q2.*skin./(Gsk.*tsk));
momentmotion = sum(q3.*skin./(Gsk.*tsk));

eycenter = Yst(2)-ymotion/momentmotion;
ezcenter = 0;

% Wing Tip Displacement / Rotation

wtip = 5

vtip = 5

Momentatpoint = @(x) (c/4-eycenter)*(z0*x+z2/3*((x.^3)/(L^2))+z4/5*((x.^5)/(L^4)))+m0*x-(c/2-eycenter)*nload*Wwing*x;
phi = integral(Momentatpoint,0,L)/GJ*180/pi

%% Echo

Table1 = [Yst; Zst; Ast; Iyyst; Izzst; Iyzst; Est*10^-6; Ytst; Utst; Ycst; Ucst]; 
Table2 = [tsk; Gsk; Tysk; Tusk];
Table3 = [EYc EZc EA EIyy EIzz EIyz GJ eycenter ezcenter]';
Table4 = [Px Py Pz Mx My Mz]';
Table5 = [Sxx; Ststar; Scstar; MSxx];
Table6 = [Tsk; Tstar; MSxy];

xlswrite(outFile,Name,1,'D7');
xlswrite(outFile,PID,1,'D8');
xlswrite(outFile,title,1,'D10');
xlswrite(outFile,[1,1]',1,'E15');
xlswrite(outFile,Table1,1,'E20');
xlswrite(outFile,Table2,1,'E34');
xlswrite(outFile,Winggeo,1,'E41');
xlswrite(outFile,Wingaero,1,'E50');
xlswrite(outFile,Table3,1,'E63');
xlswrite(outFile,Table4,1,'E75');
xlswrite(outFile,Table5,1,'E84');
xlswrite(outFile,Table6,1,'E91');
xlswrite(outFile,phi,1,'E99');