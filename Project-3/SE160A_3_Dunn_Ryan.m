%% Preface %%

clear all; close all; clc;
format long;
Name = {'Ryan Dunn'};
PID = {'A15600858'};

%% Input Processing %%

inFile  = 'SE160A_3_Composites_Input.xlsx';
outFile = 'SE160A_3_Composites_Output.xlsx';
insheet = readcell(inFile,'Sheet',1);   
title = insheet(7,4);
Innout = cell2mat(insheet(10:11,5));
unitsout = cell2mat(insheet(11,5));

% Fiber Properties
FName = insheet(19,5);
FProp = cell2mat(insheet(20:28,5));
Efl = insheet{20,5};    
Eft = insheet{21,5};
Gflt = insheet{22,5};
Gftt = insheet{23,5};
Vflt = insheet{24,5};
Vftt = insheet{25,5};
Fdens = insheet{26,5};
Fft = insheet{27,5};
Ffc = insheet{28,5};

% Resin Properties
RName = insheet(32,5);
RProp = cell2mat(insheet(33:39,5));
Er = insheet{33,5};
Gr = insheet{34,5};
Vr = insheet{35,5};
Rdens = insheet{36,5};
Fmt = insheet{37,5};
Fmc = insheet{38,5};
Fms = insheet{39,5};

% Micro Mechanics
MMech = insheet(44:46,5);
AFW = MMech{2};
if MMech{1}==1
    AFWUnits = 'oz/yd^2';
elseif MMech{1}==2
    AFWUnits = 'g/m^2';
    AFW = AFW*0.0294935;
end
Wres = MMech{3}*.01;

% Lamina z-angle
Lang = insheet{52,5};
% Laminate Behavior
LaminateBehav = cell2mat(insheet(59:64,5));
numplies = insheet{59,5};
tply = insheet{60,5};
Nx = insheet{61,5};
Ny = insheet{62,5};
Nxy = insheet{63,5};
SF = insheet{64,5};
ldefn = cell2mat(insheet((72-numplies):71,9));

%% 1. Weight Properties
Wfib = 1-Wres;
cdens = 1/((Wfib/Fdens)+(Wres/Rdens));
Vff = Wfib*cdens/Fdens;
Vfr = 1-Vff;
tplycured = AFW*(1/(3^4*2^8))/(Fdens*Vff); % units in inches

%% 2. Micro Mechanics 
El = Efl*Vff+Er*Vfr;
Et = Eft*Er/(Eft*Vfr+Er*Vff);
Glt = Gflt*Gr/(Gflt*Vfr+Gr*Vff);
Vlt = Vflt*Vff+Vr*Vfr;
Vtt = Vftt*Vff+Vr*Vfr;
Gtt = Et/(2*Vtt+2);
F1t = Fft*(Vff+Er*Vr/Efl);
F2t = Fmt*(1-(Vff^(1/2)-Vff)*(1-Er/Eft));
F1ccomp = Vff*Ffc;
F1cbuck = -Gr/(1-Vff*(1-Gr/Gflt))*10^3;
F1cdelam = -12.5*Fms;
F1c = max([F1ccomp, F1cbuck, F1cdelam]);
F2c = Fmc*(1-(Vff^(1/2)-Vff)*(1-Er/Eft));
Fsh = Fms*(1-(Vff^(1/2)-Vff)*(1-Gr/Gflt));

%% 3. Lamina Behavior
Vtl = Vlt*(Et/El);
Q11 = El/(1-Vlt*Vtl);
Q12 = Vlt*Et/(1-Vlt*Vtl);
Q22 = Et/(1-Vlt*Vtl);
Q66 = Glt;
Qmat = [Q11 Q12 0; Q12 Q22 0; 0 0 Q66];
Smat = inv(Qmat);
T1 = @(d) [cosd(d)^2 sind(d)^2 2*sind(d)*cosd(d);...
      sind(d)^2 cosd(d)^2 -2*sind(d)*cosd(d);...
      -sind(d)*cosd(d) sind(d)*cosd(d) cosd(d)^2-sind(d)^2];
T2 = @(d) [cosd(d)^2 sind(d)^2 -2*sind(d)*cosd(d);...
      sind(d)^2 cosd(d)^2 2*sind(d)*cosd(d);...
      sind(d)*cosd(d) -sind(d)*cosd(d) cosd(d)^2-sind(d)^2];
Qbar = T2(Lang)*Qmat*(T2(Lang)');
Sbar = (T1(Lang)')*Smat*T1(Lang);

Exa = 1/Sbar(1,1);
Eya = 1/Sbar(2,2);
Gxya = 1/Sbar(3,3);
Vxya = -Sbar(1,2)*Exa;
nxxya = Sbar(1,3)*Exa;
nyxya = Sbar(2,3)*Eya;
Fstar = [F1t F1c F2t F2c Fsh]'/SF;

% % Question 2 in writeup
% results = zeros(19,2);
% results(:,1) = [0:5:90];
% index=1;
% for Testang = 0:5:90
%     Sbar = (T1(Testang)')*Smat*T1(Testang);
%     nxxya = Sbar(1,3)*Exa;
%     
%     results(index,2) = nxxya;
%     index=index+1;
% end
% results
% return

%% 4. Laminate Behavior
Amat = zeros(3,3);
for n=1:numplies
    Amat=Amat+(T2(ldefn(n))*Qmat*(T2(ldefn(n))'))*tply;
end
Amat=Amat*10^6;               % M-lb/in to lb/in
Astar=inv(Amat);
Ls = Astar*[Nx Ny Nxy]'; % inch/inch
tlam = tply*numplies;
Ext = 1/(tlam*Astar(1,1))/10^6; % lb/in to Msi
Eyt = 1/(tlam*Astar(2,2))/10^6; % lb/in to Msi
Gxyt = 1/(tlam*Astar(3,3))/10^6; % lb/in to Msi
Vxyt = -Astar(1,2)/Astar(1,1);
nxxyt = Astar(1,3)/Astar(1,1);
nyxyt = Astar(1,3)/Astar(2,2);
isocheck = (Ext+Eyt)/(Gxyt*(4+4*Vxyt));

%% Chart work
Top = [tlam/2:-tply:(tlam/2-(numplies-1)*tply)]';
Bot = [(tlam/2-tply):-tply:(tlam/2-numplies*tply)]';
for n=1:numplies
    Qlam{n} = T2(ldefn(n))*Qmat*(T2(ldefn(n))');
    sxy{n} = (Qlam{n}*Ls)'*10^3;                    % Qlam is in Msi, *e3 for Ksi
    s12{n} = (T1(ldefn(n))*sxy{n}')';
end

%% Echo & Units %%
% Generate Units
Strength = 'Ksi';
Modulus = 'Msi';
Density = 'lb/in^3';
NLoad = 'lb/inch';
sunit = '1/Msi';
tunit = 'inch';
Stiffness = 'inch/inch';
echo_f1 = {Modulus Modulus Modulus Modulus}';
echo_f2 = {Density Strength Strength}';
echo_r1 = {Modulus Modulus}';
echo_r2 = {Density Strength Strength Strength}';
echo_a = {AFWUnits}';
echo_l = {tunit NLoad NLoad NLoad}';
uout_wp = {Density tunit}'; 
uout_mmech1 = {Modulus Modulus Modulus Modulus}';
uout_mmech2 = {Strength Strength Strength Strength Strength Strength Strength Strength}';
uout_lam1 = {Modulus Modulus Modulus}';
uout_lam2 = {Strength Strength Strength Strength Strength}';
uout_lamnt = {tunit Modulus Modulus Modulus}';

% Echo
xlswrite(outFile,Name,1,'D7');
xlswrite(outFile,PID,1,'D8');
xlswrite(outFile,title,1,'D10');
xlswrite(outFile,FName,1,'E20');
xlswrite(outFile,RName,1,'E33');
xlswrite(outFile,Innout,1,'E15');
xlswrite(outFile,FProp,1,'E21');
xlswrite(outFile,RProp,1,'E34');
xlswrite(outFile,MMech,1,'E44');
xlswrite(outFile,Lang,1,'E50');  
xlswrite(outFile,LaminateBehav,1,'E54');
loutput='I%d';
loutput=sprintf(loutput,(60-numplies));
xlswrite(outFile,ldefn,1,loutput);
 
% Echo Units
xlswrite(outFile,echo_f1,1,'F21');
xlswrite(outFile,echo_f2,1,'F27');
xlswrite(outFile,echo_r1,1,'F34');
xlswrite(outFile,echo_r2,1,'F37');
xlswrite(outFile,echo_a,1,'F45');
xlswrite(outFile,echo_l,1,'F55');

% Output Units
xlswrite(outFile,uout_wp,1,'F72');
xlswrite(outFile,uout_mmech1,1,'F78');
xlswrite(outFile,uout_mmech2,1,'F84');
xlswrite(outFile,{Modulus},1,'H96');
xlswrite(outFile,{sunit},1,'H100');
xlswrite(outFile,{Modulus},1,'H114');
xlswrite(outFile,{sunit},1,'H118');
xlswrite(outFile,uout_lam1,1,'F123');
xlswrite(outFile,uout_lam2,1,'F130');
xlswrite(outFile,{NLoad},1,'I141');
xlswrite(outFile,{Stiffness},1,'K148');
xlswrite(outFile,uout_lamnt,1,'F154');
xlswrite(outFile,{tunit},1,'F164');
xlswrite(outFile,{Strength},1,'I164');
xlswrite(outFile,{Strength},1,'L164');

%% Write Answers %%

% Generate Tables
Tab1 = [Vff Vfr Wfib Wres cdens tplycured]';
Tab2 = [El Et Glt Gtt Vlt Vtt  F1t F2t F1c F1ccomp F1cbuck F1cdelam F2c Fsh]';
Tab3 = [Exa Eya Gxya Vxya nxxya nyxya]';
Tab4 = [tlam Ext Eyt Gxyt Vxyt nxxyt nyxyt isocheck]';
for n=1:numplies
    Tab5(n,:) = [ldefn(n) Bot(n) Top(n) sxy{n} s12{n}];
end
% Section 1
xlswrite(outFile,Tab1,1,'E68');
% Section 2
xlswrite(outFile,Tab2,1,'E78');
% Section 3
xlswrite(outFile,Qmat,1,'E95');
xlswrite(outFile,Smat,1,'E99');
xlswrite(outFile,Lang,1,'E103');
xlswrite(outFile,T1(Lang),1,'E105');
xlswrite(outFile,T2(Lang),1,'E109');
xlswrite(outFile,Qbar,1,'E113');
xlswrite(outFile,Sbar,1,'E117');
xlswrite(outFile,Tab3,1,'E123');
xlswrite(outFile,Fstar,1,'E130');
% Section 4
xlswrite(outFile,Amat,1,'E140');
xlswrite(outFile,Astar,1,'E147');
xlswrite(outFile,Ls,1,'J147');
xlswrite(outFile,Tab4,1,'E154');
toutput='D%d';
toutput=sprintf(toutput,(178-numplies));
xlswrite(outFile,Tab5,1,toutput);