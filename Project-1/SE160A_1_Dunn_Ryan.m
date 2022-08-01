%% Preface %%

clear all; close all; clc;
format long;
Name = {'Ryan Dunn'};
PID = {'A15600858'};

%% Input Processing %%

inFile  = 'SE160A_1_AirLoads_Input.xlsx';
outFile = 'SE160A_1_AirLoads_Output.xlsx';
insheet = readcell(inFile,'Sheet',1);
title = insheet(7,4);
InnOut = insheet(10:11,5);
unitsin = 1;
unitsout = 1;

AeroDefn = insheet(16:29,5);
dclda = insheet{16,5};
a0  = insheet{17,5};
clmaxpos = insheet{18,5};
stallangpos = insheet{19,5};
clmaxneg = insheet{20,5};
stallangneg = insheet{21,5};
Cm0 = insheet{22,5};
Cd0fuse = insheet{23,5};
Cd0wing = insheet{24,5};
Cdi = insheet{25,5};
d0 = insheet{26,5};
d10 = insheet{27,5};
Pmax = insheet{28,5};
neff = insheet{29,5};

FStowd = insheet(36:45,5);
Vstallstow = insheet{36,5};
Vcruise = insheet{37,5};
Vdive = insheet{38,5};
Nmaxstow = insheet{39,5};
Nminstow = insheet{40,5};
gustupcruise = insheet{41,5};
gustupdive = insheet{42,5};
gustdowncruise = insheet{43,5};
gustdowndive = insheet{44,5};
Kg = insheet{45,5};

FExtended = insheet(49:52,5);
Vstallext = insheet{49,5};
Vmaxext = insheet{50,5};
Nmaxext = insheet{51,5};
Nminext = insheet{52,5};

Mainwingdefn = insheet(58:69,5);
MWweight = insheet{58,5};
MWspan = insheet{59,5};
MWfusespan = 10;
MWailspan = 45;
MWailchord = 20;
MWflapspan = 45;
MWflapchord = 20;
MWcroot = insheet{65,5};
MWctip = insheet{66,5};
MWincangroot = 0;
MWincangtip = 0;
MWqchord = insheet{69,5};

Horizstabdefn = insheet(73:76,5);
Hspan = 50;
Hcroot = 40;
Hctip = 40;
Hqchord = 220;

Vertstabdefn = insheet(80:83,5);
Vspan = insheet{80,5};
Vcroot = insheet{81,5};
Vctip = insheet{82,5};
Vqchord = insheet{83,5};

Stations = insheet(90:96,3);
airstat = insheet{90,3};
fuelstat = insheet{91,3};
pilotstat = insheet{92,3};
copilotstat = insheet{93,3};
passstat = insheet{94,3};
lugcstat = insheet{95,3};
lugfstat = insheet{96,3};

Weights = insheet(90:96,5);
airw = insheet{90,5};
fuelw = insheet{91,5};
pilotw = insheet{92,5};
copilotw = insheet{93,5};
passw = insheet{94,5};
lugcw = insheet{95,5};
lugfw = insheet{96,5};

config{1}=insheet(90:96,6);
config{2}=insheet(90:96,7);
config{3}=insheet(90:96,8);
config{4}=insheet(90:96,9);
config{5}=insheet(90:96,10);
config{6}=insheet(90:96,11);
config{7}=insheet(90:96,12);
config{8}=insheet(90:96,13);

Alt = insheet{101,5};
Tempa = insheet{102,5};

%% Calculations %%

index = cell(8,1);
for n=1:8
    index{n} = zeros(7,1);
end
for n=1:8
    x=1;
    for m=1:7
        if ismissing(config{n}{m}) ~= 1
            index{n}(x)=m;
            x=x+1;
        end
    end
end

netweight = zeros(1,8);
CG = zeros(1,8);
for n=1:8
    for m=1:7
        switch index{n}(m)
            case 1
                netweight(n) = netweight(n)+airw;
                CG(n) = CG(n)+airw*airstat;
            case 2
                netweight(n) = netweight(n)+fuelw;
                CG(n) = CG(n)+fuelw*fuelstat;
            case 3
                netweight(n) = netweight(n)+pilotw;
                CG(n) = CG(n)+pilotw*pilotstat;
            case 4
                netweight(n) = netweight(n)+copilotw;
                CG(n) = CG(n)+copilotw*copilotstat;
            case 5
                netweight(n) = netweight(n)+passw;
                CG(n) = CG(n)+passw*passstat;
            case 6
                netweight(n) = netweight(n)+lugcw;
                CG(n) = CG(n)+lugcw*lugcstat;
            case 7
                netweight(n) = netweight(n)+lugfw;
                CG(n) = CG(n)+lugfw*lugfstat;
        end
    end
end
CG = CG./netweight;

figure(1)
hold on
plot(CG(1),netweight(1),'m+','MarkerSize',10);
plot(CG(2),netweight(2),'m*','MarkerSize',10);
plot(CG(3),netweight(3),'c.','MarkerSize',10);
plot(CG(4),netweight(4),'r*','MarkerSize',10);
plot(CG(5),netweight(5),'gd','MarkerSize',10);
plot(CG(6),netweight(6),'bx','MarkerSize',10);
plot(CG(7),netweight(7),'ks','MarkerSize',10);
plot(CG(8),netweight(8),'k+','MarkerSize',10);
legend('1','2','3','4','5','6','7','8')
xlabel('CG location (in)')
ylabel('Net Weights (lb)')



%% V-n Diagram

% Atmosphere Calculations
R = 1716.5488;
TempaR = Tempa + 459.67;
T0stp = 518.67;
p0stp = 2116.224;
dens0 = 0.00237691;
Tastp = T0stp-(.00356616)*Alt;
Pastp = p0stp*(Tastp/T0stp)^5.255912;
densa = Pastp/(R*TempaR);
Mach = sqrt(1.4*R*TempaR);
Valtitude = @(vsea) vsea*(dedensns0/densa);
Powaltitude = @(powsea) powsea*(densa/dens0 - (1 - densa/dens0)/7.55)^(-1);
% N-loadings
maxweight = sum([Weights{:}]);
meanc = (MWcroot+MWctip)/2;
lambda = MWctip/MWcroot;
MWarea = (MWspan*meanc)*12*12; % ft^2

% Velocities (mph)

Vflutter = 1.2*Vdive*15/22;
VPHAA = sqrt(Nmaxstow)*Vstallstow*15/22;
VPLAA = Vdive*15/22;
VNHAA = sqrt(-Nminstow)*Vstallstow*15/22;
VNLAA = Vdive*15/22;
Vgcpos = Vcruise*15/22;
Vgdpos = Vdive*15/22;
Vgcneg = Vcruise*15/22;
Vgdneg = Vdive*15/22;

Vmat0 = [Vstallext*15/22 Vstallstow*15/22 Vcruise Vdive Vflutter VPHAA VPLAA VNHAA VNLAA...
         Vgcpos Vgdpos Vgcneg Vgdneg]';

Ngposcruise = 1 + Kg*(0.5*dens0*gustupcruise*(MWarea/maxweight)*dclda*...
                  (1/360))*Vcruise;
Ngposdive = 1 + Kg*(0.5*dens0*gustupdive*(MWarea/maxweight)*dclda*...
                (1/360))*Vdive;
Ngnegcruise = 1 + Kg*(0.5*dens0*gustdowncruise*(MWarea/maxweight)*dclda*...
                  (1/360))*Vcruise;
Ngnegdive = 1 + Kg*(0.5*dens0*gustdowndive*(MWarea/maxweight)*dclda*...
                (1/360))*Vdive;
Nmat0 = [1 1 1 1 1 Nmaxstow Nmaxstow Nminstow Nminstow Ngposcruise Ngposdive...
         Ngnegcruise Ngnegdive]';

% Velocity & N-loads @ altitude
Vmatalt = Vmat0*sqrt(dens0/densa);
Ngpcalt = 1 + Kg*(0.5*densa*gustupcruise*(MWarea/maxweight)*dclda*...
                  (1/360))*Vcruise;
Ngpdalt = 1 + Kg*(0.5*densa*gustupdive*(MWarea/maxweight)*dclda*...
                (1/360))*Vdive;
Ngncalt = 1 + Kg*(0.5*densa*gustdowncruise*(MWarea/maxweight)*dclda*...
                  (1/360))*Vcruise;
Ngndalt = 1 + Kg*(0.5*densa*gustdowndive*(MWarea/maxweight)*dclda*...
                (1/360))*Vdive;
Nmatalt = [1 1 1 1 1 Nmaxstow Nmaxstow Nminstow Nminstow Ngpcalt Ngpdalt...
         Ngncalt Ngndalt]';

%% V-n Diagram



%% Normal flight

Valtn = Vmatalt(1:4)';
qaltn = 0.5*densa*(Valtn*(22/15)).^2;
nloadan = [1 1 1 1];

Valtg = [Vmatalt(6:9)' Vmat0(10:13)'];
qaltg = 0.5*densa*(Valtg*22/15).^2;
nloadag = [Nmaxstow Nmaxstow Nminstow Nminstow Nmat0(10:13)'];

%% Echo %%

createFigure(outFile,1,1,'C127','M141')
xlswrite(outFile,Name,1,'D7');
xlswrite(outFile,PID,1,'D8');
xlswrite(outFile,title,1,'D10');
xlswrite(outFile,AeroDefn,1,'E21');
xlswrite(outFile,FStowd,1,'E41');
xlswrite(outFile,FExtended,1,'E54');
xlswrite(outFile,Mainwingdefn,1,'E63');
xlswrite(outFile,Horizstabdefn,1,'E78');
xlswrite(outFile,Vertstabdefn,1,'E85');
xlswrite(outFile,Stations,1,'C95');
xlswrite(outFile,Weights,1,'E95');
xlswrite(outFile,Alt,1,'E106');
xlswrite(outFile,Tempa,1,'E107');


%% Write Answers %%

xlswrite(outFile,Stations,1,'C116');
xlswrite(outFile,Weights,1,'E116');
xlswrite(outFile,netweight,1,'F124');
xlswrite(outFile,CG,1,'F125');
xlswrite(outFile,densa,1,'E147');
xlswrite(outFile,Vmat0,1,'E151');
xlswrite(outFile,Nmat0,1,'F151');
xlswrite(outFile,Vmatalt,1,'G151');
xlswrite(outFile,Nmatalt,1,'H151');
xlswrite(outFile,Valtn,1,'E193');
xlswrite(outFile,qaltn,1,'E194');
xlswrite(outFile,nloadan,1,'E195');
xlswrite(outFile,Valtg,1,'E214');
xlswrite(outFile,qaltg,1,'E215');
xlswrite(outFile,nloadag,1,'E216');