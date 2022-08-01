%% Preface %%

clear all; close all; clc;
format long;
Name = {'Ryan Dunn'};
PID = {'A15600858'};

%% Input Processing %%

inFile  = 'SE160A_4_Section_Input.xlsx';
outFile = 'SE160A_4_Section_Output.xlsx';
deleteFigure(outFile);
insheet = readcell(inFile,'Sheet',1);
title = insheet(7,4);
unitsin = insheet{10,5};
unitsout = insheet{11,5};

isection = insheet{14,5};
x0 = insheet{15,5};
y0 = insheet{16,5};

% Generate Units Echo
LUnits = 'inch';
AUnits = 'inch^2';
IUnits = 'inch^4';
MwUnits = 'lb/inch';
MiUnits = 'lb-inch';
DensUnits = 'lb/in^3';
ModUnits = 'Msi';

% Generate Cell/Tables
Units1 = {AUnits LUnits LUnits IUnits IUnits IUnits}';
Units2 = {MwUnits LUnits LUnits MiUnits MiUnits MiUnits}';
Units3 = {LUnits LUnits DensUnits ModUnits}';

if isection==1 % General case
    %% General Case (isection=1)
    NNodes = insheet{22,5};
    NSegments = insheet{23,5};
    for n=1:NNodes
        Xn(n) = insheet{26+n,4};
        Yn(n) = insheet{26+n,5};
        An(n) = insheet{26+n,6};
        Ixxn(n) = insheet{26+n,7};
        Iyyn(n) = insheet{26+n,8};
        Ixyn(n) = insheet{26+n,9};
        densn(n) = insheet{26+n,10};
        En(n) = insheet{26+n,11};
    end
    for n=1:NSegments
        seg1(n) = insheet{49+n,4};
        seg2(n) = insheet{49+n,5};
        tseg(n) = insheet{49+n,6};
        densseg(n) = insheet{49+n,7};
        Eseg(n) = insheet{49+n,8};
    end

    %% Setup Segments 
    for n=1:length(seg1)
        x1g(n) = Xn(seg1(n));
        x2g(n) = Xn(seg2(n));
        y1g(n) = Yn(seg1(n));
        y2g(n) = Yn(seg2(n));
        L(n) = ((x2g(n)-x1g(n))^2+(y2g(n)-y1g(n))^2)^(1/2);
        Aseg(n) = L(n)*tseg(n);
        Xcseg(n) = (x1g(n)+x2g(n))/2;
        Ycseg(n) = (y1g(n)+y2g(n))/2;
        Ixxseg(n) = Aseg(n)*(y1g(n)*(2*y1g(n)+y2g(n))+y2g(n)*(y1g(n)+2*y2g(n)))/6;
        Iyyseg(n) = Aseg(n)*(x1g(n)*(2*x1g(n)+x2g(n))+x2g(n)*(x1g(n)+2*x2g(n)))/6;
        Ixyseg(n) = Aseg(n)*(x1g(n)*(2*y1g(n)+y2g(n))+x2g(n)*(y1g(n)+2*y2g(n)))/6;

        Ixxseg(n) = Ixxseg(n)-Aseg(n)*Ycseg(n)^2; % Shift Back
        Iyyseg(n) = Iyyseg(n)-Aseg(n)*Xcseg(n)^2; % Shift Back
        Ixyseg(n) = Ixyseg(n)-Aseg(n)*Xcseg(n)*Ycseg(n); % Shift Back
    end

    %% General Properties
    tempx=Xn;
    tempy=Yn;
    tempxs=Xcseg;
    tempys=Ycseg;
    tempIxxn=Ixxn;
    tempIyyn=Iyyn;
    tempIxyn=Ixyn;
    tempIxxseg = Ixxseg;
    tempIyyseg = Iyyseg;
    tempIxyseg = Ixyseg;
    Agen = sum([An,Aseg]);
    wAgen = sum([densn.*An,Aseg.*densseg]);
    EAgen = sum([En.*An,Aseg.*Eseg]);
    
    Xcgen = cell(1,4);
    Ycgen = cell(1,4);
    Ixxgen = cell(1,4);
    Iyygen = cell(1,4);
    Ixygen = cell(1,4);
    wXcgen = cell(1,4);
    wYcgen = cell(1,4);
    wIxxgen = cell(1,4);
    wIyygen = cell(1,4);
    wIxygen = cell(1,4);
    EXcgen = cell(1,4);
    EYcgen = cell(1,4);
    EIxxgen = cell(1,4);
    EIyygen = cell(1,4);
    EIxygen = cell(1,4);
    % GENERAL BOI
    for n=1:4
        % Shift Axes according to loop 
        % 1 = original axes
        % 2 = centroid
        % 3 = mod-weighted centroid
        % 4 = user origin
        switch n
            case 1
                Xn=tempx;
                Yn=tempy;
                Xcseg=tempxs;
                Ycseg=tempys;
            case 2
                Xn=tempx-Xcgen{1};
                Yn=tempy-Ycgen{1};
                Xcseg=tempxs-Xcgen{1};
                Ycseg=tempys-Ycgen{1};
            case 3
                Xn=tempx-EXcgen{1};
                Yn=tempy-EYcgen{1};
                Xcseg=tempxs-EXcgen{1};
                Ycseg=tempys-EYcgen{1};
            case 4
                Xn=tempx-x0;
                Yn=tempy-y0;
                Xcseg=tempxs-x0;
                Ycseg=tempys-y0;
        end
        % Define I's in new frame
        Ixxn = tempIxxn+An.*(Yn.^2);
        Iyyn = tempIyyn+An.*(Xn.^2);
        Ixyn = tempIxyn+An.*(Xn.*Yn);
        Ixxseg = (tempIxxseg+(Ycseg.^2).*Aseg);
        Iyyseg = (tempIyyseg+(Xcseg.^2).*Aseg);
        Ixyseg = (tempIxyseg+(Ycseg.*Xcseg.*Aseg));
        % Area Properties
        Xcgen{n} = sum([Xn.*An,Xcseg.*Aseg])/Agen;
        Ycgen{n} = sum([Yn.*An,Ycseg.*Aseg])/Agen;
        Ixxgen{n} = sum([Ixxn,Ixxseg]);
        Iyygen{n} = sum([Iyyn,Iyyseg]);
        Ixygen{n} = sum([Ixyn,Ixyseg]);
        % Mass Properties
        wXcgen{n} = sum([densn.*Xn.*An,densseg.*Xcseg.*Aseg])/wAgen;
        wYcgen{n} = sum([densn.*Yn.*An,densseg.*Ycseg.*Aseg])/wAgen;
        wIxxgen{n} = sum([(densn.*Ixxn),(densseg.*Ixxseg)]);
        wIyygen{n} = sum([(densn.*Iyyn),(densseg.*Iyyseg)]);
        wIxygen{n} = sum([(densn.*Ixyn),(densseg.*Ixyseg)]);
        % Stiffness Properties
        EXcgen{n} = sum([En.*Xn.*An,Eseg.*Xcseg.*Aseg])/EAgen;
        EYcgen{n} = sum([En.*Yn.*An,Eseg.*Ycseg.*Aseg])/EAgen;
        EIxxgen{n} = sum([(En.*Ixxn),(Eseg.*Ixxseg)]);
        EIyygen{n} = sum([(En.*Iyyn),(Eseg.*Iyyseg)]);
        EIxygen{n} = sum([(En.*Ixyn),(Eseg.*Ixyseg)]);
    end

    %% General Plot
    figure(1)
    hold on
    for n=1:length(seg1)
        segplot = plot([x1g(n);x2g(n)],[y1g(n);y2g(n)],'k','LineWidth',3);
    end
    for n=1:length(tempx)
        nodeplot = plot(tempx(n),tempy(n),'r.','MarkerSize',30);
    end
    axis('equal')
    legend([segplot nodeplot],'Segments','Nodes')
    xlabel('X (inches)')
    ylabel('Y (inches)')
    %% Echo to General Section (Sheet 1)
    Tablenodes = [tempx' tempy' An' tempIxxn' tempIyyn' tempIxyn' densn' En'];
    Tablesegments = [seg1' seg2' tseg' densseg' Eseg'];
    xlswrite(outFile,Name,1,'D7');
    xlswrite(outFile,PID,1,'D8');
    xlswrite(outFile,title,1,'D10');
    xlswrite(outFile,[unitsin unitsout]',1,'E15');
    xlswrite(outFile,[isection x0 y0]',1,'E19');
    xlswrite(outFile,[NNodes NSegments]',1,'E27');
    xlswrite(outFile,Tablenodes,1,'D32');
    xlswrite(outFile,Tablesegments,1,'D55');
    
    % Units
    xlswrite(outFile,Units1,1,'F83');
    xlswrite(outFile,Units2,1,'K83');
    xlswrite(outFile,Units1,1,'F92');
    xlswrite(outFile,Units1,1,'F135');
    xlswrite(outFile,Units2,1,'K135');
    xlswrite(outFile,Units1,1,'F144');
    xlswrite(outFile,Units1,1,'F155');
    xlswrite(outFile,Units2,1,'K155');
    xlswrite(outFile,Units1,1,'F164');
    xlswrite(outFile,Units1,1,'F175');
    xlswrite(outFile,Units2,1,'K175');
    xlswrite(outFile,Units1,1,'F184');
    
    %% Write Answers
    % Plots
    createFigure(outFile,1,1,'C99','L129')
    % Sheet 1
    for n=1:4
        TableAreasgen{n} = [Agen Xcgen{n} Ycgen{n} Ixxgen{n} Iyygen{n} Ixygen{n}]';
        TableMassesgen{n} = [wAgen wXcgen{n} wYcgen{n} wIxxgen{n} wIyygen{n} wIxygen{n}]';
        TableStiffsgen{n} = [EAgen EXcgen{n} EYcgen{n} EIxxgen{n} EIyygen{n} EIxygen{n}]';
    end
    xlswrite(outFile,TableAreasgen{1},1,'E83');
    xlswrite(outFile,TableAreasgen{2},1,'E135');
    xlswrite(outFile,TableAreasgen{3},1,'E155');
    xlswrite(outFile,TableAreasgen{4},1,'E175');
    xlswrite(outFile,TableMassesgen{1},1,'J83');
    xlswrite(outFile,TableMassesgen{2},1,'J135');
    xlswrite(outFile,TableMassesgen{3},1,'J155');
    xlswrite(outFile,TableMassesgen{4},1,'J175');
    xlswrite(outFile,TableStiffsgen{1},1,'E92');
    xlswrite(outFile,TableStiffsgen{2},1,'E144');
    xlswrite(outFile,TableStiffsgen{3},1,'E164');
    xlswrite(outFile,TableStiffsgen{4},1,'E184');

elseif isection==2 % NACA Case
    %% NACA Import
    NACA = cell2mat(insheet(75,5:7));
    numnacaseg = insheet{79,5};
    c = insheet{80,5};
    tnaca = insheet{81,5};
    densnaca = insheet{82,5};
    Enaca = insheet{83,5};
    for n=1:5
        if ismissing(insheet{86+n,3})
            break
        else
            xoc(n) = insheet{86+n,3};
            tspar(n) = insheet{86+n,4};
            densspar(n) = insheet{86+n,5};
            Espar(n) = insheet{86+n,6};
        end
    end
    figure(1)
    hold on
    for n=1:20
        if ismissing(insheet{94+n,4})
            break
        else
            Xst(n) = insheet{(94+n),4};
            Yst(n) = insheet{(94+n),5};
            Ast(n) = insheet{(94+n),6};
            Ixxst(n) = insheet{(94+n),7};
            Iyyst(n) = insheet{(94+n),8};
            Ixyst(n) = insheet{(94+n),9};
            densst(n) = insheet{(94+n),10};
            Est(n) = insheet{(94+n),11};
            stringerplot = plot(Xst(n),Yst(n),'b.','MarkerSize',15);  
        end
    end
    axis('equal')
    %% NACA Generation
    M = NACA(1);
    P = NACA(2);
    TT = NACA(3);
    pt=0;
    n=1;
    while pt<=(c*P/10)
        nacapts1(n)=pt;
        pt=pt+c/numnacaseg;
        n=n+1;
    end
    n=1;
    while pt>=(c*P/10) && pt<=c
        nacapts2(n)=pt;
        pt=pt+c/numnacaseg;
        n=n+1;
    end
    ysym = @(x) 5*c*TT/100*(.2969*(x/c).^(1/2)-.126*(x/c)...
        -.3516*(x/c).^2+.2843*(x/c).^3-.1015*(x/c).^4);
    if P==0
        ycam1=@(x) 0;
        ycam2=@(x) 0;
    else
        ycam1 = @(x) c*(-M/(P^2)*(x/c).^2+M/(5*P)*(x/c));
        ycam2 = @(x) c*M/(100-20*P+P^2)*(-(x/c).^2+P/5*(x/c)+1-P/5);
    end
    ynaca11 = @(x) ycam1(x)+ysym(x);
    ynaca12 = @(x) ycam2(x)+ysym(x);
    ynaca21 = @(x) ycam1(x)-ysym(x);
    ynaca22 = @(x) ycam2(x)-ysym(x);

    Ynaca1 = [ynaca11(nacapts1), ynaca12(nacapts2)];
    Ynaca2 = [ynaca21(nacapts1), ynaca22(nacapts2)];
    Xnaca = [nacapts1 nacapts2];
    %% NACA Properties
    % NACA spars
    numspar = length(xoc);
    for n=1:numspar
        Xcspar(n) = xoc(n)*c;
        if Xcspar(n)<(c*P/10)
            hspar(n) = ynaca11(Xcspar(n))-ynaca21(Xcspar(n));
            Ycspar(n) = (ynaca11(Xcspar(n))+ynaca21(Xcspar(n)))/2;
            vsparplot = plot([Xcspar(n),Xcspar(n)],[ynaca11(Xcspar(n)),ynaca21(Xcspar(n))],'r');
        elseif Xcspar(n)>=(c*P/10)
            hspar(n) = ynaca21(Xcspar(n))-ynaca22(Xcspar(n));
            Ycspar(n) = (ynaca21(Xcspar(n))+ynaca22(Xcspar(n)))/2;
            vsparplot = plot([Xcspar(n),Xcspar(n)],[ynaca12(Xcspar(n)),ynaca22(Xcspar(n))],'r');
        end
        Aspar(n) = hspar(n)*tspar(n);
        Ixxspar(n) = tspar(n)*(hspar(n)^3)/12;
        Iyyspar(n) = hspar(n)*(tspar(n)^3)/12;
        Ixyspar(n) = 0;    
    end
    
    % NACA Segments (BOT)
    for n=1:(length(Xnaca)-1)
        x1nb(n) = Xnaca(n);
        x2nb(n) = Xnaca(n+1);
        y1nb(n) = Ynaca2(n);
        y2nb(n) = Ynaca2(n+1);
        Lb(n) = ((x2nb(n)-x1nb(n))^2+(y2nb(n)-y1nb(n))^2)^(1/2);
        Asnacab(n) = Lb(n)*tnaca;
        Xcsnacab(n) = (x1nb(n)+x2nb(n))/2;
        Ycsnacab(n) = (y1nb(n)+y2nb(n))/2;
        Ixxsnacab(n) = Asnacab(n)*(y1nb(n)*(2*y1nb(n)+y2nb(n))+y2nb(n)*(y1nb(n)+2*y2nb(n)))/6;
        Iyysnacab(n) = Asnacab(n)*(x1nb(n)*(2*x1nb(n)+x2nb(n))+x2nb(n)*(x1nb(n)+2*x2nb(n)))/6;
        Ixysnacab(n) = Asnacab(n)*(x1nb(n)*(2*y1nb(n)+y2nb(n))+x2nb(n)*(y1nb(n)+2*y2nb(n)))/6;
        
        Ixxsnacab(n) = Ixxsnacab(n)-Asnacab(n)*Ycsnacab(n)^2; % Shift back
        Iyysnacab(n) = Iyysnacab(n)-Asnacab(n)*Xcsnacab(n)^2; % Shift back
        Ixysnacab(n) = Ixysnacab(n)-Asnacab(n)*Xcsnacab(n)*Ycsnacab(n); % Shift back
    end
    % NACA Segments (TOP)
    for n=1:(length(Xnaca)-1)
        x1nt(n) = Xnaca(n);
        x2nt(n) = Xnaca(n+1);
        y1nt(n) = Ynaca1(n);
        y2nt(n) = Ynaca1(n+1);
        Lt(n) = ((x2nt(n)-x1nt(n))^2+(y2nt(n)-y1nt(n))^2)^(1/2);
        Asnacat(n) = Lt(n)*tnaca;
        Xcsnacat(n) = (x1nt(n)+x2nt(n))/2;
        Ycsnacat(n) = (y1nt(n)+y2nt(n))/2;
        Ixxsnacat(n) = Asnacat(n)*(y1nt(n)*(2*y1nt(n)+y2nt(n))+y2nt(n)*(y1nt(n)+2*y2nt(n)))/6;
        Iyysnacat(n) = Asnacat(n)*(x1nt(n)*(2*x1nt(n)+x2nt(n))+x2nt(n)*(x1nt(n)+2*x2nt(n)))/6;
        Ixysnacat(n) = Asnacat(n)*(x1nt(n)*(2*y1nt(n)+y2nt(n))+x2nt(n)*(y1nt(n)+2*y2nt(n)))/6;
        
        Ixxsnacat(n) = Ixxsnacat(n)-Asnacat(n)*Ycsnacat(n)^2; % Shift back
        Iyysnacat(n) = Iyysnacat(n)-Asnacat(n)*Xcsnacat(n)^2; % Shift back
        Ixysnacat(n) = Ixysnacat(n)-Asnacat(n)*Xcsnacat(n)*Ycsnacat(n); % Shift back
    end
    x1n = [x1nb x1nt];
    x2n = [x2nb x2nt];
    y1n = [y1nb y1nt];
    y2n = [y2nb y2nt];
    L = [Lb Lt];
    Asnaca= [Asnacab Asnacat];
    Xcsnaca = [Xcsnacab Xcsnacat];
    Ycsnaca = [Ycsnacab Ycsnacat];
    Ixxsnaca = [Ixxsnacab Ixxsnacat];
    Iyysnaca = [Iyysnacab Iyysnacat];
    Ixysnaca = [Ixysnacab Ixysnacat];
    Anaca = sum([Aspar,Asnaca,Ast]);
    wAnaca = sum([Aspar.*densspar,Asnaca.*densnaca,Ast.*densst]);
    EAnaca = sum([Aspar.*Espar,Asnaca.*Enaca,Ast.*Est]);
    tXcsnaca=Xcsnaca;
    tYcsnaca=Ycsnaca;
    tXst=Xst;
    tYst=Yst;
    tXcspar=Xcspar;
    tYcspar=Ycspar;
    Xcnaca = cell(1,4);
    Ycnaca = cell(1,4);
    Ixxnaca = cell(1,4);
    Iyynaca = cell(1,4);
    Ixynaca = cell(1,4);
    wXcnaca = cell(1,4);
    wYcnaca = cell(1,4);
    wIxxnaca = cell(1,4);
    wIyynaca = cell(1,4);
    wIxynaca = cell(1,4);
    EXcnaca = cell(1,4);
    EYcnaca = cell(1,4);
    EIxxnaca = cell(1,4);
    EIyynaca = cell(1,4);
    EIxynaca = cell(1,4);
    for n=1:4
        % Shift Axes according to loop 
        % 1 = original axes
        % 2 = centroid
        % 3 = mod-weighted centroid
        % 4 = user origin
        switch n
            case 1
                Xcsnaca=tXcsnaca;
                Ycsnaca=tYcsnaca;
                Xstt=Xst;
                Ystt=Yst;
                Xcspar=tXcspar;
                Ycspar=tYcspar;
            case 2
                Xcsnaca=tXcsnaca-Xcnaca{1};
                Ycsnaca=tYcsnaca-Ycnaca{1};
                Xstt=Xst-Xcnaca{1};
                Ystt=Yst-Ycnaca{1};
                Xcspar=tXcspar-Xcnaca{1};
                Ycspar=tYcspar-Ycnaca{1};
            case 3
                Xcsnaca=tXcsnaca-EXcnaca{1};
                Ycsnaca=tYcsnaca-EYcnaca{1};
                Xstt=Xst-EXcnaca{1};
                Ystt=Yst-EYcnaca{1};
                Xcspar=tXcspar-EXcnaca{1};
                Ycspar=tYcspar-EYcnaca{1};
            case 4
                Xcsnaca=tXcsnaca-x0;
                Ycsnaca=tYcsnaca-y0;
                Xstt=Xst-x0;
                Ystt=Yst-y0;
                Xcspar=tXcspar-x0;
                Ycspar=tYcspar-y0;
        end
        % Segments
        tempIxxsnaca = Ixxsnaca+Asnaca.*(Ycsnaca.^2);
        tempIyysnaca = Iyysnaca+Asnaca.*(Xcsnaca.^2);
        tempIxysnaca = Ixysnaca+(Asnaca.*Xcsnaca.*Ycsnaca);
        % Stringers
        tempIxxst = Ixxst+Ast.*(Ystt.^2);
        tempIyyst = Iyyst+Ast.*(Xstt.^2);
        tempIxyst = Ixyst+(Ast.*Xstt.*Ystt);
        % Vertical Spars
        tempIxxspar = Ixxspar+Aspar.*(Ycspar.^2);
        tempIyyspar = Iyyspar+Aspar.*(Xcspar.^2);
        tempIxyspar = Ixyspar+(Aspar.*Xcspar.*Ycspar);        
        % Area Properties
        Xcnaca{n} = sum([Xcsnaca.*Asnaca,Xstt.*Ast,Xcspar.*Aspar])/Anaca;
        Ycnaca{n} = sum([Ycsnaca.*Asnaca,Ystt.*Ast,Ycspar.*Aspar])/Anaca;
        Ixxnaca{n} = sum([tempIxxsnaca,tempIxxst,tempIxxspar]);
        Iyynaca{n} = sum([tempIyysnaca,tempIyyst,tempIyyspar]);
        Ixynaca{n} = sum([tempIxysnaca,tempIxyst,tempIxyspar]);
        % Mass Properties
        wXcnaca{n} = sum([densnaca.*Xcsnaca.*Asnaca,densst.*Xstt.*Ast,densspar.*Xcspar.*Aspar])/wAnaca;
        wYcnaca{n} = sum([densnaca.*Ycsnaca.*Asnaca,densst.*Ystt.*Ast,densspar.*Ycspar.*Aspar])/wAnaca;
        wIxxnaca{n} = sum([densnaca.*tempIxxsnaca,densst.*tempIxxst,densspar.*tempIxxspar]);
        wIyynaca{n} = sum([densnaca.*tempIyysnaca,densst.*tempIyyst,densspar.*tempIyyspar]);
        wIxynaca{n} = sum([densnaca.*tempIxysnaca,densst.*tempIxyst,densspar.*tempIxyspar]);
        % Stiffness Properties
        EXcnaca{n} = sum([Enaca.*Xcsnaca.*Asnaca,Est.*Xstt.*Ast,Espar.*Xcspar.*Aspar])/EAnaca;
        EYcnaca{n} = sum([Enaca.*Ycsnaca.*Asnaca,Est.*Ystt.*Ast,Espar.*Ycspar.*Aspar])/EAnaca;
        EIxxnaca{n} = sum([Enaca.*tempIxxsnaca,Est.*tempIxxst,Espar.*tempIxxspar]);
        EIyynaca{n} = sum([Enaca.*tempIyysnaca,Est.*tempIyyst,Espar.*tempIyyspar]);
        EIxynaca{n} = sum([Enaca.*tempIxysnaca,Est.*tempIxyst,Espar.*tempIxyspar]);
    end  
       %% NACA Plot
    skin1 = plot(Xnaca,Ynaca1,'k');
    skin2 = plot(Xnaca,Ynaca2,'k');
    centroid = plot(Xcnaca{1},Ycnaca{2},'r*');
    modcentroid = plot(EXcnaca{1}, EYcnaca{2},'ro');
    legend([vsparplot stringerplot skin1 centroid modcentroid],'Vertical Spars','Stringers',...
        'Skin','Centroid','Modulus-weighted centroid')
    xlabel('X (inches)')
    ylabel('Y (inches)')
    %% Echo & Units
    Tablestringersgen = [Xst' Yst' Ast' Ixxst' Iyyst' Ixyst' densst' Est'];
    % Units
    xlswrite(outFile,Units3,2,'F32');
    xlswrite(outFile,Units1,2,'F75');
    xlswrite(outFile,Units2,2,'K75');
    xlswrite(outFile,Units1,2,'F84');
    xlswrite(outFile,Units1,2,'F110');
    xlswrite(outFile,Units2,2,'K110');
    xlswrite(outFile,Units1,2,'F119');
    xlswrite(outFile,Units1,2,'F130');
    xlswrite(outFile,Units2,2,'K130');
    xlswrite(outFile,Units1,2,'F139');
    xlswrite(outFile,Units1,2,'F150');
    xlswrite(outFile,Units2,2,'K150');
    xlswrite(outFile,Units1,2,'F159');
    % Echo to Airfoil Section (Sheet 2)
    createFigure(outFile,2,1,'C91','L104')
    xlswrite(outFile,Name,2,'D7');
    xlswrite(outFile,PID,2,'D8');
    xlswrite(outFile,title,2,'D10');
    xlswrite(outFile,[unitsin unitsout]',2,'E15');
    xlswrite(outFile,[isection x0 y0]',2,'E19');
    xlswrite(outFile,NACA,2,'E27');
    xlswrite(outFile,[numnacaseg c tnaca densnaca Enaca]',2,'E31');
    xlswrite(outFile,[xoc' tspar' densspar' Espar'],2,'C39');
    xlswrite(outFile,Tablestringersgen,2,'D47');
    %% Write Answers
    for n=1:4
        TableAreasnaca{n} = [Anaca Xcnaca{n} Ycnaca{n} Ixxnaca{n} Iyynaca{n} Ixynaca{n}]';
        TableMassesnaca{n} = [wAnaca wXcnaca{n} wYcnaca{n} wIxxnaca{n} wIyynaca{n} wIxynaca{n}]';
        TableStiffsnaca{n} = [EAnaca EXcnaca{n} EYcnaca{n} EIxxnaca{n} EIyynaca{n} EIxynaca{n}]';
    end
    xlswrite(outFile,TableAreasnaca{1},2,'E75');
    xlswrite(outFile,TableAreasnaca{2},2,'E110');
    xlswrite(outFile,TableAreasnaca{3},2,'E130');
    xlswrite(outFile,TableAreasnaca{4},2,'E150');
    xlswrite(outFile,TableMassesnaca{1},2,'J75');
    xlswrite(outFile,TableMassesnaca{2},2,'J110');
    xlswrite(outFile,TableMassesnaca{3},2,'J130');
    xlswrite(outFile,TableMassesnaca{4},2,'J150');
    xlswrite(outFile,TableStiffsnaca{1},2,'E84');
    xlswrite(outFile,TableStiffsnaca{2},2,'E119');
    xlswrite(outFile,TableStiffsnaca{3},2,'E139');
    xlswrite(outFile,TableStiffsnaca{4},2,'E159');
end