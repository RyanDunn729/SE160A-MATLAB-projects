%% Preface %%

clear all; close all; clc;
format long;
Name = {'Ryan Dunn'};
PID = {'A15600858'};

%% Input Processing %%
inFile  = 'SE160A_2_Metallics_Input.xlsx';
outFile = 'SE160A_2_Metallics_Output.xlsx';
% deleteFigure(outFile);
insheet = readcell(inFile,'Sheet',1);
title = insheet(7,4);
unitsout = insheet{10,5};
Analysis = insheet{12,5};

if Analysis==1
    ABasis = cell2mat(insheet(16:23,5));
    BBasis = cell2mat(insheet(16:23,6));
    E = cell2mat(insheet(16,5:6));
    G = cell2mat(insheet(17,5:6));
    StrsYT = cell2mat(insheet(18,5:6));
    StrsUT = cell2mat(insheet(19,5:6));
    StrsYC = cell2mat(insheet(20,5:6));
    StrsUC = cell2mat(insheet(21,5:6));
    ShrY = cell2mat(insheet(22,5:6));
    ShrU = cell2mat(insheet(23,5:6));
    SFy = insheet{27,5};
    SFu = insheet{28,5};
    SState = cell2mat(insheet(35:40,5));
    Sxx = insheet{35,5};
    Syy = insheet{36,5};
    Szz = insheet{37,5};
    Syz = insheet{38,5};
    Sxz = insheet{39,5};
    Sxy = insheet{40,5};
    %% Principal Stresses
    Smat = [Sxx Sxy Sxz;
            Sxy Syy Syz;
            Sxz Syz Szz];
    [Phi,Sp] = eig(Smat);
    Sp = [Sp(1,1), Sp(2,2), Sp(3,3)];
    Tmax = (max(Sp)-min(Sp))/2;
    
    %% 3d Mohr Plot
    indexmin = 1;
    indexmid = 2;
    indexmax = 3;
    Center1 = [mean([max(Sp), min(Sp)]),0];
    Radius1 = abs((max(Sp)-min(Sp))/2);
    Center2 = [mean([Sp(indexmid), min(Sp)]),0];
    Radius2 = abs((Sp(indexmid)-min(Sp))/2);
    Center3 = [mean([max(Sp), Sp(indexmid)]),0];
    Radius3 = abs((max(Sp)-Sp(indexmid))/2);
    
    figure(1);
    plot([0 0], [-Tmax Tmax], '-k');
    hold on;
    plot([min(Sp) max(Sp)], [0 0], '-ko');
    plot([min(Sp) Sp(indexmid)], [0 0], '-ko');
    axis equal;
    grid on;
    xlabel('sigma (Ksi)'), ylabel('tau (Ksi)');
    viscircles(Center1, Radius1);
    viscircles(Center2, Radius2);
    viscircles(Center3, Radius3);  

    %% Allowables
    Ststar = [min([StrsYT(1)/SFy StrsUT(1)/SFu]) min([StrsYT(2)/SFy StrsUT(2)/SFu])];
    Scstar = [max([StrsYC(1)/SFy StrsUC(1)/SFu]) max([StrsYC(2)/SFy StrsUC(2)/SFu])];
    Shrstar = [min([ShrY(1)/SFy ShrU(1)/SFu]) min([ShrY(2)/SFy ShrU(2)/SFu])];
    Shrmixstar = [min([Shrstar(1), (Ststar(1)-Scstar(1))/4]),...
                  min([Shrstar(2), (Ststar(2)-Scstar(2))/4])];

    %% Margin of Safety
    % 1 = Rankine
    % 2 = Tresca
    % 3 = Von Mises
    MSA = zeros(1,3);
    MSB = zeros(1,3);
    RAmat = [Ststar(1)/max(Sp), abs(Scstar(1)/min(Sp))];
    RBmat = [Ststar(2)/max(Sp), abs(Scstar(2)/min(Sp))];
    MSA(1) = min(RAmat)-1;
    MSB(1) = min(RBmat)-1;
    TAmat = [Ststar(1)/max(Sp), abs(Scstar(1)/min(Sp)), (2*(Ststar(1)-Scstar(1))/4)/(max(Sp)-min(Sp))];
    TBmat = [Ststar(2)/max(Sp), abs(Scstar(2)/min(Sp)), (2*(Ststar(2)-Scstar(2))/4)/(max(Sp)-min(Sp))];
    MSA(2) = min(TAmat)-1;
    MSB(2) = min(TBmat)-1;
    eff = sqrt(((Sp(1)-Sp(2))^2+(Sp(2)-Sp(3))^2+(Sp(1)-Sp(3))^2)/2);
    VMA = Ststar(1)/eff;
    VMB = Ststar(2)/eff;
    MSA(3) = VMA-1;
    MSB(3) = VMB-1;

    %% A-Basis stress states
    % Rankine
    ARmat = min(RAmat)*Smat;
    % Tresca
    ATmat = min(TAmat)*Smat;
    % Von Mises
    AVMmat = VMA*Smat;

    %% B-Basis stress states
    BRmat = min(RBmat)*Smat;
    % Tresca
    BTmat = min(TBmat)*Smat;
    % Von Mises
    BVMmat = VMB*Smat;

    %% Echo & Answers
    Modulus = 'Msi';
    Stress = 'Ksi';
    Units1 = {Modulus Modulus Stress Stress Stress Stress Stress Stress}';
    Units2 = {Stress Stress Stress Stress Stress Stress}';
    xlswrite(outFile,Units1,1,'G21');
    xlswrite(outFile,Units2,1,'F37');
    xlswrite(outFile,{Stress},1,'H50');
    xlswrite(outFile,{Stress},1,'F56');
    xlswrite(outFile,Units2(1:4),1,'G85');
    xlswrite(outFile,Units2,1,'H98');
    xlswrite(outFile,Units2,1,'H107');
    ARoutput = [ARmat(1,1) ARmat(2,2) ARmat(3,3) ARmat(2,3) ARmat(1,3) ARmat(2,1)]';
    AToutput = [ATmat(1,1) ATmat(2,2) ATmat(3,3) ATmat(2,3) ATmat(1,3) ATmat(2,1)]';
    AVMoutput = [AVMmat(1,1) AVMmat(2,2) AVMmat(3,3) AVMmat(2,3) AVMmat(1,3) AVMmat(2,1)]';
    BRoutput = [BRmat(1,1) BRmat(2,2) BRmat(3,3) BRmat(2,3) BRmat(1,3) BRmat(2,1)]';
    BToutput = [BTmat(1,1) BTmat(2,2) BTmat(3,3) BTmat(2,3) BTmat(1,3) BTmat(2,1)]';
    BVMoutput = [BVMmat(1,1) BVMmat(2,2) BVMmat(3,3) BVMmat(2,3) BVMmat(1,3) BVMmat(2,1)]';
    % Sheet 2
    xlswrite(outFile,Name,2,'D7');
    xlswrite(outFile,PID,2,'D8');
    xlswrite(outFile,title,2,'D10');
    xlswrite(outFile,Analysis,2,'E17');
    xlswrite(outFile,[ABasis BBasis],2,'E21');
    xlswrite(outFile,[SFy SFu]',2,'E32');
    xlswrite(outFile,Units1,2,'G21');
    % Echo & Write Answers
    xlswrite(outFile,Name,1,'D7');
    xlswrite(outFile,PID,1,'D8');
    xlswrite(outFile,title,1,'D10');
    xlswrite(outFile,Analysis,1,'E17');
    xlswrite(outFile,[ABasis BBasis],1,'E21');
    xlswrite(outFile,[SFy; SFu],1,'E32');
    xlswrite(outFile,SState,1,'E37');
    xlswrite(outFile,Sp,1,'E50');
    xlswrite(outFile,Phi,1,'E51');
    xlswrite(outFile,Tmax,1,'E56');
    createFigure(outFile,1,1,'C58','H80');
    xlswrite(outFile,[Ststar; Scstar; Shrstar; Shrmixstar],1,'E85');
    xlswrite(outFile,[MSA; MSB],1,'E93');
    xlswrite(outFile,[ARoutput AToutput AVMoutput],1,'E98');
    xlswrite(outFile,[BRoutput BToutput BVMoutput],1,'E107');

elseif Analysis==2
    ABasis = cell2mat(insheet(16:23,5));
    BBasis = cell2mat(insheet(16:23,6));
    E = cell2mat(insheet(16,5:6));
    G = cell2mat(insheet(17,5:6));
    StrsYT = cell2mat(insheet(18,5:6));
    StrsUT = cell2mat(insheet(19,5:6));
    StrsYC = cell2mat(insheet(20,5:6));
    StrsUC = cell2mat(insheet(21,5:6));
    ShrY = cell2mat(insheet(22,5:6));
    ShrU = cell2mat(insheet(23,5:6));
    SFy = insheet{27,5};
    SFu = insheet{28,5};
    Aang = insheet{47,5};
    Bang = insheet{48,5};
    Cang = insheet{49,5};
    Gang = insheet{50,5};
    StrnA = insheet{54,5};
    StrnB = insheet{55,5};
    StrnC = insheet{56,5};
    %% Applied Strain Case 2
    Angmat = [cosd(Aang)^2, sind(Aang)^2, cosd(Aang)*sind(Aang);...
              cosd(Bang)^2, sind(Bang)^2, cosd(Bang)*sind(Bang);...
              cosd(Cang)^2, sind(Cang)^2, cosd(Cang)*sind(Cang)];
    T = [cosd(Gang)^2, sind(Gang)^2, cosd(Gang)*sind(Gang);...
         sind(Gang)^2, cosd(Gang)^2, -cosd(Gang)*sind(Gang);...
         -2*cosd(Gang)*sind(Gang), 2*cosd(Gang)*sind(Gang), cosd(Gang)^2-sind(Gang)^2];
    Smat = inv(Angmat)*([StrnA, StrnB, StrnC]');
    StrucFrame = inv(T)*Smat;
    Sxx = Smat(1);
    Syy = Smat(2);
    Sxy = Smat(3);
    Strnmat = [Sxx Sxy/2;
            Sxy/2 Syy];
    Strnp = eig(Strnmat);

    %% Strains
    StrnTstar = [min(StrsUT(1)/E(1)/SFu, StrsYT(1)/E(1)/SFy),...
                    min(StrsUT(2)/E(2)/SFu, StrsYT(2)/E(2)/SFy)] *10^3; % Ksi/Msi to in/in
    StrnCstar = [max(StrsUC(1)/E(1)/SFu, StrsYC(1)/E(1)/SFy),...
                    max(StrsUC(2)/E(2)/SFu, StrsYC(2)/E(2)/SFy)] *10^3; % Ksi/Msi to in/in
    StrnSstar = [min(ShrU(1)/G(1)/SFu, ShrY(1)/G(1)/SFy),...
                    min(ShrU(2)/G(2)/SFu, ShrY(2)/G(2)/SFy)] *10^3; % Ksi/Msi to in/in
    StrnSTstar = [min([StrnSstar(1), StrnTstar(1)/2, abs(StrnCstar(1)/2)]),...
                  min([StrnSstar(2), StrnTstar(2)/2, abs(StrnCstar(2)/2)])];
    MS = [min(StrnTstar(1)/max(Strnp), abs(StrnCstar(1)/min(Strnp)))-1,...
               min(StrnTstar(2)/max(Strnp), abs(StrnCstar(2)/min(Strnp)))-1];
    %% Write to Output
    Modulus = 'Msi';
    Stress = 'Ksi';
    Units1 = {Modulus Modulus Stress Stress Stress Stress Stress Stress}';
    xlswrite(outFile,Units1,2,'G21');
    Strnoutput = [StrnTstar; StrnCstar; StrnSstar; StrnSTstar];
    
    % Sheet 1
    xlswrite(outFile,Name,1,'D7');
    xlswrite(outFile,PID,1,'D8');
    xlswrite(outFile,title,1,'D10');
    xlswrite(outFile,Analysis,1,'E17');
    xlswrite(outFile,[ABasis BBasis],1,'E21');
    xlswrite(outFile,[SFy; SFu],1,'E32');
    xlswrite(outFile,Units1,1,'G21');
    % Echo & Answers
    xlswrite(outFile,Name,2,'D7');
    xlswrite(outFile,PID,2,'D8');
    xlswrite(outFile,title,2,'D10');
    xlswrite(outFile,Analysis,2,'E17');
    xlswrite(outFile,[ABasis BBasis],2,'E21');
    xlswrite(outFile,[SFy SFu]',2,'E32');
    xlswrite(outFile,[Aang Bang Cang Gang]',2,'E37');
    xlswrite(outFile,[StrnA StrnB StrnC]',2,'E44');
    xlswrite(outFile,StrucFrame,2,'E54');
    xlswrite(outFile,Strnp,2,'E61');
    xlswrite(outFile,Strnoutput,2,'E67');
    xlswrite(outFile,MS,2,'E75');

end