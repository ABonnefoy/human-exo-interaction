% Pour les muscles droits de Syrine

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ON COMMENCE PAR TOUT CALCULER
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

statique=load('StatiqueSyrine');



%% Syrine, muscle 1 : gastrocnémien médial droit

cyclesGastrocnemienDroitSyrine=cell(0);  
activationsGastrocnemienDroitSyrine=cell(0);

for i=1:18   % itération sur toutes les marches
    filename_activations=strcat(['activationsCycleDroit_marche', num2str(i)]);
    data_activations=load(filename_activations);
    filename_cycles=strcat(['cyclesDroit_marche', num2str(i)]);
    data_cycles=load(filename_cycles);
    cyclesGastrocnemienDroitSyrine=[cyclesGastrocnemienDroitSyrine; data_cycles.cyclesDroit{1,1}];
    activationsGastrocnemienDroitSyrine=[activationsGastrocnemienDroitSyrine; data_activations.activationsCycleDroit{1,1}];
end

nbCyclesDroit=size(cyclesGastrocnemienDroitSyrine,1);



%% on plot tout pour vérifier que les activations sont ok

baselineGastrocnemien=mean(statique.StatiqueSyrine{1,1}(1:3000,1))+3*std(statique.StatiqueSyrine{1,1}(1:3000,1)); % =0.0285

for i=1:nbCyclesDroit
    figure;
    plot(abs(cyclesGastrocnemienDroitSyrine{i,1}))
    yline(baselineGastrocnemien);
    for u = 1:size(activationsGastrocnemienDroitSyrine{i,1},1)
        line(activationsGastrocnemienDroitSyrine{i,1}(u,1)*[1 1],[0 max(cyclesGastrocnemienDroitSyrine{i,1})],'color','r')
        line(activationsGastrocnemienDroitSyrine{i,1}(u,2)*[1 1],[0 max(cyclesGastrocnemienDroitSyrine{i,1})],'color','k')
    end
   
end

% SI BESOIN ON MODIFIE LES ACTIVATIONS A CE NIVEAU LA

%%

% on met tout entre 0 et 100
GastrocnemienDroitSyrine=cell(nbCyclesDroit, 2);


for i=1:nbCyclesDroit
    value=cyclesGastrocnemienDroitSyrine{i,1};
    s=size(value,1);
    cleanValue=interp1(0:1/(size(value,1)-1):1,value,0:1/100:1)';
    GastrocnemienDroitSyrine{i,1}=cleanValue;
    GastrocnemienDroitSyrine{i,2}=round(activationsGastrocnemienDroitSyrine{i,1}*100/s, 2); % on arrondit à deux décimales
end

x=0:1:100;
for i=1:nbCyclesDroit
    figure;
    plot(x, abs(GastrocnemienDroitSyrine{i,1}))
    for u = 1:size(activationsGastrocnemienDroitSyrine{i,1},1)
        line(GastrocnemienDroitSyrine{i,2}(u,1)*[1 1],[0 max(cyclesGastrocnemienDroitSyrine{i,1})],'color','r')
        line(GastrocnemienDroitSyrine{i,2}(u,2)*[1 1],[0 max(cyclesGastrocnemienDroitSyrine{i,1})],'color','k')
    end
   
end

%% Syrine muscle 2 : Droit antérieur Droit

cyclesDroitAnterieurDroitSyrine=cell(0);  
activationsDroitAnterieurDroitSyrine=cell(0);

for i=1:18   % itération sur toutes les marches
    filename_activations=strcat(['activationsCycleDroit_marche', num2str(i)]);
    data_activations=load(filename_activations);
    filename_cycles=strcat(['cyclesDroit_marche', num2str(i)]);
    data_cycles=load(filename_cycles);
    cyclesDroitAnterieurDroitSyrine=[cyclesDroitAnterieurDroitSyrine; data_cycles.cyclesDroit{3,1}];
    activationsDroitAnterieurDroitSyrine=[activationsDroitAnterieurDroitSyrine; data_activations.activationsCycleDroit{3,1}];
end

%% on plot tout pour vérifier que les activations sont ok

baselineDroitAnterieur=mean(statique.StatiqueSyrine{1,1}(3000:4980,3))+3*std(statique.StatiqueSyrine{1,1}(3000:4980,3)); % =0.0363

for i=1:nbCyclesDroit
    figure;
    plot(abs(cyclesDroitAnterieurDroitSyrine{i,1}))
    yline(baselineDroitAnterieur);
    for u = 1:size(activationsDroitAnterieurDroitSyrine{i,1},1)
        line(activationsDroitAnterieurDroitSyrine{i,1}(u,1)*[1 1],[0 max(cyclesDroitAnterieurDroitSyrine{i,1})],'color','r')
        line(activationsDroitAnterieurDroitSyrine{i,1}(u,2)*[1 1],[0 max(cyclesDroitAnterieurDroitSyrine{i,1})],'color','k')
    end
   
end

% SI BESOIN ON MODIFIE LES ACTIVATIONS A CE NIVEAU LA

%% on met tout entre 0 et 100
DroitAnterieurDroitSyrine=cell(nbCyclesDroit, 2);


for i=1:nbCyclesDroit
    value=cyclesDroitAnterieurDroitSyrine{i,1};
    s=size(value,1);
    cleanValue=interp1(0:1/(size(value,1)-1):1,value,0:1/100:1)';
    DroitAnterieurDroitSyrine{i,1}=cleanValue;
    DroitAnterieurDroitSyrine{i,2}=round(activationsDroitAnterieurDroitSyrine{i,1}*100/s, 2); % on arrondit à deux décimales
end

% for i=1:nbCyclesDroit
%     figure;
%     plot(x, abs(DroitAnterieurDroitSyrine{i,1}))
%     for u = 1:size(activationsDroitAnterieurDroitSyrine{i,1},1)
%         line(DroitAnterieurDroitSyrine{i,2}(u,1)*[1 1],[0 max(cyclesDroitAnterieurDroitSyrine{i,1})],'color','r')
%         line(DroitAnterieurDroitSyrine{i,2}(u,2)*[1 1],[0 max(cyclesDroitAnterieurDroitSyrine{i,1})],'color','k')
%     end
%    
% end

%% Syrine muscle 3 : Vaste externe droit

cyclesVasteDroitSyrine=cell(0);  
activationsVasteDroitSyrine=cell(0);

for i=1:18   % itération sur toutes les marches
    filename_activations=strcat(['activationsCycleDroit_marche', num2str(i)]);
    data_activations=load(filename_activations);
    filename_cycles=strcat(['cyclesDroit_marche', num2str(i)]);
    data_cycles=load(filename_cycles);
    cyclesVasteDroitSyrine=[cyclesVasteDroitSyrine; data_cycles.cyclesDroit{5,1}];
    activationsVasteDroitSyrine=[activationsVasteDroitSyrine; data_activations.activationsCycleDroit{5,1}];
end

%% on plot tout pour vérifier que les activations sont ok

baselineVaste=mean(statique.StatiqueSyrine{1,1}(:,5))+3*std(statique.StatiqueSyrine{1,1}(:,5)); % =0.0710

for i=1:nbCyclesDroit
    figure;
    plot(abs(cyclesVasteDroitSyrine{i,1}))
    yline(baselineVaste);
    for u = 1:size(activationsVasteDroitSyrine{i,1},1)
        line(activationsVasteDroitSyrine{i,1}(u,1)*[1 1],[0 max(cyclesVasteDroitSyrine{i,1})],'color','r')
        line(activationsVasteDroitSyrine{i,1}(u,2)*[1 1],[0 max(cyclesVasteDroitSyrine{i,1})],'color','k')
    end
   
end

% SI BESOIN ON MODIFIE LES ACTIVATIONS A CE NIVEAU LA

%% on met tout entre 0 et 100
VasteDroitSyrine=cell(nbCyclesDroit, 2);


for i=1:nbCyclesDroit
    value=cyclesVasteDroitSyrine{i,1};
    s=size(value,1);
    cleanValue=interp1(0:1/(size(value,1)-1):1,value,0:1/100:1)';
    VasteDroitSyrine{i,1}=cleanValue;
    VasteDroitSyrine{i,2}=round(activationsVasteDroitSyrine{i,1}*100/s, 2); % on arrondit à deux décimales
end


%% Syrine muscle 4 : tibial antérieur droit

cyclesTADroitSyrine=cell(0);  
activationsTADroitSyrine=cell(0);

for i=1:18   % itération sur toutes les marches
    filename_activations=strcat(['activationsCycleDroit_marche', num2str(i)]);
    data_activations=load(filename_activations);
    filename_cycles=strcat(['cyclesDroit_marche', num2str(i)]);
    data_cycles=load(filename_cycles);
    cyclesTADroitSyrine=[cyclesTADroitSyrine; data_cycles.cyclesDroit{7,1}];
    activationsTADroitSyrine=[activationsTADroitSyrine; data_activations.activationsCycleDroit{7,1}];
end

%% on plot tout pour vérifier que les activations sont ok

baselineTA=mean(statique.StatiqueSyrine{1,1}(:,7))+3*std(statique.StatiqueSyrine{1,1}(:,7)); % =0.0411

for i=1:nbCyclesDroit
    figure;
    plot(abs(cyclesTADroitSyrine{i,1}))
    yline(baselineTA);
    for u = 1:size(activationsTADroitSyrine{i,1},1)
        line(activationsTADroitSyrine{i,1}(u,1)*[1 1],[0 max(cyclesTADroitSyrine{i,1})],'color','r')
        line(activationsTADroitSyrine{i,1}(u,2)*[1 1],[0 max(cyclesTADroitSyrine{i,1})],'color','k')
    end
   
end

% SI BESOIN ON MODIFIE LES ACTIVATIONS A CE NIVEAU LA

%% on met tout entre 0 et 100
TADroitSyrine=cell(nbCyclesDroit, 2);


for i=1:nbCyclesDroit
    value=cyclesTADroitSyrine{i,1};
    s=size(value,1);
    cleanValue=interp1(0:1/(size(value,1)-1):1,value,0:1/100:1)';
    TADroitSyrine{i,1}=cleanValue;
    TADroitSyrine{i,2}=round(activationsTADroitSyrine{i,1}*100/s, 2); % on arrondit à deux décimales
end

%% Syrine muscle 5: soléaire droit

cyclesSoleaireDroitSyrine=cell(0);  
activationsSoleaireDroitSyrine=cell(0);

for i=1:18   % itération sur toutes les marches
    filename_activations=strcat(['activationsCycleDroit_marche', num2str(i)]);
    data_activations=load(filename_activations);
    filename_cycles=strcat(['cyclesDroit_marche', num2str(i)]);
    data_cycles=load(filename_cycles);
    cyclesSoleaireDroitSyrine=[cyclesSoleaireDroitSyrine; data_cycles.cyclesDroit{9,1}];
    activationsSoleaireDroitSyrine=[activationsSoleaireDroitSyrine; data_activations.activationsCycleDroit{9,1}];
end

%% on plot tout pour vérifier que les activations sont ok

baselineSoleaire=mean(statique.StatiqueSyrine{1,1}(2500:4980,9))+3*std(statique.StatiqueSyrine{1,1}(2500:4980,9)); % =0.0411

for i=1:nbCyclesDroit
    figure;
    plot(abs(cyclesSoleaireDroitSyrine{i,1}))
    yline(baselineSoleaire);
    for u = 1:size(activationsSoleaireDroitSyrine{i,1},1)
        line(activationsSoleaireDroitSyrine{i,1}(u,1)*[1 1],[0 max(cyclesSoleaireDroitSyrine{i,1})],'color','r')
        line(activationsSoleaireDroitSyrine{i,1}(u,2)*[1 1],[0 max(cyclesSoleaireDroitSyrine{i,1})],'color','k')
    end
   
end

% SI BESOIN ON MODIFIE LES ACTIVATIONS A CE NIVEAU LA

%% on met tout entre 0 et 100
SoleaireDroitSyrine=cell(nbCyclesDroit, 2);


for i=1:nbCyclesDroit
    value=cyclesSoleaireDroitSyrine{i,1};
    s=size(value,1);
    cleanValue=interp1(0:1/(size(value,1)-1):1,value,0:1/100:1)';
    SoleaireDroitSyrine{i,1}=cleanValue;
    SoleaireDroitSyrine{i,2}=round(activationsSoleaireDroitSyrine{i,1}*100/s, 2); % on arrondit à deux décimales
end

%% Syrine muscle 6 : semi tendineux droit

cyclesTendineuxDroitSyrine=cell(0);  
activationsTendineuxDroitSyrine=cell(0);

for i=1:18   % itération sur toutes les marches
    filename_activations=strcat(['activationsCycleDroit_marche', num2str(i)]);
    data_activations=load(filename_activations);
    filename_cycles=strcat(['cyclesDroit_marche', num2str(i)]);
    data_cycles=load(filename_cycles);
    cyclesTendineuxDroitSyrine=[cyclesTendineuxDroitSyrine; data_cycles.cyclesDroit{11,1}];
    activationsTendineuxDroitSyrine=[activationsTendineuxDroitSyrine; data_activations.activationsCycleDroit{11,1}];
end

%% on plot tout pour vérifier que les activations sont ok

baselineTendineux=mean(statique.StatiqueSyrine{1,1}(1:3000,11))+3*std(statique.StatiqueSyrine{1,1}(1:3000,11)); % =0.0271

for i=1:nbCyclesDroit
    figure;
    plot(abs(cyclesTendineuxDroitSyrine{i,1}))
    yline(baselineTendineux);
    for u = 1:size(activationsTendineuxDroitSyrine{i,1},1)
        line(activationsTendineuxDroitSyrine{i,1}(u,1)*[1 1],[0 max(cyclesTendineuxDroitSyrine{i,1})],'color','r')
        line(activationsTendineuxDroitSyrine{i,1}(u,2)*[1 1],[0 max(cyclesTendineuxDroitSyrine{i,1})],'color','k')
    end
   
end

% SI BESOIN ON MODIFIE LES ACTIVATIONS A CE NIVEAU LA

%% on met tout entre 0 et 100
TendineuxDroitSyrine=cell(nbCyclesDroit, 2);


for i=1:nbCyclesDroit
    value=cyclesTendineuxDroitSyrine{i,1};
    s=size(value,1);
    cleanValue=interp1(0:1/(size(value,1)-1):1,value,0:1/100:1)';
    TendineuxDroitSyrine{i,1}=cleanValue;
    TendineuxDroitSyrine{i,2}=round(activationsTendineuxDroitSyrine{i,1}*100/s, 2); % on arrondit à deux décimales
end

%% Syrine muscle 7 : biceps femoral droit

cyclesBicepsDroitSyrine=cell(0);  
activationsBicepsDroitSyrine=cell(0);

for i=1:18   % itération sur toutes les marches
    filename_activations=strcat(['activationsCycleDroit_marche', num2str(i)]);
    data_activations=load(filename_activations);
    filename_cycles=strcat(['cyclesDroit_marche', num2str(i)]);
    data_cycles=load(filename_cycles);
    cyclesBicepsDroitSyrine=[cyclesBicepsDroitSyrine; data_cycles.cyclesDroit{13,1}];
    activationsBicepsDroitSyrine=[activationsBicepsDroitSyrine; data_activations.activationsCycleDroit{13,1}];
end

%% on plot tout pour vérifier que les activations sont ok

baselineBiceps=mean(statique.StatiqueSyrine{1,1}(1:3000,13))+3*std(statique.StatiqueSyrine{1,1}(1:3000,13)); % =0.0271

for i=1:nbCyclesDroit
    figure;
    plot(abs(cyclesBicepsDroitSyrine{i,1}))
    yline(baselineBiceps);
    for u = 1:size(activationsBicepsDroitSyrine{i,1},1)
        line(activationsBicepsDroitSyrine{i,1}(u,1)*[1 1],[0 max(cyclesBicepsDroitSyrine{i,1})],'color','r')
        line(activationsBicepsDroitSyrine{i,1}(u,2)*[1 1],[0 max(cyclesBicepsDroitSyrine{i,1})],'color','k')
    end
   
end

% SI BESOIN ON MODIFIE LES ACTIVATIONS A CE NIVEAU LA

%% on met tout entre 0 et 100
BicepsDroitSyrine=cell(nbCyclesDroit, 2);


for i=1:nbCyclesDroit
    value=cyclesBicepsDroitSyrine{i,1};
    s=size(value,1);
    cleanValue=interp1(0:1/(size(value,1)-1):1,value,0:1/100:1)';
    BicepsDroitSyrine{i,1}=cleanValue;
    BicepsDroitSyrine{i,2}=round(activationsBicepsDroitSyrine{i,1}*100/s, 2); % on arrondit à deux décimales
end


%% Syrine muscle 8 : grand fessier droit


cyclesGlutealDroitSyrine=cell(0);  
activationsGlutealDroitSyrine=cell(0);

for i=1:18   % itération sur toutes les marches
    filename_activations=strcat(['activationsCycleDroit_marche', num2str(i)]);
    data_activations=load(filename_activations);
    filename_cycles=strcat(['cyclesDroit_marche', num2str(i)]);
    data_cycles=load(filename_cycles);
    cyclesGlutealDroitSyrine=[cyclesGlutealDroitSyrine; data_cycles.cyclesDroit{15,1}];
    activationsGlutealDroitSyrine=[activationsGlutealDroitSyrine; data_activations.activationsCycleDroit{15,1}];
end

%% on plot tout pour vérifier que les activations sont ok

baselineGluteal=mean(statique.StatiqueSyrine{1,1}(1:3000,15))+3*std(statique.StatiqueSyrine{1,1}(1:3000,15)); % =0.0271

for i=1:nbCyclesDroit
    figure;
    plot(abs(cyclesGlutealDroitSyrine{i,1}))
    yline(baselineGluteal);
    for u = 1:size(activationsGlutealDroitSyrine{i,1},1)
        line(activationsGlutealDroitSyrine{i,1}(u,1)*[1 1],[0 max(cyclesGlutealDroitSyrine{i,1})],'color','r')
        line(activationsGlutealDroitSyrine{i,1}(u,2)*[1 1],[0 max(cyclesGlutealDroitSyrine{i,1})],'color','k')
    end
   
end

% SI BESOIN ON MODIFIE LES ACTIVATIONS A CE NIVEAU LA

%% on met tout entre 0 et 100
GlutealDroitSyrine=cell(nbCyclesDroit, 2);


for i=1:nbCyclesDroit
    value=cyclesGlutealDroitSyrine{i,1};
    s=size(value,1);
    cleanValue=interp1(0:1/(size(value,1)-1):1,value,0:1/100:1)';
    GlutealDroitSyrine{i,1}=cleanValue;
    GlutealDroitSyrine{i,2}=round(activationsGlutealDroitSyrine{i,1}*100/s, 2); % on arrondit à deux décimales
end


