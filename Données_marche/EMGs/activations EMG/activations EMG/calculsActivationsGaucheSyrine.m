%% On commence par tout charger
% On a les muscles Gauches de Syrine

biceps=load('BicepsGaucheSyrine'); biceps=biceps.BicepsGaucheSyrine;
DroitAnterieur=load('DroitAnterieurGaucheSyrine'); DroitAnterieur=DroitAnterieur.DroitAnterieurGaucheSyrine;
gastrocnemien=load('GastrocnemienGaucheSyrine'); gastrocnemien=gastrocnemien.GastrocnemienGaucheSyrine;
gluteal=load('GlutealGaucheSyrine'); gluteal=gluteal.GlutealGaucheSyrine;
soleaire=load('SoleaireGaucheSyrine'); soleaire=soleaire.SoleaireGaucheSyrine;
tibialAnt=load('TibialAnterieurGaucheSyrine'); tibialAnt=tibialAnt.TAGaucheSyrine;
tendineux=load('TendineuxGaucheSyrine'); tendineux=tendineux.TendineuxGaucheSyrine; 
vaste=load('VasteGaucheSyrine'); vaste=vaste.VasteGaucheSyrine;

%%
nbCyclesGauche=50;
x=0:1:100;

%% Biceps 

nbActivations=1;

for i=1:nbCyclesGauche
    if size(biceps{i,2},1)>nbActivations
        nbActivations=size(biceps{i,2},1);
    end
end
activationsBiceps=cell(nbActivations, 1);

for j=1:nbActivations 
    A=zeros(j,2);
    compt=1;
    for i=1:nbCyclesGauche
        if size(biceps{i,2},1)==j
            A(:,:,compt)=biceps{i,2};
            compt=compt+1;
        end
    end
    activationsBiceps{j, 1}=A;                
end

activationsMoyennesBiceps=cell(nbActivations, 2); % moyenne en position 1, dÈviation standard en position 2
for j=1:nbActivations
    activationsMoyennesBiceps{j, 1}=zeros(j,2);
    activationsMoyennesBiceps{j, 2}=zeros(j,2);
    
    for i=1:j 
        activationsMoyennesBiceps{j,1}(i,1)=mean(activationsBiceps{j,1}(i,1,:));
        activationsMoyennesBiceps{j,1}(i,2)=mean(activationsBiceps{j,1}(i,2,:));

        activationsMoyennesBiceps{j,2}(i,1)=std(activationsBiceps{j,1}(i,1,:));
        activationsMoyennesBiceps{j,2}(i,2)=std(activationsBiceps{j,1}(i,2,:));
    end
    
end  
    
            

%% Droit AntÈrieur

nbActivations=1;

for i=1:nbCyclesGauche
    if size(DroitAnterieur{i,2},1)>nbActivations
        nbActivations=size(DroitAnterieur{i,2},1);
    end
end
activationsDroitAnterieur=cell(nbActivations, 1);

for j=1:nbActivations 
    A=zeros(j,2);
    compt=1;
    for i=1:nbCyclesGauche
        if size(DroitAnterieur{i,2},1)==j
            A(:,:,compt)=DroitAnterieur{i,2};
            compt=compt+1;
        end
    end
    activationsDroitAnterieur{j, 1}=A;                
end

activationsMoyennesDroitAnterieur=cell(nbActivations, 2); % moyenne en position 1, d√©viation standard en position 2
for j=1:nbActivations
    activationsMoyennesDroitAnterieur{j, 1}=zeros(j,2);
    activationsMoyennesDroitAnterieur{j, 2}=zeros(j,2);
    
    for i=1:j 
        activationsMoyennesDroitAnterieur{j,1}(i,1)=mean(activationsDroitAnterieur{j,1}(i,1,:));
        activationsMoyennesDroitAnterieur{j,1}(i,2)=mean(activationsDroitAnterieur{j,1}(i,2,:));

        activationsMoyennesDroitAnterieur{j,2}(i,1)=std(activationsDroitAnterieur{j,1}(i,1,:));
        activationsMoyennesDroitAnterieur{j,2}(i,2)=std(activationsDroitAnterieur{j,1}(i,2,:));
    end
    
end  
    



%% GastrocnÈmien


nbActivations=1;

for i=1:nbCyclesGauche
    if size(gastrocnemien{i,2},1)>nbActivations
        nbActivations=size(gastrocnemien{i,2},1);
    end
end
activationsGastrocnemien=cell(nbActivations, 1);

for j=1:nbActivations 
    A=zeros(j,2);
    compt=1;
    for i=1:nbCyclesGauche
        if size(gastrocnemien{i,2},1)==j
            A(:,:,compt)=gastrocnemien{i,2};
            compt=compt+1;
        end
    end
    activationsGastrocnemien{j, 1}=A;                
end

activationsMoyennesGastrocnemien=cell(nbActivations, 2); % moyenne en position 1, d√©viation standard en position 2
for j=1:nbActivations
    activationsMoyennesGastrocnemien{j, 1}=zeros(j,2);
    activationsMoyennesGastrocnemien{j, 2}=zeros(j,2);
    
    for i=1:j 
        activationsMoyennesGastrocnemien{j,1}(i,1)=mean(activationsGastrocnemien{j,1}(i,1,:));
        activationsMoyennesGastrocnemien{j,1}(i,2)=mean(activationsGastrocnemien{j,1}(i,2,:));

        activationsMoyennesGastrocnemien{j,2}(i,1)=std(activationsGastrocnemien{j,1}(i,1,:));
        activationsMoyennesGastrocnemien{j,2}(i,2)=std(activationsGastrocnemien{j,1}(i,2,:));
    end
    
end  
    




%% Gluteal

nbActivations=1;

for i=1:nbCyclesGauche
    if size(gluteal{i,2},1)>nbActivations
        nbActivations=size(gluteal{i,2},1);
    end
end
activationsGluteal=cell(nbActivations, 1);

for j=1:nbActivations 
    A=zeros(j,2);
    compt=1;
    for i=1:nbCyclesGauche
        if size(gluteal{i,2},1)==j
            A(:,:,compt)=gluteal{i,2};
            compt=compt+1;
        end
    end
    activationsGluteal{j, 1}=A;                
end

activationsMoyennesGluteal=cell(nbActivations, 2); % moyenne en position 1, d√©viation standard en position 2
for j=1:nbActivations
    activationsMoyennesGluteal{j, 1}=zeros(j,2);
    activationsMoyennesGluteal{j, 2}=zeros(j,2);
    
    for i=1:j 
        activationsMoyennesGluteal{j,1}(i,1)=mean(activationsGluteal{j,1}(i,1,:));
        activationsMoyennesGluteal{j,1}(i,2)=mean(activationsGluteal{j,1}(i,2,:));

        activationsMoyennesGluteal{j,2}(i,1)=std(activationsGluteal{j,1}(i,1,:));
        activationsMoyennesGluteal{j,2}(i,2)=std(activationsGluteal{j,1}(i,2,:));
    end
    
end  
    

%% SolÈaire

nbActivations=1;

for i=1:nbCyclesGauche
    if size(soleaire{i,2},1)>nbActivations
        nbActivations=size(soleaire{i,2},1);
    end
end
activationsSoleaire=cell(nbActivations, 1);

for j=1:nbActivations 
    A=zeros(j,2);
    compt=1;
    for i=1:nbCyclesGauche
        if size(soleaire{i,2},1)==j
            A(:,:,compt)=soleaire{i,2};
            compt=compt+1;
        end
    end
    activationsSoleaire{j, 1}=A;                
end

activationsMoyennesSoleaire=cell(nbActivations, 2); % moyenne en position 1, d√©viation standard en position 2
for j=1:nbActivations
    activationsMoyennesSoleaire{j, 1}=zeros(j,2);
    activationsMoyennesSoleaire{j, 2}=zeros(j,2);
    
    for i=1:j 
        activationsMoyennesSoleaire{j,1}(i,1)=mean(activationsSoleaire{j,1}(i,1,:));
        activationsMoyennesSoleaire{j,1}(i,2)=mean(activationsSoleaire{j,1}(i,2,:));

        activationsMoyennesSoleaire{j,2}(i,1)=std(activationsSoleaire{j,1}(i,1,:));
        activationsMoyennesSoleaire{j,2}(i,2)=std(activationsSoleaire{j,1}(i,2,:));
    end
    
end  



%% Tibial antÈrieur

nbActivations=1;

for i=1:nbCyclesGauche
    if size(tibialAnt{i,2},1)>nbActivations
        nbActivations=size(tibialAnt{i,2},1);
    end
end
activationsTibialAnt=cell(nbActivations, 1);

for j=1:nbActivations 
    A=zeros(j,2);
    compt=1;
    for i=1:nbCyclesGauche
        if size(tibialAnt{i,2},1)==j
            A(:,:,compt)=tibialAnt{i,2};
            compt=compt+1;
        end
    end
    activationsTibialAnt{j, 1}=A;                
end

activationsMoyennesTibialAnt=cell(nbActivations, 2); % moyenne en position 1, d√©viation standard en position 2
for j=1:nbActivations
    activationsMoyennesTibialAnt{j, 1}=zeros(j,2);
    activationsMoyennesTibialAnt{j, 2}=zeros(j,2);
    
    for i=1:j 
        activationsMoyennesTibialAnt{j,1}(i,1)=mean(activationsTibialAnt{j,1}(i,1,:));
        activationsMoyennesTibialAnt{j,1}(i,2)=mean(activationsTibialAnt{j,1}(i,2,:));

        activationsMoyennesTibialAnt{j,2}(i,1)=std(activationsTibialAnt{j,1}(i,1,:));
        activationsMoyennesTibialAnt{j,2}(i,2)=std(activationsTibialAnt{j,1}(i,2,:));
    end
    
end  
    


%% Tendineux


nbActivations=1;

for i=1:nbCyclesGauche
    if size(tendineux{i,2},1)>nbActivations
        nbActivations=size(tendineux{i,2},1);
    end
end

activationsTendineux=cell(nbActivations, 1);


for j=1:nbActivations 
    A=zeros(j,2);
    compt=1;
    for i=1:nbCyclesGauche
        if size(tendineux{i,2},1)==j
            A(:,:,compt)=tendineux{i,2};
            compt=compt+1;
        end
    end
    activationsTendineux{j, 1}=A;                
end

%%
activationsMoyennesTendineux=cell(nbActivations, 2); % moyenne en position 1, d√©viation standard en position 2
for j=1:nbActivations
    activationsMoyennesTendineux{j, 1}=zeros(j,2);
    activationsMoyennesTendineux{j, 2}=zeros(j,2);
    
    for i=1:j 
        activationsMoyennesTendineux{j,1}(i,1)=mean(activationsTendineux{j,1}(i,1,:));
        activationsMoyennesTendineux{j,1}(i,2)=mean(activationsTendineux{j,1}(i,2,:));

        activationsMoyennesTendineux{j,2}(i,1)=std(activationsTendineux{j,1}(i,1,:));
        activationsMoyennesTendineux{j,2}(i,2)=std(activationsTendineux{j,1}(i,2,:));
    end
    
end  
    



%% Vaste

nbActivations=1;

for i=1:nbCyclesGauche
    if size(vaste{i,2},1)>nbActivations
        nbActivations=size(vaste{i,2},1);
    end
end
activationsVaste=cell(nbActivations, 1);

for j=1:nbActivations 
    A=zeros(j,2);
    compt=1;
    for i=1:nbCyclesGauche
        if size(vaste{i,2},1)==j
            A(:,:,compt)=vaste{i,2};
            compt=compt+1;
        end
    end
    activationsVaste{j, 1}=A;                
end

activationsMoyennesVaste=cell(nbActivations, 2); % moyenne en position 1, d√©viation standard en position 2
for j=1:nbActivations
    activationsMoyennesVaste{j, 1}=zeros(j,2);
    activationsMoyennesVaste{j, 2}=zeros(j,2);
    
    for i=1:j 
        activationsMoyennesVaste{j,1}(i,1)=mean(activationsVaste{j,1}(i,1,:));
        activationsMoyennesVaste{j,1}(i,2)=mean(activationsVaste{j,1}(i,2,:));

        activationsMoyennesVaste{j,2}(i,1)=std(activationsVaste{j,1}(i,1,:));
        activationsMoyennesVaste{j,2}(i,2)=std(activationsVaste{j,1}(i,2,:));
    end
    
end  




%% calcul des signaux moyens

% on calcule le signal moyen

bicepsMoyen=zeros(101,2);

for j=2:nbActivations
    nbCycles=0;
    signalDeBase=zeros(101,1);
    for i=1:nbCyclesGauche
        if size(biceps{i,2}, 1)==j
            nbCycles=nbCycles+1;
            signalDeBase=signalDeBase+abs(biceps{i,1});
        end
    end
    bicepsMoyen(:, j-1)=signalDeBase/nbCycles;
end

vasteMoyen=zeros(101,2);

for j=2:nbActivations
    nbCycles=0;
    signalDeBase=zeros(101,1);
    for i=1:nbCyclesGauche
        if size(vaste{i,2}, 1)==j
            nbCycles=nbCycles+1;
            signalDeBase=signalDeBase+abs(vaste{i,1});
        end
    end
    vasteMoyen(:, j-1)=signalDeBase/nbCycles;
end

tendineuxMoyen=zeros(101,2);

for j=2:nbActivations
    nbCycles=0;
    signalDeBase=zeros(101,1);
    for i=1:nbCyclesGauche
        if size(tendineux{i,2}, 1)==j
            nbCycles=nbCycles+1;
            signalDeBase=signalDeBase+abs(tendineux{i,1});
        end
    end
    tendineuxMoyen(:, j-1)=signalDeBase/nbCycles;
end

tibialAntMoyen=zeros(101,2);

for j=2:nbActivations
    nbCycles=0;
    signalDeBase=zeros(101,1);
    for i=1:nbCyclesGauche
        if size(tibialAnt{i,2}, 1)==j
            nbCycles=nbCycles+1;
            signalDeBase=signalDeBase+abs(tibialAnt{i,1});
        end
    end
    tibialAntMoyen(:, j-1)=signalDeBase/nbCycles;
end


DroitAnterieurMoyen=zeros(101,2);

for j=2:nbActivations
    nbCycles=0;
    signalDeBase=zeros(101,1);
    for i=1:nbCyclesGauche
        if size(DroitAnterieur{i,2}, 1)==j
            nbCycles=nbCycles+1;
            signalDeBase=signalDeBase+abs(DroitAnterieur{i,1});
        end
    end
    DroitAnterieurMoyen(:, j-1)=signalDeBase/nbCycles;
end


gastrocnemienMoyen=zeros(101,2);

for j=2:nbActivations
    nbCycles=0;
    signalDeBase=zeros(101,1);
    for i=1:nbCyclesGauche
        if size(gastrocnemien{i,2}, 1)==j
            nbCycles=nbCycles+1;
            signalDeBase=signalDeBase+abs(gastrocnemien{i,1});
        end
    end
    gastrocnemienMoyen(:, j-1)=signalDeBase/nbCycles;
end

soleaireMoyen=zeros(101,2);

for j=2:nbActivations
    nbCycles=0;
    signalDeBase=zeros(101,1);
    for i=1:nbCyclesGauche
        if size(soleaire{i,2}, 1)==j
            nbCycles=nbCycles+1;
            signalDeBase=signalDeBase+abs(soleaire{i,1});
        end
    end
    soleaireMoyen(:, j-1)=signalDeBase/nbCycles;
end

glutealMoyen=zeros(101,2);

for j=2:nbActivations
    nbCycles=0;
    signalDeBase=zeros(101,1);
    for i=1:nbCyclesGauche
        if size(gluteal{i,2}, 1)==j
            nbCycles=nbCycles+1;
            signalDeBase=signalDeBase+abs(gluteal{i,1});
        end
    end
    glutealMoyen(:, j-1)=signalDeBase/nbCycles;
end




