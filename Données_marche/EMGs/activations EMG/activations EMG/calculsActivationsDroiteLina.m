%% On commence par tout charger
% On a les muscles droits de Lina

biceps=load('BicepsDroitLina'); biceps=biceps.BicepsDroitLina;
droitAnterieur=load('DroitAnterieurDroitLina'); droitAnterieur=droitAnterieur.DroitAnterieurDroitLina;
gastrocnemien=load('GastrocnemienDroitLina'); gastrocnemien=gastrocnemien.GastrocnemienDroitLina;
gluteal=load('GlutealDroitLina'); gluteal=gluteal.GlutealDroitLina;
soleaire=load('SoleaireDroitLina'); soleaire=soleaire.SoleaireDroitLina;
tibialAnt=load('TibialAnterieurDroitLina'); tibialAnt=tibialAnt.TADroitLina;
tendineux=load('TendineuxDroitLina'); tendineux=tendineux.TendineuxDroitLina; tendineux{34,2}=tendineux{34,2}(:, 1:2);
vaste=load('VasteDroitLina'); vaste=vaste.VasteDroitLina;

%%
nbCyclesDroit=77;
x=0:1:100;

%% Biceps 

nbActivations=1;

for i=1:nbCyclesDroit
    if size(biceps{i,2},1)>nbActivations
        nbActivations=size(biceps{i,2},1);
    end
end
activationsBiceps=cell(nbActivations, 1);

for j=1:nbActivations 
    A=zeros(j,2);
    compt=1;
    for i=1:nbCyclesDroit
        if size(biceps{i,2},1)==j
            A(:,:,compt)=biceps{i,2};
            compt=compt+1;
        end
    end
    activationsBiceps{j, 1}=A;                
end

activationsMoyennesBiceps=cell(nbActivations, 2); % moyenne en position 1, déviation standard en position 2
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
            

%% Droit Antérieur

nbActivations=1;

for i=1:nbCyclesDroit
    if size(droitAnterieur{i,2},1)>nbActivations
        nbActivations=size(droitAnterieur{i,2},1);
    end
end
activationsdroitAnterieur=cell(nbActivations, 1);

for j=1:nbActivations 
    A=zeros(j,2);
    compt=1;
    for i=1:nbCyclesDroit
        if size(droitAnterieur{i,2},1)==j
            A(:,:,compt)=droitAnterieur{i,2};
            compt=compt+1;
        end
    end
    activationsdroitAnterieur{j, 1}=A;                
end

activationsMoyennesdroitAnterieur=cell(nbActivations, 2); % moyenne en position 1, dÃ©viation standard en position 2
for j=1:nbActivations
    activationsMoyennesdroitAnterieur{j, 1}=zeros(j,2);
    activationsMoyennesdroitAnterieur{j, 2}=zeros(j,2);
    
    for i=1:j 
        activationsMoyennesdroitAnterieur{j,1}(i,1)=mean(activationsdroitAnterieur{j,1}(i,1,:));
        activationsMoyennesdroitAnterieur{j,1}(i,2)=mean(activationsdroitAnterieur{j,1}(i,2,:));

        activationsMoyennesdroitAnterieur{j,2}(i,1)=std(activationsdroitAnterieur{j,1}(i,1,:));
        activationsMoyennesdroitAnterieur{j,2}(i,2)=std(activationsdroitAnterieur{j,1}(i,2,:));
    end
    
end  



%% Gastrocnémien


nbActivations=1;

for i=1:nbCyclesDroit
    if size(gastrocnemien{i,2},1)>nbActivations
        nbActivations=size(gastrocnemien{i,2},1);
    end
end
activationsGastrocnemien=cell(nbActivations, 1);

for j=1:nbActivations 
    A=zeros(j,2);
    compt=1;
    for i=1:nbCyclesDroit
        if size(gastrocnemien{i,2},1)==j
            A(:,:,compt)=gastrocnemien{i,2};
            compt=compt+1;
        end
    end
    activationsGastrocnemien{j, 1}=A;                
end

activationsMoyennesGastrocnemien=cell(nbActivations, 2); % moyenne en position 1, dÃ©viation standard en position 2
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

for i=1:nbCyclesDroit
    if size(gluteal{i,2},1)>nbActivations
        nbActivations=size(gluteal{i,2},1);
    end
end
activationsGluteal=cell(nbActivations, 1);
%%
for j=1:nbActivations 
    A=zeros(j,2);
    compt=1;
    for i=1:nbCyclesDroit
        if size(gluteal{i,2},1)==j
            A(:,:,compt)=gluteal{i,2};
            compt=compt+1;
        end
    end
    activationsGluteal{j, 1}=A;                
end

activationsGluteal{3,1}=activationsGluteal{3,1}(:,:,1:8);
%%
activationsMoyennesGluteal=cell(nbActivations, 2); % moyenne en position 1, dÃ©viation standard en position 2
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


%% Soléaire

nbActivations=1;

for i=1:nbCyclesDroit
    if size(soleaire{i,2},1)>nbActivations
        nbActivations=size(soleaire{i,2},1);
    end
end
activationsSoleaire=cell(nbActivations, 1);

for j=1:nbActivations 
    A=zeros(j,2);
    compt=1;
    for i=1:nbCyclesDroit
        if size(soleaire{i,2},1)==j
            A(:,:,compt)=soleaire{i,2};
            compt=compt+1;
        end
    end
    activationsSoleaire{j, 1}=A;                
end

activationsMoyennesSoleaire=cell(nbActivations, 2); % moyenne en position 1, dÃ©viation standard en position 2
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
    





%% Tibial antérieur

nbActivations=1;

for i=1:nbCyclesDroit
    if size(tibialAnt{i,2},1)>nbActivations
        nbActivations=size(tibialAnt{i,2},1);
    end
end
activationsTibialAnt=cell(nbActivations, 1);

for j=1:nbActivations 
    A=zeros(j,2);
    compt=1;
    for i=1:nbCyclesDroit
        if size(tibialAnt{i,2},1)==j
            A(:,:,compt)=tibialAnt{i,2};
            compt=compt+1;
        end
    end
    activationsTibialAnt{j, 1}=A;                
end

activationsMoyennesTibialAnt=cell(nbActivations, 2); % moyenne en position 1, dÃ©viation standard en position 2
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

for i=1:nbCyclesDroit
    if size(tendineux{i,2},1)>nbActivations
        nbActivations=size(tendineux{i,2},1);
    end
end

activationsTendineux=cell(nbActivations, 1);


for j=1:nbActivations 
    A=zeros(j,2);
    compt=1;
    for i=1:nbCyclesDroit
        if size(tendineux{i,2},1)==j
            A(:,:,compt)=tendineux{i,2};
            compt=compt+1;
        end
    end
    activationsTendineux{j, 1}=A;                
end


activationsMoyennesTendineux=cell(nbActivations, 2); % moyenne en position 1, dÃ©viation standard en position 2
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

for i=1:nbCyclesDroit
    if size(vaste{i,2},1)>nbActivations
        nbActivations=size(vaste{i,2},1);
    end
end
activationsVaste=cell(nbActivations, 1);

for j=1:nbActivations 
    A=zeros(j,2);
    compt=1;
    for i=1:nbCyclesDroit
        if size(vaste{i,2},1)==j
            A(:,:,compt)=vaste{i,2};
            compt=compt+1;
        end
    end
    activationsVaste{j, 1}=A;                
end

activationsMoyennesVaste=cell(nbActivations, 2); % moyenne en position 1, dÃ©viation standard en position 2
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

%% bidouillage biceps 2 activations

compt=1;
% 
K=zeros(2,2);
for i=1:77
    if size(biceps{i,2},1)==2
        K(:,:,compt)=biceps{i,2};
        compt=compt+1;
    end
end

%seuil à 65
comptL=1;
comptM=1;
comptN=1;
L=zeros(2,2); % offset 1 et onset 2 > seuil
M=zeros(2,2); % offset 1 < seuil et onset 2 > seuil
N=zeros(2,2); % offset 1 et onset 2 < seuil

for i=1:28
    if K(1,2,i)>65 && K(2,1,i)>65
        L(:,:,comptL)=K(:,:,i);
        comptL=comptL+1;
    elseif K(1,2,i)<65 && K(2,1,i)>65
        M(:,:,comptM)=K(:,:,i);
        comptM=comptM+1;   
    else 
        N(:,:,comptN)=K(:,:,i);
        comptN=comptN+1;
    end
end
        
%% bidouillage biceps 3 activations

compt=1;
% 
K=zeros(3,2);
for i=1:77
    if size(biceps{i,2},1)==3
        K(:,:,compt)=biceps{i,2};
        compt=compt+1;
    end
end

