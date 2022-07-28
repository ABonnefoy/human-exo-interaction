%% On commence par tout charger
% On a les muscles droits de Syrine

biceps=load('BicepsDroitSyrine'); biceps=biceps.BicepsDroitSyrine;
droitAnterieur=load('DroitAnterieurDroitSyrine'); droitAnterieur=droitAnterieur.DroitAnterieurDroitSyrine;
gastrocnemien=load('GastrocnemienDroitSyrine'); gastrocnemien=gastrocnemien.GastrocnemienDroitSyrine;
gluteal=load('GlutealDroitSyrine'); gluteal=gluteal.GlutealDroitSyrine;
soleaire=load('SoleaireDroitSyrine'); soleaire=soleaire.SoleaireDroitSyrine;
tibialAnt=load('TibialAnterieurDroitSyrine'); tibialAnt=tibialAnt.TADroitSyrine;
tendineux=load('TendineuxDroitSyrine'); tendineux=tendineux.TendineuxDroitSyrine;
vaste=load('VasteDroitSyrine'); vaste=vaste.VasteDroitSyrine;
%%
nbCyclesDroit=50;
x=0:1:100;

%% Biceps 
nbActivations=1;

for i=1:nbCyclesDroit
    if size(biceps{i,2},1)>nbActivations
        nbActivations=size(biceps{i,2},1);
    end
end
activations=cell(nbActivations, 1);

for j=1:nbActivations 
    A=zeros(j,2);
    compt=1;
    for i=1:nbCyclesDroit
        if size(biceps{i,2},1)==j
            A(:,:,compt)=biceps{i,2};
            compt=compt+1;
        end
    end
    activations{j, 1}=A;                
end

activations_moyennes=cell(nbActivations, 2); % moyenne en position 1, déviation standard en position 2
for j=1:nbActivations
    activations_moyennes{j, 1}=zeros(j,2);
    activations_moyennes{j, 2}=zeros(j,2);
    
    for i=1:j 
        activations_moyennes{j,1}(i,1)=mean(activations{j,1}(i,1,:));
        activations_moyennes{j,1}(i,2)=mean(activations{j,1}(i,2,:));

        activations_moyennes{j,2}(i,1)=std(activations{j,1}(i,1,:));
        activations_moyennes{j,2}(i,2)=std(activations{j,1}(i,2,:));
    end
    
end  
    
% on calcule le signal moyen

bicepsMoyen=zeros(101,2);

for j=2:nbActivations
    nbCycles=0;
    signalDeBase=zeros(101,1);
    for i=1:nbCyclesDroit
        if size(biceps{i,2}, 1)==j
            nbCycles=nbCycles+1;
            signalDeBase=signalDeBase+abs(biceps{i,1});
        end
    end
    bicepsMoyen(:, j-1)=signalDeBase/nbCycles;
end
            













































