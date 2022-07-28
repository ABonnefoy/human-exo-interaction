% Quantification des emg par cycle de marche
Lina=load('EMGLina_processed'); Lina=Lina.EMGLina_processed;
Syrine=load('EMGSyrine_processed'); Syrine=Syrine.EMGSyrine_processed;
toeOffMoyen;

fc_smooth=9;
fsAnalog=2000;
order=4;

for i=1:18
    data=Syrine{i,1};
    for j=1:16
        data1=butterworth(fsAnalog, fc_smooth, order, 'low', abs(data(:,j)));
        Syrine{i,1}(:,j)=data1;
    end
end

for i=1:9
    data=Lina{i,1};
    for j=1:16
        data1=butterworth(fsAnalog, fc_smooth, order, 'low', abs(data(:,j)));
        Lina{i,1}(:,j)=data1;
    end
end

eventsLina=load('eventsLina');
boundariesLina=load('boundariesLina');

cyclesDroitLina=cell(16, 9);
cyclesGaucheLina=cell(16, 9);

for i=1:9
    EMGChannels=Lina{i,1};
    events=eventsLina.eventsLina{i,1}; % events of the gait cycle
    rightFoot=(int32(events.Right_Foot_Strike*200)-boundariesLina.boundariesLina(i,1)+1)*10; % on les met dans le même système que fsAnalog
    leftFoot=(int32(events.Left_Foot_Strike*200)-boundariesLina.boundariesLina(i,1)+1)*10; % foot strikes
    sRight=size(rightFoot, 2);
    sLeft=size(leftFoot,2);
    for j=1:16
        cyclesDroitLina{j,i}=cell(sRight-1,1);
        cyclesGaucheLina{j,i}=cell(sLeft-1,1);
        muscle=[];
        for k=1:sRight-1
            muscle=EMGChannels(rightFoot(k):rightFoot(k+1), j); % pas les bons indices: décalage de x10
            cyclesDroitLina{j,i}{k,1}=muscle;
        end
        for k=1:sLeft-1
            muscle=EMGChannels(leftFoot(k):leftFoot(k+1), j);
            cyclesGaucheLina{j,i}{k,1}=muscle;
        end    
    end
end

eventsSyrine=load('eventsSyrine');
boundariesSyrine=load('boundariesSyrine');

cyclesDroitSyrine=cell(16, 18);
cyclesGaucheSyrine=cell(16, 18);

for i=1:18
    EMGChannels=Syrine{i,1};
    events=eventsSyrine.eventsSyrine{i,1}; % events of the gait cycle
    rightFoot=(int32(events.Right_Foot_Strike*200)-boundariesSyrine.boundariesSyrine(i,1)+1)*10; % on les met dans le même système que fsAnalog
    leftFoot=(int32(events.Left_Foot_Strike*200)-boundariesSyrine.boundariesSyrine(i,1)+1)*10; % foot strikes
    sRight=size(rightFoot, 2);
    sLeft=size(leftFoot,2);
    for j=1:16
        cyclesDroitSyrine{j,i}=cell(sRight-1,1);
        cyclesGaucheSyrine{j,i}=cell(sLeft-1,1);
        muscle=[];
        for k=1:sRight-1
            muscle=EMGChannels(rightFoot(k):rightFoot(k+1), j); % pas les bons indices: décalage de x10
            cyclesDroitSyrine{j,i}{k,1}=muscle;
        end
        for k=1:sLeft-1
            muscle=EMGChannels(leftFoot(k):leftFoot(k+1), j);
            cyclesGaucheSyrine{j,i}{k,1}=muscle;
        end    
    end
end


%% Interpolation entre 0 et 100

droitSyrineInterp=cell(16,18);
gaucheSyrineInterp=cell(16,18);

for i=1:18
    for j=1:16
        nbCycles=size(cyclesDroitSyrine{j,i}, 1);
        droitSyrineInterp{j,i}=cell(nbCycles, 1);
        for k=1:nbCycles
            value=cyclesDroitSyrine{j,i}{k,1};
            droitSyrineInterp{j,i}{k,1}=interp1(0:1/(size(value,1)-1):1,value,0:1/100:1)';
        end
        nbCycles=size(cyclesGaucheSyrine{j,i}, 1);
        gaucheSyrineInterp{j,i}=cell(nbCycles, 1);
        for k=1:nbCycles
            value=cyclesGaucheSyrine{j,i}{k,1};
            gaucheSyrineInterp{j,i}{k,1}=interp1(0:1/(size(value,1)-1):1,value,0:1/100:1)';
        end
    end
end
            



droitLinaInterp=cell(16,9);
gaucheLinaInterp=cell(16,18);

for i=1:9
    for j=1:16
        nbCycles=size(cyclesDroitLina{j,i}, 1);
        droitLinaInterp{j,i}=cell(nbCycles, 1);
        for k=1:nbCycles
            value=cyclesDroitLina{j,i}{k,1};
            droitLinaInterp{j,i}{k,1}=interp1(0:1/(size(value,1)-1):1,value,0:1/100:1)';
        end
        nbCycles=size(cyclesGaucheLina{j,i}, 1);
        gaucheLinaInterp{j,i}=cell(nbCycles, 1);
        for k=1:nbCycles
            value=cyclesGaucheLina{j,i}{k,1};
            gaucheLinaInterp{j,i}{k,1}=interp1(0:1/(size(value,1)-1):1,value,0:1/100:1)';
        end
    end
end
        

%% on fait la moyenne pour chaque muscle
% Lina droite

LDbiceps=[];
LDdroitAnt=[];
LDgastrocnemien=[];
LDgluteal=[];
LDsoleaire=[];
LDtendineux=[];
LDtibialAnt=[];
LDvaste=[];


for i=1:9
    s=size(droitLinaInterp{1,i},1);
    for j=1:s
        LDbiceps=[LDbiceps, droitLinaInterp{14, i}{j,1}];
        LDdroitAnt=[LDdroitAnt, droitLinaInterp{3, i}{j,1}];
        LDgastrocnemien=[LDgastrocnemien, droitLinaInterp{1, i}{j,1}];
        LDgluteal=[LDgluteal, droitLinaInterp{16, i}{j,1}];
        LDsoleaire=[LDsoleaire, droitLinaInterp{9, i}{j,1}];
        LDtendineux=[LDtendineux, droitLinaInterp{11, i}{j,1}];
        LDtibialAnt=[LDtibialAnt, droitLinaInterp{7, i}{j,1}];
        LDvaste=[LDvaste, droitLinaInterp{5, i}{j,1}];
    end
end

LDbicepsMoy=mean(LDbiceps,2);
LDdroitAntMoy=mean(LDdroitAnt,2);
LDgastrocnemienMoy=mean(LDgastrocnemien,2);
LDglutealMoy=mean(LDgluteal,2);
LDsoleaireMoy=mean(LDsoleaire,2);
LDtendineuxMoy=mean(LDtendineux,2);
LDtibialAntMoy=mean(LDtibialAnt,2);
LDvasteMoy=mean(LDvaste,2);


LDbicepsSd=std(LDbiceps,1,2);
LDdroitAntSd=std(LDdroitAnt,1,2);
LDgastrocnemienSd=std(LDgastrocnemien,1,2);
LDglutealSd=std(LDgluteal,1,2);
LDsoleaireSd=std(LDsoleaire,1,2);
LDtendineuxSd=std(LDtendineux,1,2);
LDtibialAntSd=std(LDtibialAnt,1,2);
LDvasteSd=std(LDvaste,1,2);


%% Lina gauche

LGbiceps=[];
LGdroitAnt=[];
LGgastrocnemien=[];
LGgluteal=[];
LGsoleaire=[];
LGtendineux=[];
LGtibialAnt=[];
LGvaste=[];


for i=1:9
    s=size(gaucheLinaInterp{1,i},1);
    for j=1:s
        LGbiceps=[LGbiceps, gaucheLinaInterp{13, i}{j,1}];
        LGdroitAnt=[LGdroitAnt, gaucheLinaInterp{4, i}{j,1}];
        LGgastrocnemien=[LGgastrocnemien, gaucheLinaInterp{2, i}{j,1}];
        LGgluteal=[LGgluteal, gaucheLinaInterp{15, i}{j,1}];
        LGsoleaire=[LGsoleaire, gaucheLinaInterp{10, i}{j,1}];
        LGtendineux=[LGtendineux, gaucheLinaInterp{12, i}{j,1}];
        LGtibialAnt=[LGtibialAnt, gaucheLinaInterp{8, i}{j,1}];
        LGvaste=[LGvaste, gaucheLinaInterp{6, i}{j,1}];
    end
end

LGbicepsMoy=mean(LGbiceps,2);
LGdroitAntMoy=mean(LGdroitAnt,2);
LGgastrocnemienMoy=mean(LGgastrocnemien,2);
LGglutealMoy=mean(LGgluteal,2);
LGsoleaireMoy=mean(LGsoleaire,2);
LGtendineuxMoy=mean(LGtendineux,2);
LGtibialAntMoy=mean(LGtibialAnt,2);
LGvasteMoy=mean(LGvaste,2);


LGbicepsSd=std(LGbiceps,1,2);
LGdroitAntSd=std(LGdroitAnt,1,2);
LGgastrocnemienSd=std(LGgastrocnemien,1,2);
LGglutealSd=std(LGgluteal,1,2);
LGsoleaireSd=std(LGsoleaire,1,2);
LGtendineuxSd=std(LGtendineux,1,2);
LGtibialAntSd=std(LGtibialAnt,1,2);
LGvasteSd=std(LGvaste,1,2);


%% Syrine Droite


SDbiceps=[];
SDdroitAnt=[];
SDgastrocnemien=[];
SDgluteal=[];
SDsoleaire=[];
SDtendineux=[];
SDtibialAnt=[];
SDvaste=[];


for i=1:18
    s=size(droitSyrineInterp{1,i},1);
    for j=1:s
        SDbiceps=[SDbiceps, droitSyrineInterp{13, i}{j,1}];
        SDdroitAnt=[SDdroitAnt, droitSyrineInterp{3, i}{j,1}];
        SDgastrocnemien=[SDgastrocnemien, droitSyrineInterp{1, i}{j,1}];
        SDgluteal=[SDgluteal, droitSyrineInterp{15, i}{j,1}];
        SDsoleaire=[SDsoleaire, droitSyrineInterp{9, i}{j,1}];
        SDtendineux=[SDtendineux, droitSyrineInterp{11, i}{j,1}];
        SDtibialAnt=[SDtibialAnt, droitSyrineInterp{7, i}{j,1}];
        SDvaste=[SDvaste, droitSyrineInterp{5, i}{j,1}];
    end
end

SDbicepsMoy=mean(SDbiceps,2);
SDdroitAntMoy=mean(SDdroitAnt,2);
SDgastrocnemienMoy=mean(SDgastrocnemien,2);
SDglutealMoy=mean(SDgluteal,2);
SDsoleaireMoy=mean(SDsoleaire,2);
SDtendineuxMoy=mean(SDtendineux,2);
SDtibialAntMoy=mean(SDtibialAnt,2);
SDvasteMoy=mean(SDvaste,2);


SDbicepsSd=std(SDbiceps,1,2);
SDdroitAntSd=std(SDdroitAnt,1,2);
SDgastrocnemienSd=std(SDgastrocnemien,1,2);
SDglutealSd=std(SDgluteal,1,2);
SDsoleaireSd=std(SDsoleaire,1,2);
SDtendineuxSd=std(SDtendineux,1,2);
SDtibialAntSd=std(SDtibialAnt,1,2);
SDvasteSd=std(SDvaste,1,2);


%% Syrine gauche

SGbiceps=[];
SGdroitAnt=[];
SGgastrocnemien=[];
SGgluteal=[];
SGsoleaire=[];
SGtendineux=[];
SGtibialAnt=[];
SGvaste=[];


for i=1:18
    s=size(gaucheSyrineInterp{1,i},1);
    for j=1:s
        SGbiceps=[SGbiceps, gaucheSyrineInterp{14, i}{j,1}];
        SGdroitAnt=[SGdroitAnt, gaucheSyrineInterp{4, i}{j,1}];
        SGgastrocnemien=[SGgastrocnemien, gaucheSyrineInterp{2, i}{j,1}];
        SGgluteal=[SGgluteal, gaucheSyrineInterp{16, i}{j,1}];
        SGsoleaire=[SGsoleaire, gaucheSyrineInterp{10, i}{j,1}];
        SGtendineux=[SGtendineux, gaucheSyrineInterp{12, i}{j,1}];
        SGtibialAnt=[SGtibialAnt, gaucheSyrineInterp{8, i}{j,1}];
        SGvaste=[SGvaste, gaucheSyrineInterp{6, i}{j,1}];
    end
end

SGbicepsMoy=mean(SGbiceps,2);
SGdroitAntMoy=mean(SGdroitAnt,2);
SGgastrocnemienMoy=mean(SGgastrocnemien,2);
SGglutealMoy=mean(SGgluteal,2);
SGsoleaireMoy=mean(SGsoleaire,2);
SGtendineuxMoy=mean(SGtendineux,2);
SGtibialAntMoy=mean(SGtibialAnt,2);
SGvasteMoy=mean(SGvaste,2);


SGbicepsSd=std(SGbiceps,1,2);
SGdroitAntSd=std(SGdroitAnt,1,2);
SGgastrocnemienSd=std(SGgastrocnemien,1,2);
SGglutealSd=std(SGgluteal,1,2);
SGsoleaireSd=std(SGsoleaire,1,2);
SGtendineuxSd=std(SGtendineux,1,2);
SGtibialAntSd=std(SGtibialAnt,1,2);
SGvasteSd=std(SGvaste,1,2);

%% comparaison biceps
x=0:1:100;
figure;
subplot(2,2,1); 
plot(x,LDbicepsMoy, 'k', 'LineWidth', 3);
hold on
plot(x,LDbicepsMoy-LDbicepsSd, 'k', 'LineWidth', 0.5);
hold on
plot(x,LDbicepsMoy+LDbicepsSd, 'k', 'LineWidth', 0.5);

title('Biceps Droit Lina')


subplot(2,2,2); 
plot(x,SDbicepsMoy, 'k', 'LineWidth', 3);
hold on
plot(x,SDbicepsMoy-SDbicepsSd, 'k', 'LineWidth', 0.5);
hold on
plot(x,SDbicepsMoy+SDbicepsSd, 'k', 'LineWidth', 0.5);
title('Biceps Droit Syrine')


subplot(2,2,3); 
plot(x,LGbicepsMoy, 'k', 'LineWidth', 3);
hold on
plot(x,LGbicepsMoy-LGbicepsSd, 'k', 'LineWidth', 0.5);
hold on
plot(x,LGbicepsMoy+LGbicepsSd, 'k', 'LineWidth', 0.5);
title('Biceps Gauche Lina')

subplot(2,2,4); 

plot(x,SGbicepsMoy, 'k', 'LineWidth', 3);
hold on
plot(x,SGbicepsMoy-SGbicepsSd, 'k', 'LineWidth', 0.5);
hold on
plot(x,SGbicepsMoy+SGbicepsSd, 'k', 'LineWidth', 0.5);
title('Biceps Gauche Syrine')


%% comparaison droit antérieur

figure;
subplot(2,2,1); 
plot(x,LDdroitAntMoy, 'k', 'LineWidth', 3);
hold on
plot(x,LDdroitAntMoy-LDdroitAntSd, 'k', 'LineWidth', 0.5);
hold on
plot(x,LDdroitAntMoy+LDdroitAntSd, 'k', 'LineWidth', 0.5);

title('Droit Antérieur Droit Lina')


subplot(2,2,2); 
plot(x,SDdroitAntMoy, 'k', 'LineWidth', 3);
hold on
plot(x,SDdroitAntMoy-SDdroitAntSd, 'k', 'LineWidth', 0.5);
hold on
plot(x,SDdroitAntMoy+SDdroitAntSd, 'k', 'LineWidth', 0.5);
title('Droit Anterieur Droit Syrine')


subplot(2,2,3); 
plot(x,LGdroitAntMoy, 'k', 'LineWidth', 3);
hold on
plot(x,LGdroitAntMoy-LGdroitAntSd, 'k', 'LineWidth', 0.5);
hold on
plot(x,LGdroitAntMoy+LGdroitAntSd, 'k', 'LineWidth', 0.5);
title('Droit Antérieur Gauche Lina')

subplot(2,2,4); 

plot(x,SGdroitAntMoy, 'k', 'LineWidth', 3);
hold on
plot(x,SGdroitAntMoy-SGdroitAntSd, 'k', 'LineWidth', 0.5);
hold on
plot(x,SGdroitAntMoy+SGdroitAntSd, 'k', 'LineWidth', 0.5);
title('Droit Antérieur Gauche Syrine')

%% comparaison gastrocnemien

figure;
subplot(2,2,1); 
plot(x,LDgastrocnemienMoy, 'k', 'LineWidth', 3);
hold on
plot(x,LDgastrocnemienMoy-LDgastrocnemienSd, 'k', 'LineWidth', 0.5);
hold on
plot(x,LDgastrocnemienMoy+LDgastrocnemienSd, 'k', 'LineWidth', 0.5);

title('Gastrocnemien Droit Lina')


subplot(2,2,2); 
plot(x,SDgastrocnemienMoy, 'k', 'LineWidth', 3);
hold on
plot(x,SDgastrocnemienMoy-SDgastrocnemienSd, 'k', 'LineWidth', 0.5);
hold on
plot(x,SDgastrocnemienMoy+SDgastrocnemienSd, 'k', 'LineWidth', 0.5);
title('Gastrocnemien Droit Syrine')


subplot(2,2,3); 
plot(x,LGgastrocnemienMoy, 'k', 'LineWidth', 3);
hold on
plot(x,LGgastrocnemienMoy-LGgastrocnemienSd, 'k', 'LineWidth', 0.5);
hold on
plot(x,LGgastrocnemienMoy+LGgastrocnemienSd, 'k', 'LineWidth', 0.5);
title('Gastrocnemien Gauche Lina')

subplot(2,2,4); 

plot(x,SGgastrocnemienMoy, 'k', 'LineWidth', 3);
hold on
plot(x,SGgastrocnemienMoy-SGgastrocnemienSd, 'k', 'LineWidth', 0.5);
hold on
plot(x,SGgastrocnemienMoy+SGgastrocnemienSd, 'k', 'LineWidth', 0.5);
title('Gastrocnemien Gauche Syrine')

%% comparaison gluteal

figure;
subplot(2,2,1); 
plot(x,LDglutealMoy, 'k', 'LineWidth', 3);
hold on
plot(x,LDglutealMoy-LDglutealSd, 'k', 'LineWidth', 0.5);
hold on
plot(x,LDglutealMoy+LDglutealSd, 'k', 'LineWidth', 0.5);

title('Gluteal Droit Lina')


subplot(2,2,2); 
plot(x,SDglutealMoy, 'k', 'LineWidth', 3);
hold on
plot(x,SDglutealMoy-SDglutealSd, 'k', 'LineWidth', 0.5);
hold on
plot(x,SDglutealMoy+SDglutealSd, 'k', 'LineWidth', 0.5);
title('Gluteal Droit Syrine')


subplot(2,2,3); 
plot(x,LGglutealMoy, 'k', 'LineWidth', 3);
hold on
plot(x,LGglutealMoy-LGglutealSd, 'k', 'LineWidth', 0.5);
hold on
plot(x,LGglutealMoy+LGglutealSd, 'k', 'LineWidth', 0.5);
title('Gluteal Gauche Lina')

subplot(2,2,4); 

plot(x,SGglutealMoy, 'k', 'LineWidth', 3);
hold on
plot(x,SGglutealMoy-SGglutealSd, 'k', 'LineWidth', 0.5);
hold on
plot(x,SGglutealMoy+SGglutealSd, 'k', 'LineWidth', 0.5);
title('Gluteal Gauche Syrine')


%% comparaison soleaire

figure;
subplot(2,2,1); 
plot(x,LDsoleaireMoy, 'k', 'LineWidth', 3);
hold on
plot(x,LDsoleaireMoy-LDsoleaireSd, 'k', 'LineWidth', 0.5);
hold on
plot(x,LDsoleaireMoy+LDsoleaireSd, 'k', 'LineWidth', 0.5);

title('Soleaire Droit Lina')


subplot(2,2,2); 
plot(x,SDsoleaireMoy, 'k', 'LineWidth', 3);
hold on
plot(x,SDsoleaireMoy-SDsoleaireSd, 'k', 'LineWidth', 0.5);
hold on
plot(x,SDsoleaireMoy+SDsoleaireSd, 'k', 'LineWidth', 0.5);
title('Soleaire Droit Syrine')


subplot(2,2,3); 
plot(x,LGsoleaireMoy, 'k', 'LineWidth', 3);
hold on
plot(x,LGsoleaireMoy-LGsoleaireSd, 'k', 'LineWidth', 0.5);
hold on
plot(x,LGsoleaireMoy+LGsoleaireSd, 'k', 'LineWidth', 0.5);
title('Soleaire Gauche Lina')

subplot(2,2,4); 
plot(x,SGsoleaireMoy, 'k', 'LineWidth', 3);
hold on
plot(x,SGsoleaireMoy-SGsoleaireSd, 'k', 'LineWidth', 0.5);
hold on
plot(x,SGsoleaireMoy+SGsoleaireSd, 'k', 'LineWidth', 0.5);
title('Soleaire Gauche Syrine')


%% comparaison tendineux

figure;
subplot(2,2,1); 
plot(x,LDtendineuxMoy, 'k', 'LineWidth', 3);
hold on
plot(x,LDtendineuxMoy-LDtendineuxSd, 'k', 'LineWidth', 0.5);
hold on
plot(x,LDtendineuxMoy+LDtendineuxSd, 'k', 'LineWidth', 0.5);

title('Tendineux Droit Lina')


subplot(2,2,2); 
plot(x,SDtendineuxMoy, 'k', 'LineWidth', 3);
hold on
plot(x,SDtendineuxMoy-SDtendineuxSd, 'k', 'LineWidth', 0.5);
hold on
plot(x,SDtendineuxMoy+SDtendineuxSd, 'k', 'LineWidth', 0.5);
title('Tendineux Droit Syrine')


subplot(2,2,3); 
plot(x,LGtendineuxMoy, 'k', 'LineWidth', 3);
hold on
plot(x,LGtendineuxMoy-LGtendineuxSd, 'k', 'LineWidth', 0.5);
hold on
plot(x,LGtendineuxMoy+LGtendineuxSd, 'k', 'LineWidth', 0.5);
title('Tendineux Gauche Lina')

subplot(2,2,4); 
plot(x,SGtendineuxMoy, 'k', 'LineWidth', 3);
hold on
plot(x,SGtendineuxMoy-SGtendineuxSd, 'k', 'LineWidth', 0.5);
hold on
plot(x,SGtendineuxMoy+SGtendineuxSd, 'k', 'LineWidth', 0.5);
title('Tendineux Gauche Syrine')


%% comparaison tibial anterieur

figure;
subplot(2,2,1); 
plot(x,LDtibialAntMoy, 'k', 'LineWidth', 3);
hold on
plot(x,LDtibialAntMoy-LDtibialAntSd, 'k', 'LineWidth', 0.5);
hold on
plot(x,LDtibialAntMoy+LDtibialAntSd, 'k', 'LineWidth', 0.5);

title('Tibial Anterieur Droit Lina')


subplot(2,2,2); 
plot(x,SDtibialAntMoy, 'k', 'LineWidth', 3);
hold on
plot(x,SDtibialAntMoy-SDtibialAntSd, 'k', 'LineWidth', 0.5);
hold on
plot(x,SDtibialAntMoy+SDtibialAntSd, 'k', 'LineWidth', 0.5);
title('Tibial Anterieur Droit Syrine')


subplot(2,2,3); 
plot(x,LGtibialAntMoy, 'k', 'LineWidth', 3);
hold on
plot(x,LGtibialAntMoy-LGtibialAntSd, 'k', 'LineWidth', 0.5);
hold on
plot(x,LGtibialAntMoy+LGtibialAntSd, 'k', 'LineWidth', 0.5);
title('Tibial Anterieur Gauche Lina')

subplot(2,2,4); 
plot(x,SGtibialAntMoy, 'k', 'LineWidth', 3);
hold on
plot(x,SGtibialAntMoy-SGtibialAntSd, 'k', 'LineWidth', 0.5);
hold on
plot(x,SGtibialAntMoy+SGtibialAntSd, 'k', 'LineWidth', 0.5);
title('Tibial Anterieur Gauche Syrine')


%% comparaison vaste

figure;
subplot(2,2,1); 
plot(x,LDvasteMoy, 'k', 'LineWidth', 3);
hold on
plot(x,LDvasteMoy-LDvasteSd, 'k', 'LineWidth', 0.5);
hold on
plot(x,LDvasteMoy+LDvasteSd, 'k', 'LineWidth', 0.5);

title('Vaste Droit Lina')


subplot(2,2,2); 
plot(x,SDvasteMoy, 'k', 'LineWidth', 3);
hold on
plot(x,SDvasteMoy-SDvasteSd, 'k', 'LineWidth', 0.5);
hold on
plot(x,SDvasteMoy+SDvasteSd, 'k', 'LineWidth', 0.5);
title('Vaste Droit Syrine')


subplot(2,2,3); 
plot(x,LGvasteMoy, 'k', 'LineWidth', 3);
hold on
plot(x,LGvasteMoy-LGvasteSd, 'k', 'LineWidth', 0.5);
hold on
plot(x,LGvasteMoy+LGvasteSd, 'k', 'LineWidth', 0.5);
title('Vaste Gauche Lina')

subplot(2,2,4); 
plot(x,SGvasteMoy, 'k', 'LineWidth', 3);
hold on
plot(x,SGvasteMoy-SGvasteSd, 'k', 'LineWidth', 0.5);
hold on
plot(x,SGvasteMoy+SGvasteSd, 'k', 'LineWidth', 0.5);
title('Vaste Gauche Syrine')

