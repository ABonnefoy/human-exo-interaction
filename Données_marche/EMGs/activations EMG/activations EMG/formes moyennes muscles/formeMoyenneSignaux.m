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
        

%% on normalise par rapport à l'amplitude et on fait la moyenne pour chaque muscle
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

%%

LDbicepsNormalized=zeros(101,50);
for i=1:77
    signal=LDbiceps(:,i);
    moy=mean(signal);
    for j=1:101
        LDbicepsNormalized(j,i)=100*signal(j)/moy;
    end
end

LDdroitAntNormalized=zeros(101,50);
for i=1:77
    signal=LDdroitAnt(:,i);
    moy=mean(signal);
    for j=1:101
        LDdroitAntNormalized(j,i)=100*signal(j)/moy;
    end
end

LDgastrocnemienNormalized=zeros(101,50);
for i=1:77
    signal=LDgastrocnemien(:,i);
    moy=mean(signal);
    for j=1:101
        LDgastrocnemienNormalized(j,i)=100*signal(j)/moy;
    end
end

LDglutealNormalized=zeros(101,50);
for i=1:77
    signal=LDgluteal(:,i);
    moy=mean(signal);
    for j=1:101
        LDglutealNormalized(j,i)=100*signal(j)/moy;
    end
end

LDsoleaireNormalized=zeros(101,50);
for i=1:77
    signal=LDsoleaire(:,i);
    moy=mean(signal);
    for j=1:101
        LDsoleaireNormalized(j,i)=100*signal(j)/moy;
    end
end

LDtendineuxNormalized=zeros(101,50);
for i=1:77
    signal=LDtendineux(:,i);
    moy=mean(signal);
    for j=1:101
        LDtendineuxNormalized(j,i)=100*signal(j)/moy;
    end
end

LDtibialAntNormalized=zeros(101,50);
for i=1:77
    signal=LDtibialAnt(:,i);
    moy=mean(signal);
    for j=1:101
        LDtibialAntNormalized(j,i)=100*signal(j)/moy;
    end
end

LDvasteNormalized=zeros(101,50);
for i=1:77
    signal=LDvaste(:,i);
    moy=mean(signal);
    for j=1:101
        LDvasteNormalized(j,i)=100*signal(j)/moy;
    end
end


%%

LDbicepsMoy=mean(LDbicepsNormalized,2);
LDdroitAntMoy=mean(LDdroitAntNormalized,2);
LDgastrocnemienMoy=mean(LDgastrocnemienNormalized,2);
LDglutealMoy=mean(LDglutealNormalized,2);
LDsoleaireMoy=mean(LDsoleaireNormalized,2);
LDtendineuxMoy=mean(LDtendineuxNormalized,2);
LDtibialAntMoy=mean(LDtibialAntNormalized,2);
LDvasteMoy=mean(LDvasteNormalized,2);


LDbicepsSd=std(LDbicepsNormalized,1,2);
LDdroitAntSd=std(LDdroitAntNormalized,1,2);
LDgastrocnemienSd=std(LDgastrocnemienNormalized,1,2);
LDglutealSd=std(LDglutealNormalized,1,2);
LDsoleaireSd=std(LDsoleaireNormalized,1,2);
LDtendineuxSd=std(LDtendineuxNormalized,1,2);
LDtibialAntSd=std(LDtibialAntNormalized,1,2);
LDvasteSd=std(LDvasteNormalized,1,2);


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

LGbicepsNormalized=zeros(101,50);
for i=1:77
    signal=LGbiceps(:,i);
    moy=mean(signal);
    for j=1:101
        LGbicepsNormalized(j,i)=100*signal(j)/moy;
    end
end

LGdroitAntNormalized=zeros(101,50);
for i=1:77
    signal=LGdroitAnt(:,i);
    moy=mean(signal);
    for j=1:101
        LGdroitAntNormalized(j,i)=100*signal(j)/moy;
    end
end

LGgastrocnemienNormalized=zeros(101,50);
for i=1:77
    signal=LGgastrocnemien(:,i);
    moy=mean(signal);
    for j=1:101
        LGgastrocnemienNormalized(j,i)=100*signal(j)/moy;
    end
end

LGglutealNormalized=zeros(101,50);
for i=1:77
    signal=LGgluteal(:,i);
    moy=mean(signal);
    for j=1:101
        LGglutealNormalized(j,i)=100*signal(j)/moy;
    end
end

LGsoleaireNormalized=zeros(101,50);
for i=1:77
    signal=LGsoleaire(:,i);
    moy=mean(signal);
    for j=1:101
        LGsoleaireNormalized(j,i)=100*signal(j)/moy;
    end
end

LGtendineuxNormalized=zeros(101,50);
for i=1:77
    signal=LGtendineux(:,i);
    moy=mean(signal);
    for j=1:101
        LGtendineuxNormalized(j,i)=100*signal(j)/moy;
    end
end

LGtibialAntNormalized=zeros(101,50);
for i=1:77
    signal=LGtibialAnt(:,i);
    moy=mean(signal);
    for j=1:101
        LGtibialAntNormalized(j,i)=100*signal(j)/moy;
    end
end

LGvasteNormalized=zeros(101,50);
for i=1:77
    signal=LGvaste(:,i);
    moy=mean(signal);
    for j=1:101
        LGvasteNormalized(j,i)=100*signal(j)/moy;
    end
end

LGbicepsMoy=mean(LGbicepsNormalized,2);
LGdroitAntMoy=mean(LGdroitAntNormalized,2);
LGgastrocnemienMoy=mean(LGgastrocnemienNormalized,2);
LGglutealMoy=mean(LGglutealNormalized,2);
LGsoleaireMoy=mean(LGsoleaireNormalized,2);
LGtendineuxMoy=mean(LGtendineuxNormalized,2);
LGtibialAntMoy=mean(LGtibialAntNormalized,2);
LGvasteMoy=mean(LGvasteNormalized,2);


LGbicepsSd=std(LGbicepsNormalized,1,2);
LGdroitAntSd=std(LGdroitAntNormalized,1,2);
LGgastrocnemienSd=std(LGgastrocnemienNormalized,1,2);
LGglutealSd=std(LGglutealNormalized,1,2);
LGsoleaireSd=std(LGsoleaireNormalized,1,2);
LGtendineuxSd=std(LGtendineuxNormalized,1,2);
LGtibialAntSd=std(LGtibialAntNormalized,1,2);
LGvasteSd=std(LGvasteNormalized,1,2);


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


SDbicepsNormalized=zeros(101,50);
for i=1:50
    signal=SDbiceps(:,i);
    moy=mean(signal);
    for j=1:101
        SDbicepsNormalized(j,i)=100*signal(j)/moy;
    end
end

SDdroitAntNormalized=zeros(101,50);
for i=1:50
    signal=SDdroitAnt(:,i);
    moy=mean(signal);
    for j=1:101
        SDdroitAntNormalized(j,i)=100*signal(j)/moy;
    end
end

SDgastrocnemienNormalized=zeros(101,50);
for i=1:50
    signal=SDgastrocnemien(:,i);
    moy=mean(signal);
    for j=1:101
        SDgastrocnemienNormalized(j,i)=100*signal(j)/moy;
    end
end

SDglutealNormalized=zeros(101,50);
for i=1:50
    signal=SDgluteal(:,i);
    moy=mean(signal);
    for j=1:101
        SDglutealNormalized(j,i)=100*signal(j)/moy;
    end
end

SDsoleaireNormalized=zeros(101,50);
for i=1:50
    signal=SDsoleaire(:,i);
    moy=mean(signal);
    for j=1:101
        SDsoleaireNormalized(j,i)=100*signal(j)/moy;
    end
end

SDtendineuxNormalized=zeros(101,50);
for i=1:50
    signal=SDtendineux(:,i);
    moy=mean(signal);
    for j=1:101
        SDtendineuxNormalized(j,i)=100*signal(j)/moy;
    end
end

SDtibialAntNormalized=zeros(101,50);
for i=1:50
    signal=SDtibialAnt(:,i);
    moy=mean(signal);
    for j=1:101
        SDtibialAntNormalized(j,i)=100*signal(j)/moy;
    end
end

SDvasteNormalized=zeros(101,50);
for i=1:50
    signal=SDvaste(:,i);
    moy=mean(signal);
    for j=1:101
        SDvasteNormalized(j,i)=100*signal(j)/moy;
    end
end





SDbicepsMoy=mean(SDbicepsNormalized,2);
SDdroitAntMoy=mean(SDdroitAntNormalized,2);
SDgastrocnemienMoy=mean(SDgastrocnemienNormalized,2);
SDglutealMoy=mean(SDglutealNormalized,2);
SDsoleaireMoy=mean(SDsoleaireNormalized,2);
SDtendineuxMoy=mean(SDtendineuxNormalized,2);
SDtibialAntMoy=mean(SDtibialAntNormalized,2);
SDvasteMoy=mean(SDvasteNormalized,2);


SDbicepsSd=std(SDbicepsNormalized,1,2);
SDdroitAntSd=std(SDdroitAntNormalized,1,2);
SDgastrocnemienSd=std(SDgastrocnemienNormalized,1,2);
SDglutealSd=std(SDglutealNormalized,1,2);
SDsoleaireSd=std(SDsoleaireNormalized,1,2);
SDtendineuxSd=std(SDtendineuxNormalized,1,2);
SDtibialAntSd=std(SDtibialAntNormalized,1,2);
SDvasteSd=std(SDvasteNormalized,1,2);


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


SGbicepsNormalized=zeros(101,50);
for i=1:50
    signal=SGbiceps(:,i);
    moy=mean(signal);
    for j=1:101
        SGbicepsNormalized(j,i)=100*signal(j)/moy;
    end
end

SGdroitAntNormalized=zeros(101,50);
for i=1:50
    signal=SGdroitAnt(:,i);
    moy=mean(signal);
    for j=1:101
        SGdroitAntNormalized(j,i)=100*signal(j)/moy;
    end
end

SGgastrocnemienNormalized=zeros(101,50);
for i=1:50
    signal=SGgastrocnemien(:,i);
    moy=mean(signal);
    for j=1:101
        SGgastrocnemienNormalized(j,i)=100*signal(j)/moy;
    end
end

SGglutealNormalized=zeros(101,50);
for i=1:50
    signal=SGgluteal(:,i);
    moy=mean(signal);
    for j=1:101
        SGglutealNormalized(j,i)=100*signal(j)/moy;
    end
end

SGsoleaireNormalized=zeros(101,50);
for i=1:50
    signal=SGsoleaire(:,i);
    moy=mean(signal);
    for j=1:101
        SGsoleaireNormalized(j,i)=100*signal(j)/moy;
    end
end

SGtendineuxNormalized=zeros(101,50);
for i=1:50
    signal=SGtendineux(:,i);
    moy=mean(signal);
    for j=1:101
        SGtendineuxNormalized(j,i)=100*signal(j)/moy;
    end
end

SGtibialAntNormalized=zeros(101,50);
for i=1:50
    signal=SGtibialAnt(:,i);
    moy=mean(signal);
    for j=1:101
        SGtibialAntNormalized(j,i)=100*signal(j)/moy;
    end
end

SGvasteNormalized=zeros(101,50);
for i=1:50
    signal=SGvaste(:,i);
    moy=mean(signal);
    for j=1:101
        SGvasteNormalized(j,i)=100*signal(j)/moy;
    end
end




SGbicepsMoy=mean(SGbicepsNormalized,2);
SGdroitAntMoy=mean(SGdroitAntNormalized,2);
SGgastrocnemienMoy=mean(SGgastrocnemienNormalized,2);
SGglutealMoy=mean(SGglutealNormalized,2);
SGsoleaireMoy=mean(SGsoleaireNormalized,2);
SGtendineuxMoy=mean(SGtendineuxNormalized,2);
SGtibialAntMoy=mean(SGtibialAntNormalized,2);
SGvasteMoy=mean(SGvasteNormalized,2);


SGbicepsSd=std(SGbicepsNormalized,1,2);
SGdroitAntSd=std(SGdroitAntNormalized,1,2);
SGgastrocnemienSd=std(SGgastrocnemienNormalized,1,2);
SGglutealSd=std(SGglutealNormalized,1,2);
SGsoleaireSd=std(SGsoleaireNormalized,1,2);
SGtendineuxSd=std(SGtendineuxNormalized,1,2);
SGtibialAntSd=std(SGtibialAntNormalized,1,2);
SGvasteSd=std(SGvasteNormalized,1,2);

%% comparaison biceps
x=0:1:100;
figure;
subplot(2,2,1); 
plot(x,LDbicepsMoy, 'k', 'LineWidth', 3);
hold on
plot(x,LDbicepsMoy-LDbicepsSd, 'k', 'LineWidth', 0.5);
hold on
plot(x,LDbicepsMoy+LDbicepsSd, 'k', 'LineWidth', 0.5);
xline(meanToeOffDroiteLina)
yline(100)

title('Biceps Droit Lina')


subplot(2,2,2); 
plot(x,SDbicepsMoy, 'k', 'LineWidth', 3);
hold on
plot(x,SDbicepsMoy-SDbicepsSd, 'k', 'LineWidth', 0.5);
hold on
plot(x,SDbicepsMoy+SDbicepsSd, 'k', 'LineWidth', 0.5);
xline(meanToeOffDroiteSyrine)
title('Biceps Droit Syrine')
yline(100)


subplot(2,2,3); 
plot(x,LGbicepsMoy, 'k', 'LineWidth', 3);
hold on
plot(x,LGbicepsMoy-LGbicepsSd, 'k', 'LineWidth', 0.5);
hold on
plot(x,LGbicepsMoy+LGbicepsSd, 'k', 'LineWidth', 0.5);
xline(meanToeOffGaucheLina)
title('Biceps Gauche Lina')
yline(100)

subplot(2,2,4); 

plot(x,SGbicepsMoy, 'k', 'LineWidth', 3);
hold on
plot(x,SGbicepsMoy-SGbicepsSd, 'k', 'LineWidth', 0.5);
hold on
plot(x,SGbicepsMoy+SGbicepsSd, 'k', 'LineWidth', 0.5);
xline(meanToeOffGaucheSyrine)
title('Biceps Gauche Syrine')
yline(100)


%% comparaison droit antérieur

figure;
subplot(2,2,1); 
plot(x,LDdroitAntMoy, 'k', 'LineWidth', 3);
hold on
plot(x,LDdroitAntMoy-LDdroitAntSd, 'k', 'LineWidth', 0.5);
hold on
plot(x,LDdroitAntMoy+LDdroitAntSd, 'k', 'LineWidth', 0.5);
xline(meanToeOffDroiteLina)
title('Droit Antérieur Droit Lina')
yline(100)


subplot(2,2,2); 
plot(x,SDdroitAntMoy, 'k', 'LineWidth', 3);
hold on
plot(x,SDdroitAntMoy-SDdroitAntSd, 'k', 'LineWidth', 0.5);
hold on
plot(x,SDdroitAntMoy+SDdroitAntSd, 'k', 'LineWidth', 0.5);
xline(meanToeOffDroiteSyrine)
title('Droit Anterieur Droit Syrine')
yline(100)


subplot(2,2,3); 
plot(x,LGdroitAntMoy, 'k', 'LineWidth', 3);
hold on
plot(x,LGdroitAntMoy-LGdroitAntSd, 'k', 'LineWidth', 0.5);
hold on
plot(x,LGdroitAntMoy+LGdroitAntSd, 'k', 'LineWidth', 0.5);
xline(meanToeOffGaucheLina)
title('Droit Antérieur Gauche Lina')
yline(100)

subplot(2,2,4); 

plot(x,SGdroitAntMoy, 'k', 'LineWidth', 3);
hold on
plot(x,SGdroitAntMoy-SGdroitAntSd, 'k', 'LineWidth', 0.5);
hold on
plot(x,SGdroitAntMoy+SGdroitAntSd, 'k', 'LineWidth', 0.5);
xline(meanToeOffGaucheSyrine)
title('Droit Antérieur Gauche Syrine')
yline(100)

%% comparaison gastrocnemien

figure;
subplot(2,2,1); 
plot(x,LDgastrocnemienMoy, 'k', 'LineWidth', 3);
hold on
plot(x,LDgastrocnemienMoy-LDgastrocnemienSd, 'k', 'LineWidth', 0.5);
hold on
plot(x,LDgastrocnemienMoy+LDgastrocnemienSd, 'k', 'LineWidth', 0.5);
xline(meanToeOffDroiteLina)
title('Gastrocnemien Droit Lina')
yline(100)


subplot(2,2,2); 
plot(x,SDgastrocnemienMoy, 'k', 'LineWidth', 3);
hold on
plot(x,SDgastrocnemienMoy-SDgastrocnemienSd, 'k', 'LineWidth', 0.5);
hold on
plot(x,SDgastrocnemienMoy+SDgastrocnemienSd, 'k', 'LineWidth', 0.5);
xline(meanToeOffDroiteSyrine)
title('Gastrocnemien Droit Syrine')
yline(100)


subplot(2,2,3); 
plot(x,LGgastrocnemienMoy, 'k', 'LineWidth', 3);
hold on
plot(x,LGgastrocnemienMoy-LGgastrocnemienSd, 'k', 'LineWidth', 0.5);
hold on
plot(x,LGgastrocnemienMoy+LGgastrocnemienSd, 'k', 'LineWidth', 0.5);
xline(meanToeOffGaucheLina)
title('Gastrocnemien Gauche Lina')
yline(100)

subplot(2,2,4); 

plot(x,SGgastrocnemienMoy, 'k', 'LineWidth', 3);
hold on
plot(x,SGgastrocnemienMoy-SGgastrocnemienSd, 'k', 'LineWidth', 0.5);
hold on
plot(x,SGgastrocnemienMoy+SGgastrocnemienSd, 'k', 'LineWidth', 0.5);
xline(meanToeOffGaucheSyrine)
title('Gastrocnemien Gauche Syrine')
yline(100)

%% comparaison gluteal

figure;
subplot(2,2,1); 
plot(x,LDglutealMoy, 'k', 'LineWidth', 3);
hold on
plot(x,LDglutealMoy-LDglutealSd, 'k', 'LineWidth', 0.5);
hold on
plot(x,LDglutealMoy+LDglutealSd, 'k', 'LineWidth', 0.5);
xline(meanToeOffDroiteLina)
title('Gluteal Droit Lina')
yline(100)


subplot(2,2,2); 
plot(x,SDglutealMoy, 'k', 'LineWidth', 3);
hold on
plot(x,SDglutealMoy-SDglutealSd, 'k', 'LineWidth', 0.5);
hold on
plot(x,SDglutealMoy+SDglutealSd, 'k', 'LineWidth', 0.5);
xline(meanToeOffDroiteSyrine)
title('Gluteal Droit Syrine')
yline(100)


subplot(2,2,3); 
plot(x,LGglutealMoy, 'k', 'LineWidth', 3);
hold on
plot(x,LGglutealMoy-LGglutealSd, 'k', 'LineWidth', 0.5);
hold on
plot(x,LGglutealMoy+LGglutealSd, 'k', 'LineWidth', 0.5);
xline(meanToeOffGaucheLina)
title('Gluteal Gauche Lina')
yline(100)

subplot(2,2,4); 

plot(x,SGglutealMoy, 'k', 'LineWidth', 3);
hold on
plot(x,SGglutealMoy-SGglutealSd, 'k', 'LineWidth', 0.5);
hold on
plot(x,SGglutealMoy+SGglutealSd, 'k', 'LineWidth', 0.5);
xline(meanToeOffGaucheSyrine)
title('Gluteal Gauche Syrine')
yline(100)


%% comparaison soleaire

figure;
subplot(2,2,1); 
plot(x,LDsoleaireMoy, 'k', 'LineWidth', 3);
hold on
plot(x,LDsoleaireMoy-LDsoleaireSd, 'k', 'LineWidth', 0.5);
hold on
plot(x,LDsoleaireMoy+LDsoleaireSd, 'k', 'LineWidth', 0.5);
xline(meanToeOffDroiteLina)
title('Soleaire Droit Lina')
yline(100)


subplot(2,2,2); 
plot(x,SDsoleaireMoy, 'k', 'LineWidth', 3);
hold on
plot(x,SDsoleaireMoy-SDsoleaireSd, 'k', 'LineWidth', 0.5);
hold on
plot(x,SDsoleaireMoy+SDsoleaireSd, 'k', 'LineWidth', 0.5);
xline(meanToeOffDroiteSyrine)
title('Soleaire Droit Syrine')
yline(100)

subplot(2,2,3); 
plot(x,LGsoleaireMoy, 'k', 'LineWidth', 3);
hold on
plot(x,LGsoleaireMoy-LGsoleaireSd, 'k', 'LineWidth', 0.5);
hold on
plot(x,LGsoleaireMoy+LGsoleaireSd, 'k', 'LineWidth', 0.5);
xline(meanToeOffGaucheLina)
title('Soleaire Gauche Lina')
yline(100)

subplot(2,2,4); 
plot(x,SGsoleaireMoy, 'k', 'LineWidth', 3);
hold on
plot(x,SGsoleaireMoy-SGsoleaireSd, 'k', 'LineWidth', 0.5);
hold on
plot(x,SGsoleaireMoy+SGsoleaireSd, 'k', 'LineWidth', 0.5);
xline(meanToeOffGaucheSyrine)
title('Soleaire Gauche Syrine')
yline(100)


%% comparaison tendineux

figure;
subplot(2,2,1); 
plot(x,LDtendineuxMoy, 'k', 'LineWidth', 3);
hold on
plot(x,LDtendineuxMoy-LDtendineuxSd, 'k', 'LineWidth', 0.5);
hold on
plot(x,LDtendineuxMoy+LDtendineuxSd, 'k', 'LineWidth', 0.5);
xline(meanToeOffDroiteLina)
title('Tendineux Droit Lina')
yline(100)


subplot(2,2,2); 
plot(x,SDtendineuxMoy, 'k', 'LineWidth', 3);
hold on
plot(x,SDtendineuxMoy-SDtendineuxSd, 'k', 'LineWidth', 0.5);
hold on
plot(x,SDtendineuxMoy+SDtendineuxSd, 'k', 'LineWidth', 0.5);
xline(meanToeOffDroiteSyrine)
title('Tendineux Droit Syrine')
yline(100)


subplot(2,2,3); 
plot(x,LGtendineuxMoy, 'k', 'LineWidth', 3);
hold on
plot(x,LGtendineuxMoy-LGtendineuxSd, 'k', 'LineWidth', 0.5);
hold on
plot(x,LGtendineuxMoy+LGtendineuxSd, 'k', 'LineWidth', 0.5);
xline(meanToeOffGaucheLina)
title('Tendineux Gauche Lina')
yline(100)

subplot(2,2,4); 
plot(x,SGtendineuxMoy, 'k', 'LineWidth', 3);
hold on
plot(x,SGtendineuxMoy-SGtendineuxSd, 'k', 'LineWidth', 0.5);
hold on
plot(x,SGtendineuxMoy+SGtendineuxSd, 'k', 'LineWidth', 0.5);
xline(meanToeOffGaucheSyrine)
title('Tendineux Gauche Syrine')
yline(100)


%% comparaison tibial anterieur

figure;
subplot(2,2,1); 
plot(x,LDtibialAntMoy, 'k', 'LineWidth', 3);
hold on
plot(x,LDtibialAntMoy-LDtibialAntSd, 'k', 'LineWidth', 0.5);
hold on
plot(x,LDtibialAntMoy+LDtibialAntSd, 'k', 'LineWidth', 0.5);
xline(meanToeOffDroiteLina)
title('Tibial Anterieur Droit Lina')
yline(100)

subplot(2,2,2); 
plot(x,SDtibialAntMoy, 'k', 'LineWidth', 3);
hold on
plot(x,SDtibialAntMoy-SDtibialAntSd, 'k', 'LineWidth', 0.5);
hold on
plot(x,SDtibialAntMoy+SDtibialAntSd, 'k', 'LineWidth', 0.5);
xline(meanToeOffDroiteSyrine)
title('Tibial Anterieur Droit Syrine')
yline(100)

subplot(2,2,3); 
plot(x,LGtibialAntMoy, 'k', 'LineWidth', 3);
hold on
plot(x,LGtibialAntMoy-LGtibialAntSd, 'k', 'LineWidth', 0.5);
hold on
plot(x,LGtibialAntMoy+LGtibialAntSd, 'k', 'LineWidth', 0.5);
xline(meanToeOffGaucheLina)
title('Tibial Anterieur Gauche Lina')
yline(100)

subplot(2,2,4); 
plot(x,SGtibialAntMoy, 'k', 'LineWidth', 3);
hold on
plot(x,SGtibialAntMoy-SGtibialAntSd, 'k', 'LineWidth', 0.5);
hold on
plot(x,SGtibialAntMoy+SGtibialAntSd, 'k', 'LineWidth', 0.5);
xline(meanToeOffGaucheSyrine)
title('Tibial Anterieur Gauche Syrine')
yline(100)

%% comparaison vaste

figure;
subplot(2,2,1); 
plot(x,LDvasteMoy, 'k', 'LineWidth', 3);
hold on
plot(x,LDvasteMoy-LDvasteSd, 'k', 'LineWidth', 0.5);
hold on
plot(x,LDvasteMoy+LDvasteSd, 'k', 'LineWidth', 0.5);
xline(meanToeOffDroiteLina)
title('Vaste Droit Lina')
yline(100)

subplot(2,2,2); 
plot(x,SDvasteMoy, 'k', 'LineWidth', 3);
hold on
plot(x,SDvasteMoy-SDvasteSd, 'k', 'LineWidth', 0.5);
hold on
plot(x,SDvasteMoy+SDvasteSd, 'k', 'LineWidth', 0.5);
xline(meanToeOffDroiteSyrine)
title('Vaste Droit Syrine')
yline(100)

subplot(2,2,3); 
plot(x,LGvasteMoy, 'k', 'LineWidth', 3);
hold on
plot(x,LGvasteMoy-LGvasteSd, 'k', 'LineWidth', 0.5);
hold on
plot(x,LGvasteMoy+LGvasteSd, 'k', 'LineWidth', 0.5);
xline(meanToeOffGaucheLina)
title('Vaste Gauche Lina')
yline(100)

subplot(2,2,4); 
plot(x,SGvasteMoy, 'k', 'LineWidth', 3);
hold on
plot(x,SGvasteMoy-SGvasteSd, 'k', 'LineWidth', 0.5);
hold on
plot(x,SGvasteMoy+SGvasteSd, 'k', 'LineWidth', 0.5);
xline(meanToeOffGaucheSyrine)
title('Vaste Gauche Syrine')
yline(100)



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Analyse statistique
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% On commence par tester la normalité des echantillons

% Pour Lina, on sélectionne uniquement 50 parmi les 77 cycles disponibles 

r=randsample(77,50);

LDbicepsNormalized=LDbicepsNormalized(:,r);
LDdroitAntNormalized=LDdroitAntNormalized(:,r);
LDgastrocnemienNormalized=LDgastrocnemienNormalized(:,r);
LDglutealNormalized=LDglutealNormalized(:,r);
LDsoleaireNormalized=LDsoleaireNormalized(:,r);
LDtendineuxNormalized=LDtendineuxNormalized(:,r);
LDtibialAntNormalized=LDtibialAntNormalized(:,r);
LDvasteNormalized=LDvasteNormalized(:,r);

LGbicepsNormalized=LGbicepsNormalized(:,r);
LGdroitAntNormalized=LGdroitAntNormalized(:,r);
LGgastrocnemienNormalized=LGgastrocnemienNormalized(:,r);
LGglutealNormalized=LGglutealNormalized(:,r);
LGsoleaireNormalized=LGsoleaireNormalized(:,r);
LGtendineuxNormalized=LGtendineuxNormalized(:,r);
LGtibialAntNormalized=LGtibialAntNormalized(:,r);
LGvasteNormalized=LGvasteNormalized(:,r);
%%
LDNormalized=cell(8,1);
LDNormalized{1,1}=LDbicepsNormalized;
LDNormalized{2,1}=LDdroitAntNormalized;
LDNormalized{3,1}=LDgastrocnemienNormalized;
LDNormalized{4,1}=LDglutealNormalized;
LDNormalized{5,1}=LDsoleaireNormalized;
LDNormalized{6,1}=LDtendineuxNormalized;
LDNormalized{7,1}=LDtibialAntNormalized;
LDNormalized{8,1}=LDvasteNormalized;

LGNormalized=cell(8,1);
LGNormalized{1,1}=LGbicepsNormalized;
LGNormalized{2,1}=LGdroitAntNormalized;
LGNormalized{3,1}=LGgastrocnemienNormalized;
LGNormalized{4,1}=LGglutealNormalized;
LGNormalized{5,1}=LGsoleaireNormalized;
LGNormalized{6,1}=LGtendineuxNormalized;
LGNormalized{7,1}=LGtibialAntNormalized;
LGNormalized{8,1}=LGvasteNormalized;

SDNormalized=cell(8,1);
SDNormalized{1,1}=SDbicepsNormalized;
SDNormalized{2,1}=SDdroitAntNormalized;
SDNormalized{3,1}=SDgastrocnemienNormalized;
SDNormalized{4,1}=SDglutealNormalized;
SDNormalized{5,1}=SDsoleaireNormalized;
SDNormalized{6,1}=SDtendineuxNormalized;
SDNormalized{7,1}=SDtibialAntNormalized;
SDNormalized{8,1}=SDvasteNormalized;

SGNormalized=cell(8,1);
SGNormalized{1,1}=SGbicepsNormalized;
SGNormalized{2,1}=SGdroitAntNormalized;
SGNormalized{3,1}=SGgastrocnemienNormalized;
SGNormalized{4,1}=SGglutealNormalized;
SGNormalized{5,1}=SGsoleaireNormalized;
SGNormalized{6,1}=SGtendineuxNormalized;
SGNormalized{7,1}=SGtibialAntNormalized;
SGNormalized{8,1}=SGvasteNormalized;

% Test de normalité : Shapiro Wilk
pvalNormLD=zeros(101, 8); % ordre alphabétique
pvalNormLG=zeros(101, 8);

pvalNormSD=zeros(101,8);
pvalNormSG=zeros(101,8);
%%
for i=1:8
    for j=1:101
        l=LDNormalized{i,1}(j,:);
        results=ShapiroWilk(l);
        pvNormLD(j,i)=results(2);
        
        l=SDNormalized{i,1}(j,:);
        results=ShapiroWilk(l);
        pvNormSD(j,i)=results(2);
    end
end

% p valeurs > 0.05 dans +50% des valeurs : on ne peut pas faire un test de
% Student parce qu'on a même pas la normalité des échantillons


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% On utilise donc Wilcoxon rank sum test

droitepval=zeros(101, 8); % ordre alphabétique
gauchepval=zeros(101, 8);

% comparaison du même muscle du même côté chez les filles
for i=1:8
    for j=1:101
        linaDroite=LDNormalized{i,1}(j,:);
        syrineDroite=SDNormalized{i,1}(j,:);
        [p,h]=ranksum(linaDroite, syrineDroite);
        droitepval(j,i)=p;
        
        linaGauche=LGNormalized{i,1}(j,:);
        syrineGauche=SGNormalized{i,1}(j,:);
        [p,h]=ranksum(linaGauche, syrineGauche);
        gauchepval(j,i)=p;
    end
end

linapval=zeros(101, 8); % ordre alphabétique
linapval=zeros(101, 8);

% comparaison du même muscle du même côté chez les filles
for i=1:8
    for j=1:101
        linaDroite=LDNormalized{i,1}(j,:);
        linaGauche=LGNormalized{i,1}(j,:);
        [p,h]=ranksum(linaDroite, linaGauche);
        linapval(j,i)=p;
        
        syrineDroite=SDNormalized{i,1}(j,:);
        syrineGauche=SGNormalized{i,1}(j,:);
        [p,h]=ranksum(syrineDroite, syrineGauche);
        syrinepval(j,i)=p;
    end
end
        
        
%%

figure;
subplot(2,1,1)
p1=plot(x,LDvasteMoy, 'r', 'LineWidth', 2);
hold on
plot(x,LDvasteMoy-LDvasteSd, 'r', 'LineWidth', 0.5);
hold on
plot(x,LDvasteMoy+LDvasteSd, 'r', 'LineWidth', 0.5);
xline(meanToeOffDroiteLina, 'r')
title('Vaste Droit ')
yline(100)
        
hold on
p2=plot(x,SDvasteMoy, 'b', 'LineWidth', 2);
hold on
plot(x,SDvasteMoy-SDvasteSd, 'b', 'LineWidth', 0.5);
hold on
plot(x,SDvasteMoy+SDvasteSd, 'b', 'LineWidth', 0.5);
xline(meanToeOffDroiteSyrine, 'b')
hold on
legend('Lina', 'Syrine');
       
pval=droitepval(:,8);
for k=1:101
    if pval(k)<0.05
        plot(k-1, 10,  'k*')
        hold on
    end
end
legend([p1 p2], {'Lina', 'Syrine'});

subplot(2,1,2)

p1=plot(x,LGvasteMoy, 'r', 'LineWidth', 2);
hold on
plot(x,LGvasteMoy-LGvasteSd, 'r', 'LineWidth', 0.5);
hold on
plot(x,LGvasteMoy+LGvasteSd, 'r', 'LineWidth', 0.5);
xline(meanToeOffGaucheLina, 'r')
title('Vaste Gauche ')
yline(100)
        
hold on
p2=plot(x,SGvasteMoy, 'b', 'LineWidth', 2);
hold on
plot(x,SGvasteMoy-SGvasteSd, 'b', 'LineWidth', 0.5);
hold on
plot(x,SGvasteMoy+SGvasteSd, 'b', 'LineWidth', 0.5);
xline(meanToeOffGaucheSyrine, 'b')
hold on
       
pval=gauchepval(:,8);
for k=1:101
    if pval(k)<0.05
        plot(k-1, 10,  'k*')
        hold on
    end
end        
        
legend([p1 p2], {'Lina', 'Syrine'});      
        
        
        
        

