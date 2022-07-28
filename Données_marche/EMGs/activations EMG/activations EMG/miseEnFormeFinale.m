SyrineDroite=load('activationsDroiteSyrine'); SyrineDroite=SyrineDroite.activationsDroiteSyrine;
SyrineGauche=load('activationsGaucheSyrine'); SyrineGauche=SyrineGauche.activationsGaucheSyrine;
LinaDroite=load('activationsDroiteLina'); LinaDroite=LinaDroite.activationsDroiteLina;
LinaGauche=load('activationsGaucheLina'); LinaGauche=LinaGauche.activationsGaucheLina;

x=0:1:100;

% utiliser la fonction rectangle
% rajouter toe off et le nom du graphique hors de la fenetre graphique
% rajouter les déviations standards

% problème glutéal gauche / droite lina : à redéfinir
%% Biceps

% comparaison droite/gauche Lina Syrine

%% une activation 
figure;
title('biceps')
subplot(3,1,1)
axis([0 100 0 1]);

% lina droite
a=LinaDroite.biceps{1,1}(1);
b=LinaDroite.biceps{1,1}(2);

rLD1=rectangle('Position', [a 0.8 b-a 0.1]);
rLD1.FaceColor=[0.5 0.5 0.5];


% lina gauche 
a=LinaGauche.biceps{1,1}(1);
b=LinaGauche.biceps{1,1}(2);

rLG1=rectangle('Position', [a 0.6 b-a 0.1]);
rLG1.FaceColor=[0.5 0.5 0.5];

text(45, 0.85, 'Lina Droite');
text(45, 0.65, 'Lina Gauche');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% deux activations
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subplot(3,1,2)
axis([0 100 0 1]);

a=LinaDroite.biceps{2,1}(1,1);
b=LinaDroite.biceps{2,1}(1,2);

c=LinaDroite.biceps{2,1}(2,1);
d=LinaDroite.biceps{2,1}(2,2);

rLD1=rectangle('Position', [a 0.8 b-a 0.1]);
rLD1.FaceColor=[0.5 0.5 0.5];

rLD2=rectangle('Position', [c 0.8 d-c 0.1]);
rLD2.FaceColor=[0.5 0.5 0.5];

a=LinaGauche.biceps{2,1}(1,1);
b=LinaGauche.biceps{2,1}(1,2);

c=LinaGauche.biceps{2,1}(2,1);
d=LinaGauche.biceps{2,1}(2,2);

rLG1=rectangle('Position', [a 0.6 b-a 0.1]);
rLG1.FaceColor=[0.5 0.5 0.5];

rLG2=rectangle('Position', [c 0.6 d-c 0.1]);
rLG2.FaceColor=[0.5 0.5 0.5];
text(45, 0.85, 'Lina Droite');
text(45, 0.65, 'Lina Gauche');
% Syrine

a=SyrineDroite.biceps{2,1}(1,1);
b=SyrineDroite.biceps{2,1}(1,2);

c=SyrineDroite.biceps{2,1}(2,1);
d=SyrineDroite.biceps{2,1}(2,2);

rSD1=rectangle('Position', [a 0.4 b-a 0.1]);
rSD1.FaceColor=[0.5 0.5 0.5];

rSD2=rectangle('Position', [c 0.4 d-c 0.1]);
rSD2.FaceColor=[0.5 0.5 0.5];

a=SyrineGauche.biceps{2,1}(1,1);
b=SyrineGauche.biceps{2,1}(1,2);

c=SyrineGauche.biceps{2,1}(2,1);
d=SyrineGauche.biceps{2,1}(2,2);

rSG1=rectangle('Position', [a 0.2 b-a 0.1]);
rSG1.FaceColor=[0.5 0.5 0.5];

rSG2=rectangle('Position', [c 0.2 d-c 0.1]);
rSG2.FaceColor=[0.5 0.5 0.5];

text(45, 0.45, 'Syrine Droite');
text(45, 0.25, 'Syrine Gauche');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% trois activations
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subplot(3,1,3)
axis([0 100 0 1]);

a=LinaDroite.biceps{3,1}(1,1);
b=LinaDroite.biceps{3,1}(1,2);

c=LinaDroite.biceps{3,1}(2,1);
d=LinaDroite.biceps{3,1}(2,2);

e=LinaDroite.biceps{3,1}(3,1);
f=LinaDroite.biceps{3,1}(3,2);

rLD1=rectangle('Position', [a 0.8 b-a 0.1]);
rLD1.FaceColor=[0.5 0.5 0.5];

rLD2=rectangle('Position', [c 0.8 d-c 0.1]);
rLD2.FaceColor=[0.5 0.5 0.5];

rLD3=rectangle('Position', [e 0.8 f-e 0.1]);
rLD3.FaceColor=[0.5 0.5 0.5];

a=LinaGauche.biceps{3,1}(1,1);
b=LinaGauche.biceps{3,1}(1,2);

c=LinaGauche.biceps{3,1}(2,1);
d=LinaGauche.biceps{3,1}(2,2);

e=LinaGauche.biceps{3,1}(3,1);
f=LinaGauche.biceps{3,1}(3,2);

rLG1=rectangle('Position', [a 0.6 b-a 0.1]);
rLG1.FaceColor=[0.5 0.5 0.5];

rLG2=rectangle('Position', [c 0.6 d-c 0.1]);
rLG2.FaceColor=[0.5 0.5 0.5];

rLG3=rectangle('Position', [e 0.6 f-e 0.1]);
rLG3.FaceColor=[0.5 0.5 0.5];

text(45, 0.85, 'Lina Droite');
text(45, 0.65, 'Lina Gauche');

% Syrine

a=SyrineDroite.biceps{3,1}(1,1);
b=SyrineDroite.biceps{3,1}(1,2);

c=SyrineDroite.biceps{3,1}(2,1);
d=SyrineDroite.biceps{3,1}(2,2);

e=SyrineDroite.biceps{3,1}(3,1);
f=SyrineDroite.biceps{3,1}(3,2);

rSD1=rectangle('Position', [a 0.4 b-a 0.1]);
rSD1.FaceColor=[0.5 0.5 0.5];

rSD2=rectangle('Position', [c 0.4 d-c 0.1]);
rSD2.FaceColor=[0.5 0.5 0.5];

rSD3=rectangle('Position', [e 0.4 f-e 0.1]);
rSD3.FaceColor=[0.5 0.5 0.5];

a=SyrineGauche.biceps{3,1}(1,1);
b=SyrineGauche.biceps{3,1}(1,2);

c=SyrineGauche.biceps{3,1}(2,1);
d=SyrineGauche.biceps{3,1}(2,2);

e=SyrineGauche.biceps{3,1}(3,1);
f=SyrineGauche.biceps{3,1}(3,2);

rSG1=rectangle('Position', [a 0.2 b-a 0.1]);
rSG1.FaceColor=[0.5 0.5 0.5];

rSG2=rectangle('Position', [c 0.2 d-c 0.1]);
rSG2.FaceColor=[0.5 0.5 0.5];

rSG3=rectangle('Position', [e 0.2 f-e 0.1]);
rSG3.FaceColor=[0.5 0.5 0.5];

text(45, 0.45, 'Syrine Droite');
text(45, 0.25, 'Syrine Gauche');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Droit Antérieur
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 2 activations
figure;
title('Droit Antérieur')

subplot(2,1,1)
axis([0 100 0 1]);

a=LinaDroite.droitAnt{2,1}(1,1);
b=LinaDroite.droitAnt{2,1}(1,2);

c=LinaDroite.droitAnt{2,1}(2,1);
d=LinaDroite.droitAnt{2,1}(2,2);

rLD1=rectangle('Position', [a 0.8 b-a 0.1]);
rLD1.FaceColor=[0.5 0.5 0.5];

rLD2=rectangle('Position', [c 0.8 d-c 0.1]);
rLD2.FaceColor=[0.5 0.5 0.5];

a=LinaGauche.droitAnt{2,1}(1,1);
b=LinaGauche.droitAnt{2,1}(1,2);

c=LinaGauche.droitAnt{2,1}(2,1);
d=LinaGauche.droitAnt{2,1}(2,2);

rLG1=rectangle('Position', [a 0.6 b-a 0.1]);
rLG1.FaceColor=[0.5 0.5 0.5];

rLG2=rectangle('Position', [c 0.6 d-c 0.1]);
rLG2.FaceColor=[0.5 0.5 0.5];
text(45, 0.85, 'Lina Droite');
text(45, 0.65, 'Lina Gauche');

% Syrine

a=SyrineDroite.droitAnt{2,1}(1,1);
b=SyrineDroite.droitAnt{2,1}(1,2);

c=SyrineDroite.droitAnt{2,1}(2,1);
d=SyrineDroite.droitAnt{2,1}(2,2);

rSD1=rectangle('Position', [a 0.4 b-a 0.1]);
rSD1.FaceColor=[0.5 0.5 0.5];

rSD2=rectangle('Position', [c 0.4 d-c 0.1]);
rSD2.FaceColor=[0.5 0.5 0.5];

a=SyrineGauche.droitAnt{2,1}(1,1);
b=SyrineGauche.droitAnt{2,1}(1,2);

c=SyrineGauche.droitAnt{2,1}(2,1);
d=SyrineGauche.droitAnt{2,1}(2,2);

rSG1=rectangle('Position', [a 0.2 b-a 0.1]);
rSG1.FaceColor=[0.5 0.5 0.5];

rSG2=rectangle('Position', [c 0.2 d-c 0.1]);
rSG2.FaceColor=[0.5 0.5 0.5];

text(45, 0.45, 'Syrine Droite');
text(45, 0.25, 'Syrine Gauche');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 3 activations
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subplot(2,1,2)
axis([0 100 0 1]);

% pas de cycles suffisamment significatifs pour Lina

% Syrine

a=SyrineDroite.droitAnt{3,1}(1,1);
b=SyrineDroite.droitAnt{3,1}(1,2);

c=SyrineDroite.droitAnt{3,1}(2,1);
d=SyrineDroite.droitAnt{3,1}(2,2);

e=SyrineDroite.droitAnt{3,1}(3,1);
f=SyrineDroite.droitAnt{3,1}(3,2);

rSD1=rectangle('Position', [a 0.4 b-a 0.1]);
rSD1.FaceColor=[0.5 0.5 0.5];

rSD2=rectangle('Position', [c 0.4 d-c 0.1]);
rSD2.FaceColor=[0.5 0.5 0.5];

rSD3=rectangle('Position', [e 0.4 f-e 0.1]);
rSD3.FaceColor=[0.5 0.5 0.5];

a=SyrineGauche.droitAnt{3,1}(1,1);
b=SyrineGauche.droitAnt{3,1}(1,2);

c=SyrineGauche.droitAnt{3,1}(2,1);
d=SyrineGauche.droitAnt{3,1}(2,2);

e=SyrineGauche.droitAnt{3,1}(3,1);
f=SyrineGauche.droitAnt{3,1}(3,2);

rSG1=rectangle('Position', [a 0.2 b-a 0.1]);
rSG1.FaceColor=[0.5 0.5 0.5];

rSG2=rectangle('Position', [c 0.2 d-c 0.1]);
rSG2.FaceColor=[0.5 0.5 0.5];

rSG3=rectangle('Position', [e 0.2 f-e 0.1]);
rSG3.FaceColor=[0.5 0.5 0.5];

text(45, 0.45, 'Syrine Droite');
text(45, 0.25, 'Syrine Gauche');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Gastrocnémien
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% 1 et 2 activations

% 1 activation
figure;
title('Gastrocnémien')
subplot(2,1,1)
axis([0 100 0 1]);

% lina droite
a=LinaDroite.gastrocnemien{1,1}(1);
b=LinaDroite.gastrocnemien{1,1}(2);

rLD1=rectangle('Position', [a 0.8 b-a 0.1]);
rLD1.FaceColor=[0.5 0.5 0.5];


% lina gauche 
a=LinaGauche.gastrocnemien{1,1}(1);
b=LinaGauche.gastrocnemien{1,1}(2);

rLG1=rectangle('Position', [a 0.6 b-a 0.1]);
rLG1.FaceColor=[0.5 0.5 0.5];

text(45, 0.85, 'Lina Droite');
text(45, 0.65, 'Lina Gauche');


a=SyrineDroite.gastrocnemien{1,1}(1);
b=SyrineDroite.gastrocnemien{1,1}(2);

rSD1=rectangle('Position', [a 0.4 b-a 0.1]);
rSD1.FaceColor=[0.5 0.5 0.5];


% Syrine gauche 
a=SyrineGauche.gastrocnemien{1,1}(1);
b=SyrineGauche.gastrocnemien{1,1}(2);

rSG1=rectangle('Position', [a 0.2 b-a 0.1]);
rSG1.FaceColor=[0.5 0.5 0.5];

text(45, 0.45, 'Syrine Droite');
text(45, 0.25, 'Syrine Gauche');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 2 activations

subplot(2,1,2)
axis([0 100 0 1]);

a=LinaDroite.gastrocnemien{2,1}(1,1);
b=LinaDroite.gastrocnemien{2,1}(1,2);

c=LinaDroite.gastrocnemien{2,1}(2,1);
d=LinaDroite.gastrocnemien{2,1}(2,2);

rLD1=rectangle('Position', [a 0.8 b-a 0.1]);
rLD1.FaceColor=[0.5 0.5 0.5];

rLD2=rectangle('Position', [c 0.8 d-c 0.1]);
rLD2.FaceColor=[0.5 0.5 0.5];

a=LinaGauche.gastrocnemien{2,1}(1,1);
b=LinaGauche.gastrocnemien{2,1}(1,2);

c=LinaGauche.gastrocnemien{2,1}(2,1);
d=LinaGauche.gastrocnemien{2,1}(2,2);

rLG1=rectangle('Position', [a 0.6 b-a 0.1]);
rLG1.FaceColor=[0.5 0.5 0.5];

rLG2=rectangle('Position', [c 0.6 d-c 0.1]);
rLG2.FaceColor=[0.5 0.5 0.5];
text(45, 0.85, 'Lina Droite');
text(45, 0.65, 'Lina Gauche');

% Syrine

a=SyrineDroite.gastrocnemien{2,1}(1,1);
b=SyrineDroite.gastrocnemien{2,1}(1,2);

c=SyrineDroite.gastrocnemien{2,1}(2,1);
d=SyrineDroite.gastrocnemien{2,1}(2,2);

rSD1=rectangle('Position', [a 0.4 b-a 0.1]);
rSD1.FaceColor=[0.5 0.5 0.5];

rSD2=rectangle('Position', [c 0.4 d-c 0.1]);
rSD2.FaceColor=[0.5 0.5 0.5];

a=SyrineGauche.gastrocnemien{2,1}(1,1);
b=SyrineGauche.gastrocnemien{2,1}(1,2);

c=SyrineGauche.gastrocnemien{2,1}(2,1);
d=SyrineGauche.gastrocnemien{2,1}(2,2);

rSG1=rectangle('Position', [a 0.2 b-a 0.1]);
rSG1.FaceColor=[0.5 0.5 0.5];

rSG2=rectangle('Position', [c 0.2 d-c 0.1]);
rSG2.FaceColor=[0.5 0.5 0.5];

text(45, 0.45, 'Syrine Droite');
text(45, 0.25, 'Syrine Gauche');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Glutéal
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% 1,2,3 activations

figure;
title('Glutéal')
subplot(3,1,1)
axis([0 100 0 1]);

% lina droite
a=LinaDroite.gluteal{1,1}(1);
b=LinaDroite.gluteal{1,1}(2);

rLD1=rectangle('Position', [a 0.8 b-a 0.1]);
rLD1.FaceColor=[0.5 0.5 0.5];


% lina gauche 
a=LinaGauche.gluteal{1,1}(1);
b=LinaGauche.gluteal{1,1}(2);

rLG1=rectangle('Position', [a 0.6 b-a 0.1]);
rLG1.FaceColor=[0.5 0.5 0.5];

text(45, 0.85, 'Lina Droite');
text(45, 0.65, 'Lina Gauche');


a=SyrineDroite.gluteal{1,1}(1);
b=SyrineDroite.gluteal{1,1}(2);

rSD1=rectangle('Position', [a 0.4 b-a 0.1]);
rSD1.FaceColor=[0.5 0.5 0.5];


% Syrine gauche 
a=SyrineGauche.gluteal{1,1}(1);
b=SyrineGauche.gluteal{1,1}(2);

rSG1=rectangle('Position', [a 0.2 b-a 0.1]);
rSG1.FaceColor=[0.5 0.5 0.5];

text(45, 0.45, 'Syrine Droite');
text(45, 0.25, 'Syrine Gauche');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 2 activations

subplot(3,1,2)
axis([0 100 0 1]);

a=LinaDroite.gluteal{2,1}(1,1);
b=LinaDroite.gluteal{2,1}(1,2);

c=LinaDroite.gluteal{2,1}(2,1);
d=LinaDroite.gluteal{2,1}(2,2);

rLD1=rectangle('Position', [a 0.8 b-a 0.1]);
rLD1.FaceColor=[0.5 0.5 0.5];

rLD2=rectangle('Position', [c 0.8 d-c 0.1]);
rLD2.FaceColor=[0.5 0.5 0.5];

a=LinaGauche.gluteal{2,1}(1,1);
b=LinaGauche.gluteal{2,1}(1,2);

c=LinaGauche.gluteal{2,1}(2,1);
d=LinaGauche.gluteal{2,1}(2,2);

rLG1=rectangle('Position', [a 0.6 b-a 0.1]);
rLG1.FaceColor=[0.5 0.5 0.5];

rLG2=rectangle('Position', [c 0.6 d-c 0.1]);
rLG2.FaceColor=[0.5 0.5 0.5];
text(45, 0.85, 'Lina Droite');
text(45, 0.65, 'Lina Gauche');

% Syrine

a=SyrineDroite.gluteal{2,1}(1,1);
b=SyrineDroite.gluteal{2,1}(1,2);

c=SyrineDroite.gluteal{2,1}(2,1);
d=SyrineDroite.gluteal{2,1}(2,2);

rSD1=rectangle('Position', [a 0.4 b-a 0.1]);
rSD1.FaceColor=[0.5 0.5 0.5];

rSD2=rectangle('Position', [c 0.4 d-c 0.1]);
rSD2.FaceColor=[0.5 0.5 0.5];

a=SyrineGauche.gluteal{2,1}(1,1);
b=SyrineGauche.gluteal{2,1}(1,2);

c=SyrineGauche.gluteal{2,1}(2,1);
d=SyrineGauche.gluteal{2,1}(2,2);

rSG1=rectangle('Position', [a 0.2 b-a 0.1]);
rSG1.FaceColor=[0.5 0.5 0.5];

rSG2=rectangle('Position', [c 0.2 d-c 0.1]);
rSG2.FaceColor=[0.5 0.5 0.5];

text(45, 0.45, 'Syrine Droite');
text(45, 0.25, 'Syrine Gauche');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 3 activations

subplot(3,1,3)
axis([0 100 0 1]);

a=LinaDroite.gluteal{3,1}(1,1);
b=LinaDroite.gluteal{3,1}(1,2);

c=LinaDroite.gluteal{3,1}(2,1);
d=LinaDroite.gluteal{3,1}(2,2);

e=LinaDroite.gluteal{3,1}(3,1);
f=LinaDroite.gluteal{3,1}(3,2);

rLD1=rectangle('Position', [a 0.8 b-a 0.1]);
rLD1.FaceColor=[0.5 0.5 0.5];

rLD2=rectangle('Position', [c 0.8 d-c 0.1]);
rLD2.FaceColor=[0.5 0.5 0.5];

rLD3=rectangle('Position', [e 0.8 f-e 0.1]);
rLD3.FaceColor=[0.5 0.5 0.5];

a=LinaGauche.gluteal{3,1}(1,1);
b=LinaGauche.gluteal{3,1}(1,2);

c=LinaGauche.gluteal{3,1}(2,1);
d=LinaGauche.gluteal{3,1}(2,2);

e=LinaGauche.gluteal{3,1}(3,1);
f=LinaGauche.gluteal{3,1}(3,2);

rLG1=rectangle('Position', [a 0.6 b-a 0.1]);
rLG1.FaceColor=[0.5 0.5 0.5];

rLG2=rectangle('Position', [c 0.6 d-c 0.1]);
rLG2.FaceColor=[0.5 0.5 0.5];

% rLG3=rectangle('Position', [e 0.6 f-e 0.1]);
% rLG3.FaceColor=[0.5 0.5 0.5];

text(45, 0.85, 'Lina Droite');
text(45, 0.65, 'Lina Gauche');

% Syrine

a=SyrineDroite.gluteal{3,1}(1,1);
b=SyrineDroite.gluteal{3,1}(1,2);

c=SyrineDroite.gluteal{3,1}(2,1);
d=SyrineDroite.gluteal{3,1}(2,2);

e=SyrineDroite.gluteal{3,1}(3,1);
f=SyrineDroite.gluteal{3,1}(3,2);

rSD1=rectangle('Position', [a 0.4 b-a 0.1]);
rSD1.FaceColor=[0.5 0.5 0.5];

rSD2=rectangle('Position', [c 0.4 d-c 0.1]);
rSD2.FaceColor=[0.5 0.5 0.5];

rSD3=rectangle('Position', [e 0.4 f-e 0.1]);
rSD3.FaceColor=[0.5 0.5 0.5];

a=SyrineGauche.gluteal{3,1}(1,1);
b=SyrineGauche.gluteal{3,1}(1,2);

c=SyrineGauche.gluteal{3,1}(2,1);
d=SyrineGauche.gluteal{3,1}(2,2);

e=SyrineGauche.gluteal{3,1}(3,1);
f=SyrineGauche.gluteal{3,1}(3,2);

rSG1=rectangle('Position', [a 0.2 b-a 0.1]);
rSG1.FaceColor=[0.5 0.5 0.5];

rSG2=rectangle('Position', [c 0.2 d-c 0.1]);
rSG2.FaceColor=[0.5 0.5 0.5];

rSG3=rectangle('Position', [e 0.2 f-e 0.1]);
rSG3.FaceColor=[0.5 0.5 0.5];

text(45, 0.45, 'Syrine Droite');
text(45, 0.25, 'Syrine Gauche');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Soléaire
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% 1 et 2 activations

figure;
title('soléaire')
subplot(2,1,1)
axis([0 100 0 1]);

% lina droite
a=LinaDroite.soleaire{1,1}(1);
b=LinaDroite.soleaire{1,1}(2);

rLD1=rectangle('Position', [a 0.8 b-a 0.1]);
rLD1.FaceColor=[0.5 0.5 0.5];


% lina gauche 
a=LinaGauche.soleaire{1,1}(1);
b=LinaGauche.soleaire{1,1}(2);

rLG1=rectangle('Position', [a 0.6 b-a 0.1]);
rLG1.FaceColor=[0.5 0.5 0.5];

text(45, 0.85, 'Lina Droite');
text(45, 0.65, 'Lina Gauche');


a=SyrineDroite.soleaire{1,1}(1);
b=SyrineDroite.soleaire{1,1}(2);

rSD1=rectangle('Position', [a 0.4 b-a 0.1]);
rSD1.FaceColor=[0.5 0.5 0.5];


% Syrine gauche 
a=SyrineGauche.soleaire{1,1}(1);
b=SyrineGauche.soleaire{1,1}(2);

rSG1=rectangle('Position', [a 0.2 b-a 0.1]);
rSG1.FaceColor=[0.5 0.5 0.5];

text(45, 0.45, 'Syrine Droite');
text(45, 0.25, 'Syrine Gauche');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 2 activations

subplot(2,1,2)
axis([0 100 0 1]);

a=LinaDroite.soleaire{2,1}(1,1);
b=LinaDroite.soleaire{2,1}(1,2);

c=LinaDroite.soleaire{2,1}(2,1);
d=LinaDroite.soleaire{2,1}(2,2);

rLD1=rectangle('Position', [a 0.8 b-a 0.1]);
rLD1.FaceColor=[0.5 0.5 0.5];

rLD2=rectangle('Position', [c 0.8 d-c 0.1]);
rLD2.FaceColor=[0.5 0.5 0.5];

a=LinaGauche.soleaire{2,1}(1,1);
b=LinaGauche.soleaire{2,1}(1,2);

c=LinaGauche.soleaire{2,1}(2,1);
d=LinaGauche.soleaire{2,1}(2,2);

rLG1=rectangle('Position', [a 0.6 b-a 0.1]);
rLG1.FaceColor=[0.5 0.5 0.5];

rLG2=rectangle('Position', [c 0.6 d-c 0.1]);
rLG2.FaceColor=[0.5 0.5 0.5];
text(45, 0.85, 'Lina Droite');
text(45, 0.65, 'Lina Gauche');

% Syrine

a=SyrineDroite.soleaire{2,1}(1,1);
b=SyrineDroite.soleaire{2,1}(1,2);

c=SyrineDroite.soleaire{2,1}(2,1);
d=SyrineDroite.soleaire{2,1}(2,2);

rSD1=rectangle('Position', [a 0.4 b-a 0.1]);
rSD1.FaceColor=[0.5 0.5 0.5];

rSD2=rectangle('Position', [c 0.4 d-c 0.1]);
rSD2.FaceColor=[0.5 0.5 0.5];

a=SyrineGauche.soleaire{2,1}(1,1);
b=SyrineGauche.soleaire{2,1}(1,2);

c=SyrineGauche.soleaire{2,1}(2,1);
d=SyrineGauche.soleaire{2,1}(2,2);

rSG1=rectangle('Position', [a 0.2 b-a 0.1]);
rSG1.FaceColor=[0.5 0.5 0.5];

rSG2=rectangle('Position', [c 0.2 d-c 0.1]);
rSG2.FaceColor=[0.5 0.5 0.5];

text(45, 0.45, 'Syrine Droite');
text(45, 0.25, 'Syrine Gauche');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Tibial Antérieur
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% 1 activation pour Lina, 2,3,4 activations pour Syrine


figure;
title('tibialis anterior')
subplot(4,1,1)
axis([0 100 0 1]);

% lina droite
a=LinaDroite.tibialAnt{1,1}(1);
b=LinaDroite.tibialAnt{1,1}(2);

rLD1=rectangle('Position', [a 0.8 b-a 0.1]);
rLD1.FaceColor=[0.5 0.5 0.5];


% lina gauche 
a=LinaGauche.tibialAnt{1,1}(1);
b=LinaGauche.tibialAnt{1,1}(2);

rLG1=rectangle('Position', [a 0.6 b-a 0.1]);
rLG1.FaceColor=[0.5 0.5 0.5];

text(45, 0.85, 'Lina Droite');
text(45, 0.65, 'Lina Gauche');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% 2 activations pour Syrine

subplot(4,1,2)
axis([0 100 0 1]);
a=SyrineDroite.tibialAnt{2,1}(1,1);
b=SyrineDroite.tibialAnt{2,1}(1,2);

c=SyrineDroite.tibialAnt{2,1}(2,1);
d=SyrineDroite.tibialAnt{2,1}(2,2);

rSD1=rectangle('Position', [a 0.4 b-a 0.1]);
rSD1.FaceColor=[0.5 0.5 0.5];

rSD2=rectangle('Position', [c 0.4 d-c 0.1]);
rSD2.FaceColor=[0.5 0.5 0.5];

a=SyrineGauche.tibialAnt{2,1}(1,1);
b=SyrineGauche.tibialAnt{2,1}(1,2);

c=SyrineGauche.tibialAnt{2,1}(2,1);
d=SyrineGauche.tibialAnt{2,1}(2,2);

rSG1=rectangle('Position', [a 0.2 b-a 0.1]);
rSG1.FaceColor=[0.5 0.5 0.5];

rSG2=rectangle('Position', [c 0.2 d-c 0.1]);
rSG2.FaceColor=[0.5 0.5 0.5];

text(45, 0.45, 'Syrine Droite');
text(45, 0.25, 'Syrine Gauche');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 3 activations

subplot(4,1,3)
axis([0 100 0 1]);
a=SyrineDroite.tibialAnt{3,1}(1,1);
b=SyrineDroite.tibialAnt{3,1}(1,2);

c=SyrineDroite.tibialAnt{3,1}(2,1);
d=SyrineDroite.tibialAnt{3,1}(2,2);

e=SyrineDroite.tibialAnt{3,1}(3,1);
f=SyrineDroite.tibialAnt{3,1}(3,2);

rSD1=rectangle('Position', [a 0.4 b-a 0.1]);
rSD1.FaceColor=[0.5 0.5 0.5];

rSD2=rectangle('Position', [c 0.4 d-c 0.1]);
rSD2.FaceColor=[0.5 0.5 0.5];

rSD3=rectangle('Position', [e 0.4 f-e 0.1]);
rSD3.FaceColor=[0.5 0.5 0.5];

a=SyrineGauche.tibialAnt{3,1}(1,1);
b=SyrineGauche.tibialAnt{3,1}(1,2);

c=SyrineGauche.tibialAnt{3,1}(2,1);
d=SyrineGauche.tibialAnt{3,1}(2,2);

e=SyrineGauche.tibialAnt{3,1}(3,1);
f=SyrineGauche.tibialAnt{3,1}(3,2);

rSG1=rectangle('Position', [a 0.2 b-a 0.1]);
rSG1.FaceColor=[0.5 0.5 0.5];

rSG2=rectangle('Position', [c 0.2 d-c 0.1]);
rSG2.FaceColor=[0.5 0.5 0.5];

rSG3=rectangle('Position', [e 0.2 f-e 0.1]);
rSG3.FaceColor=[0.5 0.5 0.5];

text(45, 0.45, 'Syrine Droite');
text(45, 0.25, 'Syrine Gauche');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 4 activations

subplot(4,1,4)
axis([0 100 0 1]);
a=SyrineDroite.tibialAnt{4,1}(1,1);
b=SyrineDroite.tibialAnt{4,1}(1,2);

c=SyrineDroite.tibialAnt{4,1}(2,1);
d=SyrineDroite.tibialAnt{4,1}(2,2);

e=SyrineDroite.tibialAnt{4,1}(3,1);
f=SyrineDroite.tibialAnt{4,1}(3,2);

g=SyrineDroite.tibialAnt{4,1}(3,1);
h=SyrineDroite.tibialAnt{4,1}(3,2);

rSD1=rectangle('Position', [a 0.4 b-a 0.1]);
rSD1.FaceColor=[0.5 0.5 0.5];

rSD2=rectangle('Position', [c 0.4 d-c 0.1]);
rSD2.FaceColor=[0.5 0.5 0.5];

rSD3=rectangle('Position', [e 0.4 f-e 0.1]);
rSD3.FaceColor=[0.5 0.5 0.5];

rSD4=rectangle('Position', [g 0.4 h-g 0.1]);
rSD4.FaceColor=[0.5 0.5 0.5];

a=SyrineGauche.tibialAnt{4,1}(1,1);
b=SyrineGauche.tibialAnt{4,1}(1,2);

c=SyrineGauche.tibialAnt{4,1}(2,1);
d=SyrineGauche.tibialAnt{4,1}(2,2);

e=SyrineGauche.tibialAnt{4,1}(3,1);
f=SyrineGauche.tibialAnt{4,1}(3,2);

g=SyrineGauche.tibialAnt{4,1}(4,1);
h=SyrineGauche.tibialAnt{4,1}(4,2);

rSG1=rectangle('Position', [a 0.2 b-a 0.1]);
rSG1.FaceColor=[0.5 0.5 0.5];

rSG2=rectangle('Position', [c 0.2 d-c 0.1]);
rSG2.FaceColor=[0.5 0.5 0.5];

rSG3=rectangle('Position', [e 0.2 f-e 0.1]);
rSG3.FaceColor=[0.5 0.5 0.5];

rSG4=rectangle('Position', [g 0.2 h-g 0.1]);
rSG4.FaceColor=[0.5 0.5 0.5];

text(45, 0.45, 'Syrine Droite');
text(45, 0.25, 'Syrine Gauche');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Tendineux
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% 1,2,3 activations pour Lina et 2,3 activations pour Syrine

figure;
title('tendineux')
subplot(3,1,1)
axis([0 100 0 1]);

% lina droite
a=LinaDroite.tendineux{1,1}(1);
b=LinaDroite.tendineux{1,1}(2);

rLD1=rectangle('Position', [a 0.8 b-a 0.1]);
rLD1.FaceColor=[0.5 0.5 0.5];

% lina gauche 
a=LinaGauche.tendineux{1,1}(1);
b=LinaGauche.tendineux{1,1}(2);

rLG1=rectangle('Position', [a 0.6 b-a 0.1]);
rLG1.FaceColor=[0.5 0.5 0.5];

text(45, 0.85, 'Lina Droite');
text(45, 0.65, 'Lina Gauche');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 2 activations

subplot(3,1,2)
axis([0 100 0 1]);

a=LinaDroite.tendineux{2,1}(1,1);
b=LinaDroite.tendineux{2,1}(1,2);

c=LinaDroite.tendineux{2,1}(2,1);
d=LinaDroite.tendineux{2,1}(2,2);

rLD1=rectangle('Position', [a 0.8 b-a 0.1]);
rLD1.FaceColor=[0.5 0.5 0.5];

rLD2=rectangle('Position', [c 0.8 d-c 0.1]);
rLD2.FaceColor=[0.5 0.5 0.5];

a=LinaGauche.tendineux{2,1}(1,1);
b=LinaGauche.tendineux{2,1}(1,2);

c=LinaGauche.tendineux{2,1}(2,1);
d=LinaGauche.tendineux{2,1}(2,2);

rLG1=rectangle('Position', [a 0.6 b-a 0.1]);
rLG1.FaceColor=[0.5 0.5 0.5];

rLG2=rectangle('Position', [c 0.6 d-c 0.1]);
rLG2.FaceColor=[0.5 0.5 0.5];
text(45, 0.85, 'Lina Droite');
text(45, 0.65, 'Lina Gauche');
% Syrine

a=SyrineDroite.tendineux{2,1}(1,1);
b=SyrineDroite.tendineux{2,1}(1,2);

c=SyrineDroite.tendineux{2,1}(2,1);
d=SyrineDroite.tendineux{2,1}(2,2);

rSD1=rectangle('Position', [a 0.4 b-a 0.1]);
rSD1.FaceColor=[0.5 0.5 0.5];

rSD2=rectangle('Position', [c 0.4 d-c 0.1]);
rSD2.FaceColor=[0.5 0.5 0.5];

a=SyrineGauche.tendineux{2,1}(1,1);
b=SyrineGauche.tendineux{2,1}(1,2);

c=SyrineGauche.tendineux{2,1}(2,1);
d=SyrineGauche.tendineux{2,1}(2,2);

rSG1=rectangle('Position', [a 0.2 b-a 0.1]);
rSG1.FaceColor=[0.5 0.5 0.5];

rSG2=rectangle('Position', [c 0.2 d-c 0.1]);
rSG2.FaceColor=[0.5 0.5 0.5];

text(45, 0.45, 'Syrine Droite');
text(45, 0.25, 'Syrine Gauche');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 3 activations

subplot(3,1,3)
axis([0 100 0 1]);

a=LinaDroite.tendineux{3,1}(1,1);
b=LinaDroite.tendineux{3,1}(1,2);

c=LinaDroite.tendineux{3,1}(2,1);
d=LinaDroite.tendineux{3,1}(2,2);

e=LinaDroite.tendineux{3,1}(3,1);
f=LinaDroite.tendineux{3,1}(3,2);

rLD1=rectangle('Position', [a 0.8 b-a 0.1]);
rLD1.FaceColor=[0.5 0.5 0.5];

rLD2=rectangle('Position', [c 0.8 d-c 0.1]);
rLD2.FaceColor=[0.5 0.5 0.5];

rLD3=rectangle('Position', [e 0.8 f-e 0.1]);
rLD3.FaceColor=[0.5 0.5 0.5];

a=LinaGauche.tendineux{3,1}(1,1);
b=LinaGauche.tendineux{3,1}(1,2);

c=LinaGauche.tendineux{3,1}(2,1);
d=LinaGauche.tendineux{3,1}(2,2);

e=LinaGauche.tendineux{3,1}(3,1);
f=LinaGauche.tendineux{3,1}(3,2);

rLG1=rectangle('Position', [a 0.6 b-a 0.1]);
rLG1.FaceColor=[0.5 0.5 0.5];

rLG2=rectangle('Position', [c 0.6 d-c 0.1]);
rLG2.FaceColor=[0.5 0.5 0.5];

rLG3=rectangle('Position', [e 0.6 f-e 0.1]);
rLG3.FaceColor=[0.5 0.5 0.5];

text(45, 0.85, 'Lina Droite');
text(45, 0.65, 'Lina Gauche');

% Syrine

a=SyrineDroite.tendineux{3,1}(1,1);
b=SyrineDroite.tendineux{3,1}(1,2);

c=SyrineDroite.tendineux{3,1}(2,1);
d=SyrineDroite.tendineux{3,1}(2,2);

e=SyrineDroite.tendineux{3,1}(3,1);
f=SyrineDroite.tendineux{3,1}(3,2);

rSD1=rectangle('Position', [a 0.4 b-a 0.1]);
rSD1.FaceColor=[0.5 0.5 0.5];

rSD2=rectangle('Position', [c 0.4 d-c 0.1]);
rSD2.FaceColor=[0.5 0.5 0.5];

rSD3=rectangle('Position', [e 0.4 f-e 0.1]);
rSD3.FaceColor=[0.5 0.5 0.5];

a=SyrineGauche.tendineux{3,1}(1,1);
b=SyrineGauche.tendineux{3,1}(1,2);

c=SyrineGauche.tendineux{3,1}(2,1);
d=SyrineGauche.tendineux{3,1}(2,2);

e=SyrineGauche.tendineux{3,1}(3,1);
f=SyrineGauche.tendineux{3,1}(3,2);

rSG1=rectangle('Position', [a 0.2 b-a 0.1]);
rSG1.FaceColor=[0.5 0.5 0.5];

rSG2=rectangle('Position', [c 0.2 d-c 0.1]);
rSG2.FaceColor=[0.5 0.5 0.5];

rSG3=rectangle('Position', [e 0.2 f-e 0.1]);
rSG3.FaceColor=[0.5 0.5 0.5];

text(45, 0.45, 'Syrine Droite');
text(45, 0.25, 'Syrine Gauche');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Vaste

% 2 activations pour Lina, 2 & 3 pour Syrine

figure;
title('Vastus medialis')
subplot(2,1,1)
axis([0 100 0 1]);

a=LinaDroite.vaste{2,1}(1,1);
b=LinaDroite.vaste{2,1}(1,2);

c=LinaDroite.vaste{2,1}(2,1);
d=LinaDroite.vaste{2,1}(2,2);

rLD1=rectangle('Position', [a 0.8 b-a 0.1]);
rLD1.FaceColor=[0.5 0.5 0.5];

rLD2=rectangle('Position', [c 0.8 d-c 0.1]);
rLD2.FaceColor=[0.5 0.5 0.5];

a=LinaGauche.vaste{2,1}(1,1);
b=LinaGauche.vaste{2,1}(1,2);

c=LinaGauche.vaste{2,1}(2,1);
d=LinaGauche.vaste{2,1}(2,2);

rLG1=rectangle('Position', [a 0.6 b-a 0.1]);
rLG1.FaceColor=[0.5 0.5 0.5];

rLG2=rectangle('Position', [c 0.6 d-c 0.1]);
rLG2.FaceColor=[0.5 0.5 0.5];
text(45, 0.85, 'Lina Droite');
text(45, 0.65, 'Lina Gauche');
% Syrine

a=SyrineDroite.vaste{2,1}(1,1);
b=SyrineDroite.vaste{2,1}(1,2);

c=SyrineDroite.vaste{2,1}(2,1);
d=SyrineDroite.vaste{2,1}(2,2);

rSD1=rectangle('Position', [a 0.4 b-a 0.1]);
rSD1.FaceColor=[0.5 0.5 0.5];

rSD2=rectangle('Position', [c 0.4 d-c 0.1]);
rSD2.FaceColor=[0.5 0.5 0.5];

a=SyrineGauche.vaste{2,1}(1,1);
b=SyrineGauche.vaste{2,1}(1,2);

c=SyrineGauche.vaste{2,1}(2,1);
d=SyrineGauche.vaste{2,1}(2,2);

rSG1=rectangle('Position', [a 0.2 b-a 0.1]);
rSG1.FaceColor=[0.5 0.5 0.5];

rSG2=rectangle('Position', [c 0.2 d-c 0.1]);
rSG2.FaceColor=[0.5 0.5 0.5];

text(45, 0.45, 'Syrine Droite');
text(45, 0.25, 'Syrine Gauche');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 3 activations
subplot(2,1,2)
axis([0 100 0 1]);
a=SyrineDroite.vaste{3,1}(1,1);
b=SyrineDroite.vaste{3,1}(1,2);

c=SyrineDroite.vaste{3,1}(2,1);
d=SyrineDroite.vaste{3,1}(2,2);

e=SyrineDroite.vaste{3,1}(3,1);
f=SyrineDroite.vaste{3,1}(3,2);

rSD1=rectangle('Position', [a 0.4 b-a 0.1]);
rSD1.FaceColor=[0.5 0.5 0.5];

rSD2=rectangle('Position', [c 0.4 d-c 0.1]);
rSD2.FaceColor=[0.5 0.5 0.5];

rSD3=rectangle('Position', [e 0.4 f-e 0.1]);
rSD3.FaceColor=[0.5 0.5 0.5];

a=SyrineGauche.vaste{3,1}(1,1);
b=SyrineGauche.vaste{3,1}(1,2);

c=SyrineGauche.vaste{3,1}(2,1);
d=SyrineGauche.vaste{3,1}(2,2);

e=SyrineGauche.vaste{3,1}(3,1);
f=SyrineGauche.vaste{3,1}(3,2);

rSG1=rectangle('Position', [a 0.2 b-a 0.1]);
rSG1.FaceColor=[0.5 0.5 0.5];

rSG2=rectangle('Position', [c 0.2 d-c 0.1]);
rSG2.FaceColor=[0.5 0.5 0.5];

rSG3=rectangle('Position', [e 0.2 f-e 0.1]);
rSG3.FaceColor=[0.5 0.5 0.5];

text(45, 0.45, 'Syrine Droite');
text(45, 0.25, 'Syrine Gauche');



