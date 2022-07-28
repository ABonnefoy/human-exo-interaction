function indice= NotaNaNor0(vecteur)
    i=1;
    while( (~isnan(vecteur(i)) && vecteur(i)~=0.000) && i<length(vecteur)) 
        i=i+1;
    end
    if(isnan(vecteur(i)) || (vecteur(i)==0.000 && i~=1)) 
        fprintf('Valeur: %d', vecteur(i));
        indice = i-1;
    elseif(i>length(vecteur))
        disp('Pas de NaN dans ce vecteur');
        indice=length(vecteur);
    else
        %fprintf('Valeur: %d', vecteur(i));
        indice=length(vecteur);
    end
        
end