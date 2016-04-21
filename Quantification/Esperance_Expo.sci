function [res] = Esperance_Expo(intensite,sup,inf)
    intensite = abs(intensite);
    inf = max(inf,0);
    sup = max(sup,0);
    res = inf.*exp(-1*intensite*inf) + (exp(-1*intensite*inf))/(intensite) - sup.*exp(-1*intensite*sup) + (exp(-1*intensite*sup))/(intensite);
endfunction

function [res] = PHI_EXPO(intensite,x)
    //[P,Q]=cdfgam("PQ",x,Shape,Rate)
    p = length(x);
    UN = ones(1,p);
    res = cdfgam("PQ",x,UN,(1/intensite)*UN);
endfunction

function [res] = Esperance_Gauss(sup,inf)
res = (exp(-0.5.*inf.^2) -  exp(-0.5.*sup.^2)) /(sqrt(%pi * 2));
endfunction

function [res] = PHI(x)// fonction de repartition
    res = cdfnor("PQ",x,zeros(x),ones(x));
endfunction

function [res] = CallBS(grill,poid,S0,K,r,v,t)// fonction de repartition
tmp = S0*exp((r-v*v/2)*t + (v*sqrt(t))*grill) - K;
tmp = max(tmp,zeros(grill));
//disp(tmp)
res = tmp*poid';
res = exp(-r*t)*res; 
endfunction

function [res] = PutBS(grill,poid,S0,K,r,v,t)
tmp = K - S0*exp((r-v*v/2)*t + (v*sqrt(t))*grill);
tmp = max(tmp,zeros(grill));
//disp(tmp)
res = tmp*poid';
res = exp(-r*t)*res;  
endfunction


function [Grilles,Poids] = Lloyd_1d(nb_quant,init,nb_iter)
    compteur = 0;
    Grilles = init;
    Poids = ones(1,nb_quant);
  
    for compteur =1:nb_iter  // || les centres ne bougent plus*/)
    //disp('iteration numero '+string(compteur))
    tmp_memoire = ( cat(2,Grilles,%inf) + cat(2,-%inf,Grilles) );
    tmp_memoire = 0.5*tmp_memoire;// -inf;demi somme;+inf
    tmp_memoirem = tmp_memoire(1:nb_quant);//x+1/2
    tmp_memoirep = tmp_memoire(2:nb_quant+1);//x-1/2
    espe = Esperance_Gauss(tmp_memoirep,tmp_memoirem);

    Poids = PHI(tmp_memoirep) - PHI(tmp_memoirem);
    Grilles = espe ./Poids;
   end

    Grilles = Grilles';
    Poids = Poids';
    plot(Grilles,Poids)
endfunction



 function [res] = Distorsion(Grille,xpp,xmm) 
     K = length(xmm);
 A = xmm.* exp(-xmm.*xmm/2) - xpp.* exp(-xpp.*xpp/2) + PHI(xpp) - PHI (xmm);
 A(1) = - xpp(1)*exp(-xpp(1)*xpp(1)/2) + PHI(xpp(1));
 A(K) = xmm(K)*exp(-xmm(K)*xmm(K)/2) + 1 - PHI(xmm(K));
 //disp(A, "A")
 B = 2*Grille.*(exp(-xpp.*xpp/2) - exp(-xmm.*xmm/2));
// disp(B , "B")
 C = Grille.*Grille.*(PHI(xpp)-PHI(xmm))
 // disp(C , "C")
 tmp = (1/sqrt(2*%pi)) * (A + B + C);
 res = sum(tmp);
endfunction
 
 
 
function [Grilles,Poids] = Lloyd_expo_1d(nb_quant,init,nb_iter,intensite)
    compteur = 0;
    Grilles = zeros(1,nb_quant);
    Poids = zeros(1,nb_quant);
   disp("debut dla boucle")
    while compteur <= nb_iter  // || les centres ne bougent plus*/)
    //disp(compteur ,'iteration numero ')
    tmp_memoire = 0.5*( cat(2,init,%inf) + cat(2,-%inf,init) );
    tmp_memoirem = tmp_memoire(1:nb_quant);
    tmp_memoirep = tmp_memoire(2:nb_quant+1);
    disp(" tmp memore ")
     espe = Esperance_Expo(intensite,tmp_memoirep,tmp_memoirem);
     disp(" espe ")
     Poids = diff(PHI_EXPO(intensite,tmp_memoire));
     disp(" phi ")
 //    disp(length(Poids))
   //  disp(length(espe))     
     Grilles = espe ./Poids;
     tmp_memoire = Grilles; 
     compteur = compteur + 1;
   end
endfunction


 function [Grilles,Poids] = Lloyd2BS(nb_quant,nb_step,nb_iter,Mu,Sigma,x0,init,T)
        dt = T/nb_step;
        tmp = ones(1,nb_quant);
        Grilles = x0*ones(nb_step,nb_quant);
        Poids = ones(nb_step,nb_quant);
        for t=2:nb_step
        tmp_memoire = 0.5*( cat(2,Grilles(t-1,:),%inf) + cat(2,-%inf,Grilles(t-1,:)) );
        tmp_plus  = tmp_memoire(1:nb_quant);
        tmp_moins = tmp_memoire(2:nb_quant+1);
         Xpp = (tmp_plus - M_k(Grilles(t-1,:),Mu,dt))./(Sigma*sqrt(dt));
         Xmm = (tmp_moins - M_k(Grilles(t-1,:),Mu,dt))./(Sigma*sqrt(dt));
         taille = PHI(Xpp) - PHI(Xmm);
         Poids(t,:) = Poids(t,:)*taille;
         
        
        //Grilles(1,1) = x0; // INITIALISATION seule une seule valeur est possible
       // Poids(1,1) = 1;     // Seul le poids de X0 est non nul et on garde bein la somme de la ligne vaut 1
            [tmp,Poids(t,:)] = Lloyd_1d(nb_quant,init,nb_iter);
            Grilles(t,:) = M_k(Grilles(t-1,:),Mu,dt) + S_k(Grilles(t-1,:),Sigma,dt).*tmp;
            t = t+1;
            //init = gsort(rand(1,nb_quant,"normal"),'g','i');
            //init = Sigma(1)*sqrt(dt)*init + x0;
            
           // disp('grill')
           // disp(Grilles)
      //      disp('poid')
        //    disp(Poids)
          //  disp(' ')
        //  if  t== 2 then
        //       plot(Grilles(t,:)',Poids(t,:)',"<")
        //  end
        //   plot(Grilles(t,:)',Poids(t,:)')
        //    if t== nb_step then
        //       plot(Grilles(t,:)',Poids(t,:)',"*")
          end
        //end
 endfunction




/////////////////////////////           Code Sagna

function GP = Newtonp2(nbre_point,nbre_iter)
//nbre_point=513;
//nbre_iter=50;
stacksize(20000000);
rac2pi=1/(sqrt(%pi * 2));   

a=-0.5;   b=0.5; 
grille = linspace(a,b,nbre_point)';
Moy=zeros(nbre_point,1);   Var=ones(nbre_point,1);
Moy1=Moy(1:nbre_point-1);  Var1=Var(1:nbre_point-1);
for i=1:nbre_iter;
   C = 0.5*(grille(1:nbre_point-1)+grille(2:nbre_point));
   P=cdfnor("PQ",grille,Moy,Var);
   Q=cdfnor("PQ",C,Moy1,Var1);
   expo=rac2pi.*exp(-0.5.*C.^2); expdiag=rac2pi.*exp(-0.5.*grille.^2);
   grillem=grille(1:nbre_point-1) ; grillep=grille(2:nbre_point); 
   gradient=grille.*([Q;1]-[0;Q]) + [expo;0]-[0;expo];
   
   D=-0.25.*(grille-[0;grillem]).*([0;expo])+ 0.25.*(grille-[grillep;0]).*([expo;0])+ [Q;1]-[0;Q]
   D=diag(D);
   V1=-0.25.*([0;expo]).*(grille-[0;grillem]);
   U1=-0.25.*([expo;0]).*([grillep;0]-grille);
   U=diag(U1)*diag(ones(nbre_point-1,1),1);
   V=diag(V1)*diag(ones(nbre_point-1,1),-1);
   hessien=U+D+V;
   A=hessien\gradient;
   grille=grille - A; 
end
 C = 0.5*(grille(1:nbre_point-1)+grille(2:nbre_point));
 Q=cdfnor("PQ",C,Moy1,Var1);
 expo=rac2pi.*exp(-0.5.*C.^2); 
 distor_local = ([0;expo]).*([0;C]-2*grille)+([expo;0]).*(-[C;0]+2*grille)+(1+grille.^2).*([Q;1]-[0;Q]);
 disp(distor_local)
   distor= sqrt(sum(distor_local));
   
   
//C = 0.5*(grille(1:nbre_point-1)+grille(2:nbre_point));
proba_cum=cdfnor("PQ",C,Moy1,Var1);
proba=[proba_cum(1);(proba_cum(2:nbre_point-1)-proba_cum(1:nbre_point-2));(1-proba_cum(nbre_point-1))];

grille(nbre_point+1)=0;
proba(nbre_point+1)=distor;

GP=[grille,proba];
endfunction

//////////////////////////////////////////////////////////////////:::CODE TROP COMMENTER
function [G_t,P_t] = Lloyd_RMQ(G_moins,P_moins,nb_iter,Mu,Sigma,dt)// Recursive Marginal Quantization
compteur = 0;
nb_quant = length(G_moins);
G_t = G_moins;
P_t = ones(1,nb_quant);

for compteur = 1:nb_iter
    // disp("P_")
   // disp(P_t)
    //disp(sum(P_t))
  //  disp("gmoins")
  //  disp(G_moins)
    tmp_memoire = cat(2,-%inf,G_moins,%inf);
    //disp("tmp_memoire")
    //disp(tmp_memoire)
    tmp_memm = (G_t + tmp_memoire(1:nb_quant)); // Xj-1/2 ; au temps k+1;
    //disp("tmp_memm")
    //disp(tmp_memm)
    tmp_memm = 0.5*tmp_memm;
  //  disp("tmp_memm")
    //disp(tmp_memm)
    //disp("tmp_memp")
    tmp_memp = (G_t + tmp_memoire(3:nb_quant+2));// Xj+1/2 ; au temps k+1;
    //disp(tmp_memp)
    tmp_memp = 0.5*tmp_memp;
    //disp("tmp_memp")
    //disp(tmp_memp)
    //-repmat((M_k(G_moins,Mu,dt))',1,nb_quant);// matrice nn des -M_k(Xi) constante sur une meme LIGNE!!!
    Xmm = repmat(tmp_memm,nb_quant,1) - repmat((M_k(G_moins,Mu,dt))',1,nb_quant)  // Xj-1/2 -M_k(Xi) matrice nn
    //disp("m_k g moins")
    //disp(repmat((M_k(G_moins,Mu,dt))',1,nb_quant))
    Xpp = repmat(tmp_memp,nb_quant,1) - repmat((M_k(G_moins,Mu,dt))',1,nb_quant)  //  Xj+1/2 -M_k(Xi) matrice nn
    //disp("xpp")
    //disp(Xpp)
    Xmm = Xmm/(Sigma*sqrt(dt))  // V(i)j-1/2  matrice nn
  //  disp("xmm")
    //disp(Xmm)    
    Xpp = Xpp/(Sigma*sqrt(dt))  //  V(i)j+1/2 matrice nn
    //disp("xpp")
    //disp(Xpp)
    DPhi = PHI(Xpp) - PHI(Xmm) // taille des vj+1/2;i
    //disp("dphi")
    //disp(DPhi)
    //disp(sum(DPhi,1))
    P_t = P_moins*DPhi
    //disp("Pt")
    //disp(P_t)
    A = P_moins.*M_k(G_moins,Mu,dt)*(PHI(Xpp) - PHI(Xmm));
    B = P_moins*Esperance_Gauss(Xpp,Xmm);
    B = Sigma*sqrt(dt)*B;
    B = B.*G_moins;
    G_t = A + B;
    G_t = G_t./P_t;
    //disp("grill")
    //disp(G_t)
    disp(compteur, " iteration ")
end
disp(G_t , P_t, "poids , grille")
G_t = G_t;
P_t = P_t;
endfunction

