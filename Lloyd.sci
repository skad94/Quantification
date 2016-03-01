function [res] = Esperance_Gauss(sup,inf)
res = (exp(-0.5.*inf.^2) -  exp(-0.5.*sup.^2)) /(sqrt(%pi * 2));
endfunction


function [res] = PHI(x)// fonction de repartition
 //   p = length(x);
 //   mu = zeros(1,p);
  // sigma = ones(1,p);
    res = cdfnor("PQ",x,zeros(x),ones(x));
endfunction

function [res] = M_k(x,Mu,dt)// drift aussi egale à l'esperance de la marginal
    res = x + Mu.*x*dt;
endfunction

function [res] = S_k_CEV(x,Nu,Delta,dt)// drift aussi egale à l'esperance de la marginal
    res = sqrt(dt)*(Nu* x.^(Delta + 1))/(sqrt(1+x.^2));
endfunction 
   
function [res] = S_k(x,Sigma,dt)// diffusion  aussi egale à l'ecart type de la marginal
        res = Sigma.*x*sqrt(dt);
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
    Grilles = espe ./(Poids + 0.000000001);
   end

    Grilles = Grilles;
    Poids = Poids;
endfunction

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


function [G_t,P_t] = Lloyd_RMQ(G_moins,P_moins,Mu,Sigma,dt)// Recursive Marginal Quantization
compteur = 0;
nb_quant = length(G_moins);
G_t = G_moins;
P_t = ones(1,nb_quant);
//  disp(G_t, "gt")
//disp(G_moins, "gmoin")
    tmp_memoire = cat(2,-%inf,G_t,%inf);
    tmp_memm = (G_t + tmp_memoire(1:nb_quant)); // Xj-1/2 ; au temps k+1;
    tmp_memm = 0.5*tmp_memm;
    tmp_memp = (G_t + tmp_memoire(3:nb_quant+2));// Xj+1/2 ; au temps k+1;
    tmp_memp = 0.5*tmp_memp;
    //-repmat((M_k(G_moins,Mu,dt))',1,nb_quant);// matrice nn des -M_k(Xi) constante sur une meme LIGNE!!!
    Xmm = repmat(tmp_memm,nb_quant,1) - repmat((M_k(G_moins,Mu,dt))',1,nb_quant)  // Xj-1/2 -M_k(Xi) matrice nn
    Xpp = repmat(tmp_memp,nb_quant,1) - repmat((M_k(G_moins,Mu,dt))',1,nb_quant)  //  Xj+1/2 -M_k(Xi) matrice nn
    //disp(Xpp, "xpp")
    Xmm = Xmm/(Sigma*sqrt(dt))  // V(i)j-1/2  matrice nn
    Xpp = Xpp/(Sigma*sqrt(dt))  //  V(i)j+1/2 matrice nn

    DPhi = PHI(Xpp) - PHI(Xmm) // taille des vj+1/2;i
    
    P_t = P_moins*DPhi;
    P_t = P_t;
    // grilles
    A = P_moins.*M_k(G_moins,Mu,dt)*(PHI(Xpp) - PHI(Xmm));
    B = P_moins*Esperance_Gauss(Xpp,Xmm);
    B = Sigma*sqrt(dt)*B;
    B = B.*G_moins;
    G_t = A + B;
   // disp(P_t, "pt")
    G_t = G_t./(P_t+0.000001);
   //disp(G_t, " grill ")

//disp(G_t , P_t, "poids , grille")
G_t = G_t;
P_t = P_t;
endfunction

function [G_t,P_t] = Lloyd_RMQ_CEV(G_moins,P_moins,Mu,Nu,Delta,dt)// Recursive Marginal Quantization
compteur = 0;
nb_quant = length(G_moins);
G_t = G_moins;
P_t = ones(1,nb_quant);
//  disp(G_t, "gt")
//disp(G_moins, "gmoin")
    tmp_memoire = cat(2,-%inf,G_t,%inf);
    tmp_memm = (G_t + tmp_memoire(1:nb_quant)); // Xj-1/2 ; au temps k+1;
    tmp_memm = 0.5*tmp_memm;
    tmp_memp = (G_t + tmp_memoire(3:nb_quant+2));// Xj+1/2 ; au temps k+1;
    tmp_memp = 0.5*tmp_memp;
    //-repmat((M_k(G_moins,Mu,dt))',1,nb_quant);// matrice nn des -M_k(Xi) constante sur une meme LIGNE!!!
    Xmm = repmat(tmp_memm,nb_quant,1) - repmat((M_k(G_moins,Mu,dt))',1,nb_quant)  // Xj-1/2 -M_k(Xi) matrice nn
    Xpp = repmat(tmp_memp,nb_quant,1) - repmat((M_k(G_moins,Mu,dt))',1,nb_quant)  //  Xj+1/2 -M_k(Xi) matrice nn
    den = ones(G_moins)./sqrt(1+(G_moins).^2)//Nu*(G_moins.^(Delta+1))
    Xmm = Xmm./S_k_CEV(G_moins,Nu,Delta,dt)  // V(i)j-1/2  matrice nn
    Xpp = Xpp./S_k_CEV(G_moins,Nu,Delta,dt)  //  V(i)j+1/2 matrice nn

    DPhi = PHI(Xpp) - PHI(Xmm) // taille des vj+1/2;i
    
    P_t = P_moins*DPhi;
    P_t = P_t;
    // grilles
    A = P_moins.*M_k(G_moins,Mu,dt)*(PHI(Xpp) - PHI(Xmm));
    B = P_moins*Esperance_Gauss(Xpp,Xmm);
    //B = Sigma*sqrt(dt)*B;
    B = B*S_k_CEV(G_moins,Nu,Delta,dt);
    B = B.*G_moins;
    G_t = A + B;
   // disp(P_t, "pt")
    G_t = G_t./(P_t+0.000001);
   //disp(G_t, " grill ")

//disp(G_t , P_t, "poids , grille")
a_t = G_t';
b_t = P_t';
//plot(a_t,b_t)
endfunction


function [Grilles,Poids] = Lloyd2BS(nb_quant,nb_step,nb_iter,x0,Mu,Sigma,T,x1g,x1p)
dt = T/nb_step;
Grilles = ones(nb_step+1,nb_quant);
Poids = ones(nb_step+1,nb_quant);
Grilles (1,:)= x0*Grilles (1,:); // x0 et 
init = gsort(rand(1,N,"normal"),'g','i');//loi normal naif
//[x1g,x1p] = Lloyd_1d(nb_quant,init,nb_iter);
Grilles(2,:) = M_k(x0,Mu,dt)+S_k(x0,Sigma,dt)*x1g;
Poids(2,:) = x1p;
plot(Grilles(2,:),Poids(2,:),">")
plot(Grilles(2,:),Poids(2,:))
for t = 3:nb_step+1

    [tmp,Poids(t,:)] = Lloyd_RMQ(Grilles(t-1,:),Poids(t-1,:),Mu,Sigma,dt);
    Grilles(t,:) = M_k(Grilles(t-1,:),Mu,dt) + S_k(Grilles(t-1,:),Sigma,dt).*x1g;
  //  disp(Grilles,Poids, "poids , grilles ")
   // plot(Grilles(t,:),Poids(t,:))
end    
plot(Grilles(nb_step,:),Poids(nb_step,:),"*")
plot(Grilles(nb_step,:),Poids(nb_step,:))
Grilles;
Poids;  
endfunction    


 function [Grilles,Poids] = Lloyd2CEV(nb_quant,nb_step,nb_iter,x0,Mu,Nu,Delta,T)
dt = T/nb_step;
Grilles = ones(nb_step+1,nb_quant);
Poids = ones(nb_step+1,nb_quant);
Grilles (1,:)= x0*Grilles (1,:); // x0 et 
init = gsort(rand(1,N,"normal"),'g','i');//loi normal naif
moy = x0 + Mu*x0*dt;
var = (Nu*sqrt(dt/(1+x0^2))*x0.^(Delta+1));
init = moy + var*init;
[x1g,x1p] = Lloyd_1d(nb_quant,init,nb_iter);
Grilles(2,:) = M_k(x0,Mu,dt)+S_k_CEV(x0,Nu,Delta,dt).*x1g;
Poids(2,:) = x1p;

for t = 3:nb_step+1

    [tmp,Poids(t,:)] = Lloyd_RMQ_CEV(Grilles(t-1,:),Poids(t-1,:),Mu,Nu,Delta,dt);
    Grilles(t,:) = M_k(Grilles(t-1,:),Mu,dt) + S_k_CEV(Grilles(t-1,:),Nu,Delta,dt).*x1g;

end    

Grilles = Grilles';
Poids = Poids';  
clf;
plot(Grilles(:,2),Poids(:,2),"<",Grilles(:,nb_step+1),Poids(:,nb_step+1),"*")
plot(Grilles(1:nb_quant,2:nb_step+1),Poids(1:nb_quant,2:nb_step+1))
endfunction    


/////////////       §§§§§§§§§§§        1111111111   !!!!!!!!!!!!!!!!!!!!!!!!!!!
        //      2D2D2D2D22D2D2D2D2D2D22D2D2D
       // QUANTIFICATION PRODUIT
function [Grille,Poids] = Prod_2_Quantif(g1,p1,g2,p2)
n1 = length(g1);
n2 = length(g2);
Grille = ones(2,n1*n2);
Poids = ones(1,n1*n2);
Grille (1,:) = repmat(g1,1,n2);
 tmp2 = [];
 Poids = [];
 for i = 1:n2
    tmp1 = repmat(g2(i),1,n1);  
    tmp2 = cat(2,tmp2,tmp1);
    tmp3 = p2(i)*p1   
    Poids = cat(2,Poids,tmp3);
end
Grille(2,:) = tmp2;
endfunction           
       
function [Grilles,Poids] = Lloyd2BS_2d(nb_quant,nb_step,nb_iter,x0,y0,Mu,Sigma,Mu0,Sigma0,T,Zg,ZP)// uncorrelated
dt = T/nb_step;
Grilles = ones(nb_step+1,nb_quant,2);
Grilles (1,:,1)= x0*Grilles (1,:,1); // x0 et 
Grilles (1,:,2)= y0*Grilles (1,:,2);
init = gsort(rand(1,N,"normal"),'g','i');//loi normal naif
//disp(Grilles(2,:,1), "grille")
//disp(M_k(x0,Mu,dt), "mk(x0,mu,dt)")
//disp(S_k(x0,Sigma,dt), "sk(x0...)")
//disp(Zg, "zg")
//pause
Grilles(2,:,1) = M_k(x0,Mu,dt)+S_k(x0,Sigma,dt)*Zg

Poids(2,:,1) = ZP;

for t = 3:nb_step+1

    [tmp,Poids(t,:,1)] = Lloyd_RMQ(Grilles(t-1,:,1),Poids(t-1,:,1),Mu,Sigma,dt);
    [tmp2,Poids(t,:,2)] = Lloyd_RMQ(Grilles(t-1,:,1),Poids(t-1,:,1),Mu0,Sigma0,dt);    
    Grilles(t,:,1) = M_k(Grilles(t-1,:,1),Mu,dt) + S_k(Grilles(t-1,:,1),Sigma,dt).*Zg;
    Grilles(t,:,2) = M_k(Grilles(t-1,:,2),Mu,dt) + S_k(Grilles(t-1,:,2),Sigma,dt).*Zg;    
  //  disp(Grilles,Poids, "poids , grilles ")
    //plot(Grilles(t,:),Poids(t,:))
end 
[Grilles,Poids] = Prod_2_Quantif(Grilles(t,:,1),Poids(t,:,1),Grilles(t,:,2),Poids(t,:,2));  

Grilles;
Poids;  
endfunction
       
function [res] = Basket_Call(g12,p12,weig1,weig2,K)// uncorrelated 

[nb_iter,nb_quant] = size(g12);
disp(g12, "g12")
g12(1,:) = weig1*g12(1,:) //grille pondere du sous jacent 0
g12(2,:) = weig2*g12(2,:) //grille pondere du sous jacent 1 
disp(g12, "g12")
tmp = sum(g12,1);
disp(tmp,"sum")
tmp = tmp - K;
//disp(tmp, "s-k")
tmp = max(zeros(tmp),tmp);
disp(tmp, "(s-k )+")
res = tmp*p12';
//disp(res, "E(s-k)")
disp(res, "res =")
res = exp(-0.02)*res
endfunction    
       
       
       
       
       
