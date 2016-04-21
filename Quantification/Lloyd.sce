exec('C:\Users\stéphane\Documents\Scilab\Quantification\Lloyd.sci');
//exec('C:\Users\stéphane\Documents\Scilab\Quantification\Esperance_Expo.sci');
//exec('C:\Users\stéphane\Downloads\Newtonp2BS.sci');
N = 15;
nb_step = 250;
nb_iter = 3000;
Mu = 0.02;
Mu0 = Mu;
Sigma = 0.4;
Sigma2 = 0.25;
x0 = 100;
y0 = 100;
T = 1;

initg = gsort(rand(1,N,"normal"),'g','i');//quantif naif
initg = initg;//*Sigma + Mu;
//[a,b] = Lloyd_1d(N,initg,nb_iter);
//initg = linspace(-1,1,N);
//initp = zeros(1,N); initp(1)=1;

//Lloyd_RMQ(G_moins,P_moins,nb_iter,Mu,Sigma,dt)
    //Lloyd_RMQ(initg,initp,5,Mu,Sigma,1/100)
 //[g,p]=Lloyd2BS(N,120,3000,86.3,0.15,0.05,1);
 
 
 // SYNT Lloyd2CEV(nb_quant,nb_step,nb_iter,x0,Mu,Nu,Delta,T)
// [gc,pc] = Lloyd2CEV(N,nb_step,nb_iter,100,0.15,0.7,0.5,0.5)


//Lloyd2BS(nb_quant,nb_step,nb_iter,x0,Mu,Sigma,T,x1g,x1p)
     //    Lloyd2BS(N,nb_step,nb_iter,85,0.03,0.02,1,a,b);

 //LLoyd2BS  syntaxe (nb_quant,nb_step,nb_iter,Mu,Sigma,x0,init,T) et il faut dessiner les COLONNES ET NON LES LIGNES!!!!!!
//[a,b] = Lloyd_1d(N,initg,nb_iter);

 //[Grilles,Poids] = Lloyd2BS(N,nb_step,nb_iter,x0,Mu,Sigma,T,a,b)// ATTENTION COLONNES != LIGNES
//ex = Prod_2_Quantif(g);   
//[g2,p2] = NewtonBS(N,20,15,15,0.2,0.1,5,1); Lloyd2BS(nb_quant,nb_step,nb_iter,x0,Mu,Sigma,T)

 //Lloyd_expo_1d(nb_quant,init,nb_iter,intensite)
 //[a,b] = Lloyd_expo_1d(N,init,N,intensite);
prix_n = zeros(100,1);
for i = 15:115
    initg = gsort(rand(1,i,"normal"),'g','i');//quantif naif
    [a,b] = Lloyd_1d(i,initg,nb_iter); 
    prix_n(i-14) = CallBS(a,b,100,100,0.04,0.05,1);
    disp(prix_n(i-14),i,  "prix ")
end

//[Grilles,Poids] = Lloyd_1d(nb_quant,init,nb_iter)
// [a,b] = Lloyd_1d(N,initg,nb_iter);
// SYNTAXE  GP = Newtonp2(nbre_point,nbre_iter)
//[G,P] = Newtonp2(N,50)

//Syntaxe [Grilles,Poids] = Lloyd2BS_2d(nb_quant,nb_step,nb_iter,x0,y0,Mu,Sigma,Mu0,Sigma0,T)// uncorrelated
//[g,p] = Lloyd2BS_2d(N,nb_step,nb_iter,x0,y0,Mu,Sigma,Mu0,Sigma2,T,a,b);
//Syntaxe   res = Basket_Call(g1,p1,g2,p2,weig1,weig2,K)
//prix_basket_call = Basket_Call(g,p,0.5,0.5,100)
