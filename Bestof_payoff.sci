function [res] = Bestof_payoff(g1,g2,K)
    x = max(g1,g2)
    res = max(x-K,zeros(x));
endfunction

function [res] = Bestof_payoff_all(p1,p2,K)
n1 = length(p1);
n2 = length(p2);
res = ones(1,n1*n2);
 res = [];
 for i = 1:n2
    tmp3 = Bestof_payoff(p2(i),p1,K);   
    res = cat(1,res,tmp3);
end
res;
endfunction 


function [res] = Worstof_payoff(g1,g2,K)
    x = min(g1,g2)
    res = max(K-x,zeros(x));
endfunction

function [res] = Worstof_payoff_all(p1,p2,K)
n1 = length(p1);
n2 = length(p2);
res = ones(1,n1*n2);
 res = [];
 for i = 1:n2
    tmp3 = Worstof_payoff(p2(i),p1,K);   
    res = cat(1,res,tmp3);
end
res;
endfunction 


function [res] = prod_tenso_i(p1,p2)
n1 = length(p1);
n2 = length(p2);
res = ones(1,n1*n2);
 res = [];
 for i = 1:n2
    tmp3 = p2(i)*p1;   
    res = cat(1,res,tmp3);
end
res;
endfunction 

function [res] = prod_tenso_j(p1,p2)
n1 = length(p1);
n2 = length(p2);
res = ones(1,n1*n2);
 res = [];
 for i = 1:n2
    tmp3 = p2(i)*p1;   
    res = cat(2,res,tmp3);
end
res;
endfunction 


function [res] = Call_Bestof_BS_RMQ_2d(nb_quant,nb_step,nb_iter,nbre_iter1,mu1,mu2,sigma1,sigma2,x0,y0,K,T)//2 uncorrelated asset

[g1,p1] = NewtonBS(nb_quant,nb_step,nb_iter,nbre_iter1,mu1,sigma1,x0,T);   
[g2,p2] = NewtonBS(nb_quant,nb_step,nb_iter,nbre_iter1,mu2,sigma2,y0,T);
//pause
[g,p] = Prod_2_Quantif(g1(:,nb_step+1),p1(:,nb_step+1),g2(:,nb_step+1),p2(:,nb_step+1));
res = p'*Bestof_payoff(g(:,1),g(:,2),K);
endfunction

function [resf] = Bestof_BS_US_RMQ_2d(nb_quant,nb_step,nbre_iter,nbre_iter1,mu1,mu2,sigma1,sigma2,x0,y0,K,T)
dt = T/nb_step;
[g1,p1] = NewtonBS(nb_quant,nb_step,nbre_iter,nbre_iter1,mu1,sigma1,x0,T);   
[g2,p2] = NewtonBS(nb_quant,nb_step,nbre_iter,nbre_iter1,mu2,sigma2,y0,T);
// CENTRES DES CELLULES
GPUnDemi1 = zeros(nb_quant,nb_step+1)
GMUnDemi1 = zeros(nb_quant,nb_step+1)
GPUnDemi2 = zeros(nb_quant,nb_step+1)
GMUnDemi2 = zeros(nb_quant,nb_step+1)
GPUnDemi1 = 0.5*( g1 + cat(1,g1(2:nb_quant,:),%inf*ones(1,nb_step+1)))//frontiere sup
GMUnDemi1 = 0.5*( g1 + cat(1,-%inf*ones(1,nb_step+1),g1(1:nb_quant-1,:)))//frontiere inf
GPUnDemi2 = 0.5*( g2 + cat(1,g2(2:nb_quant,:),%inf*ones(1,nb_step+1)))//frontiere sup
GMUnDemi2 = 0.5*( g2 + cat(1,-%inf*ones(1,nb_step+1),g2(1:nb_quant-1,:)))//frontiere inf
//pause
// ESPERANCE CONDITIONNELLE
    espeG = zeros(g1);
    espeG = zeros(g2);
    a1 = zeros(nb_quant,nb_quant);
    b1 = zeros(nb_quant,nb_quant);
    a2 = zeros(nb_quant,nb_quant);
    b2 = zeros(nb_quant,nb_quant);
    //pause   
    p = nb_quant*nb_quant;
    res = zeros(p,nb_step+1);
    res(:,nb_step+1) = exp(-mu1*T)*Bestof_payoff_all(g1(:,nb_step+1),g2(:,nb_step+1),K);
for t = nb_step:-1:1 // recurrence rétrograde
    grillefut1 = g1(:,t+1);//g1_t+1
    grillefut2 = g2(:,t+1);//g2_t+1    
    
    for i = 1:nb_quant // parcours du support de X_{k+1}
        //pause
        a1(:,i) = (GPUnDemi1(:,t+1) - M_k(g1(i,t),mu1,dt))./S_k(g1(i,t),sigma1,dt);
        b1(:,i) = (GMUnDemi1(:,t+1) - M_k(g1(i,t),mu1,dt))./S_k(g1(i,t),sigma1,dt);   
        a2(:,i) = (GPUnDemi2(:,t+1) - M_k(g2(i,t),mu1,dt))./S_k(g2(i,t),sigma2,dt);
        b2(:,i) = (GMUnDemi2(:,t+1) - M_k(g2(i,t),mu1,dt))./S_k(g2(i,t),sigma2,dt);                   
    end
    //pause
    mat_proba = zeros(nb_quant,nb_quant,2);
    mat = zeros(p,p);
    mat_proba(:,:,1) = PHI(a1) - PHI(b1);
    mat_proba(:,:,2) = PHI(a2) - PHI(b2);
    mat
    for k = 1:nb_quant
        mat(:,k) = prod_tenso_i(mat_proba(:,k,1),mat_proba(:,k,2)) 
    end
    //(PHI(a1) - PHI(b1)./(PHI(a2) - PHI(b2));
    //tmpres = (res(:,t+1))'*(PHI(a1) - PHI(b1)./(PHI(a2) - PHI(b2));
    tmpres = (res(:,t+1))'*mat;    
    //disp("a")
    //pause
   //     disp("b")
   // espeG(:,t) = tmpres'; 
    res(:,t) = max(tmpres',exp(-mu1*t*dt)*Bestof_payoff_all(g1(:,t),g2(:,t),K));
        //pause
end
res;
resf = res(1,1);
endfunction





function [resf] = Worstof_BS_US_RMQ_2d(nb_quant,nb_step,nbre_iter,nbre_iter1,mu1,mu2,sigma1,sigma2,x0,y0,K,T)
dt = T/nb_step;
[g1,p1] = NewtonBS(nb_quant,nb_step,nbre_iter,nbre_iter1,mu1,sigma1,x0,T);   
[g2,p2] = NewtonBS(nb_quant,nb_step,nbre_iter,nbre_iter1,mu2,sigma2,y0,T);
// CENTRES DES CELLULES
GPUnDemi1 = zeros(nb_quant,nb_step+1)
GMUnDemi1 = zeros(nb_quant,nb_step+1)
GPUnDemi2 = zeros(nb_quant,nb_step+1)
GMUnDemi2 = zeros(nb_quant,nb_step+1)
GPUnDemi1 = 0.5*( g1 + cat(1,g1(2:nb_quant,:),%inf*ones(1,nb_step+1)))//frontiere sup
GMUnDemi1 = 0.5*( g1 + cat(1,-%inf*ones(1,nb_step+1),g1(1:nb_quant-1,:)))//frontiere inf
GPUnDemi2 = 0.5*( g2 + cat(1,g2(2:nb_quant,:),%inf*ones(1,nb_step+1)))//frontiere sup
GMUnDemi2 = 0.5*( g2 + cat(1,-%inf*ones(1,nb_step+1),g2(1:nb_quant-1,:)))//frontiere inf
//pause
// ESPERANCE CONDITIONNELLE
    espeG = zeros(g1);
    espeG = zeros(g2);
    a1 = zeros(nb_quant,nb_quant);
    b1 = zeros(nb_quant,nb_quant);
    a2 = zeros(nb_quant,nb_quant);
    b2 = zeros(nb_quant,nb_quant);
    //pause   
    p = nb_quant*nb_quant;
    res = zeros(p,nb_step+1);
    res(:,nb_step+1) = exp(-mu1*T)*Worstof_payoff_all(g1(:,nb_step+1),g2(:,nb_step+1),K);
for t = nb_step:-1:1 // recurrence rétrograde
    grillefut1 = g1(:,t+1);//g1_t+1
    grillefut2 = g2(:,t+1);//g2_t+1    
    
    for i = 1:nb_quant // parcours du support de X_{k+1}
        //pause
        a1(:,i) = (GPUnDemi1(:,t+1) - M_k(g1(i,t),mu1,dt))./S_k(g1(i,t),sigma1,dt);
        b1(:,i) = (GMUnDemi1(:,t+1) - M_k(g1(i,t),mu1,dt))./S_k(g1(i,t),sigma1,dt);   
        a2(:,i) = (GPUnDemi2(:,t+1) - M_k(g2(i,t),mu2,dt))./S_k(g2(i,t),sigma2,dt);
        b2(:,i) = (GMUnDemi2(:,t+1) - M_k(g2(i,t),mu2,dt))./S_k(g2(i,t),sigma2,dt);                   
    end
    //pause
    mat_proba = zeros(nb_quant,nb_quant,2);
    mat = zeros(p,p);
    mat_proba(:,:,1) = PHI(a1) - PHI(b1);
    mat_proba(:,:,2) = PHI(a2) - PHI(b2);

    for k = 1:nb_quant
     mat(:,k) = prod_tenso_i(mat_proba(:,k,1),mat_proba(:,k,2)) 
    end
    //(PHI(a1) - PHI(b1)./(PHI(a2) - PHI(b2));
    //tmpres = (res(:,t+1))'*(PHI(a1) - PHI(b1)./(PHI(a2) - PHI(b2));
    tmpres = (res(:,t+1))'*mat;    
    //pause
    //espeG(:,t) = tmpres';
    res(:,t) = max(tmpres',exp(-mu1*t*dt)*Worstof_payoff_all(g1(:,t),g2(:,t),K));
        //pause
end
resf = res(1,1);
endfunction


