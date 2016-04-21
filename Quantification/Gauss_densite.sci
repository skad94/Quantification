double Gauss_densite(double x)
function [res] = Gauss_densite(nb_quant,init,nb_iter)
{
   // printf("g\n");
    double res;
    res = PI;
   // printf("pi = %f\n",res);
    res = exp(-x*x/2)/sqrt(2*res);
    //printf("gauss(x) = %f\n",res);
    return res;
}
