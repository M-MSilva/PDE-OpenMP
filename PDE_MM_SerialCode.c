#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <time.h>
#define tolerance 1e-9


//Definindo funcoes vai ate a linha 62
//Funcao Inicial Da EDP
double g( double x, double y)
{
    return x*exp(y);
}
//condicoes de contorno para a EDP
double fc( double x, double y)
{
    return x;
}
double fd( double x, double y)
{
    return (exp(1)*x) ;
}
double fa( double x, double y)
{
 return 0;
}
double fb( double x, double y)
{
    return 2*exp(y);
}
// alocar vetor e preencher
double *iniciaVetor(int z,double o, double p) 
{
  double *v;
  int l;
  
 v = (double*)malloc(z*sizeof(double)+1);

 for(l=1;l<z;l++)
 {
    v[l] = p + l*o;
    //o faz o papel do h,k, enquanto p faz o papel do a,b
 }
return(v); 
}
double **IniciaMatriz(int m, int n)
{
  double **M;
  int i,j;
  
  M = malloc(n*sizeof(double *) + 1);
  for ( i = 0; i < n; i++ ) M[i] = malloc(m*sizeof(double) + 1);
  
  for ( i = 1; i < n; i++ )
  {
    for ( j = 1; j < m; j++ )
    {
       M[i][j]=0;
    }  
  }
  return M;
}
void main ()
{
   
    double z,k,h,mu, Norma,lambda ,*x ,*y ,**omega,a,b,c,d;
    int itera,i, j,m,n;
    a =0;
    b =2;
    c =0;
    d =1;

    
    printf ("Programa:  Resolvendo EDP de Poisson por Diferenciais Finitas v.07 \n");
    printf ("Aluno :  Marcos Matheus de Paiva Silva\n");
    printf ("Data: 24/02/2021\n");
    printf ("Executando.....\n");


    // a quantidade de particoes 2D
    n = 760;
    m = 760;

    //tamanho dos intervaloros das m,n partes
    h = (b - a)/(1.0*n);
    k = (d - c)/(1.0*m);

    // Alocar a matriz omega que faz uma aproximacao para f(x,y)
    //e posteriormente preencher omega com zero para evitar erros
    omega = IniciaMatriz(m,n);

    //alocar vetor para x e preencher vetor correspondente
    x = iniciaVetor(n,h,a);
    //faz similar com y
    y = iniciaVetor(m,k,c);

    //Para discretizar, usa-se as constantes lambda mu
    lambda = h*h/(1.0*k*k);
    mu = 2*(1 + lambda );
    itera=0;

    //repeticao feita ate a norma de f(x,y) ser menor que tolernacia
    do{
        
        //na ultima coluna da matriz omega,efetua o primeiro elemento
        z = (-(h*h)*g(x[1] ,y[m -1]) +fa(a,y[m-1]) + lambda*fd(x[1] , d)+ lambda*omega[1][m -2]+
        omega[2][m -1])/mu;

        Norma = fabs(z - omega[1][m -1]) ;//primeira norma do loop

        omega[1][m - 1] = z;

        //na ultima coluna da matriz omega,efetua os elementos do centro
        for(i=2;i<=n -2;i++)
        { 
            z=(- (h*h)*g(x[i],y[m -1]) + lambda *fd (x[i],d)+omega[i -1][m -1] + omega[i +1][m -1]+ lambda*omega[i][m -2])/mu;
            if( fabs (omega[i][m -1] -z)> Norma)
                Norma = fabs(omega[i][m -1] -z); //recalcula a norma
            omega[i][m -1]= z;
        }

        /*na ultima coluna da matriz omega,determina o ultimo elemento */

        z=(- (h*h)*g(x[n -1] ,y[m -1]) +fb(b,y[m -1]) + lambda*fd(x[n -1] ,d)+ omega[n -2][m -1] + lambda*omega[n -1][m -2])/mu;
        if( fabs(omega[n -1][m -1] -z)>Norma )
            Norma = fabs(omega[n -1][m -1] -z);//recalcula a norma
        omega[n -1][m -1]= z;

        /*computa o centro da matriz,iniciando pela penultima coluna
        ,finalizando apenas na segunda coluna da matriz*/
        for(j=m-2;j>1;j --)
        {

            //efetua o primeira elemento da matriz do centro analizada
            z=(-(h*h)*g(x[1] ,y[j])+fa(a,y[j])+ lambda*omega[1][j +1]+ lambda*omega[1][j -1]+ omega[2][ j])/mu;

            if( fabs (omega[1][ j]-z) > Norma )
                Norma = fabs (omega[1][ j]-z);//recalcula a norma
            omega[1][j]=z;


            /*na coluna analizada,calcula o centro da coluna
             de omega[2][j] ate omega[N-2][ j]*/
            for (i=2;i<=n-2;i++)
            {
                z=(-(h*h)*g(x[i],y[j]) + omega[i -1][ j]+ lambda*omega[i][j +1]+ omega[i +1][ j]+ lambda*omega[i][j -1])/mu;
                if( fabs (omega[i][j]-z)>Norma )
                    Norma= fabs(omega[i][j]-z); //recalcula a norma
                omega[i][j]=z;
            }

            /*Na coluna analizada,computa o ultimo elemento de omega*/
            z=(-(h*h)*g(x[n -1] ,y[j]) + fb(b,y[j]) + omega[n -2][ j]+ lambda*omega[n -1][ j +1]+ lambda*omega[n -1][j -1])/mu;
            if( fabs (omega[n -1][ j]-z)>Norma )
                Norma = fabs (omega[n -1][ j]-z);//recalcula a norma
            omega[n -1][ j]=z;
        }//finaliza o centro da matriz

        /* Na primeira coluna da matriz omega, efetua seu primeiro elemento*/
        z=(- (h*h)*g(x[1] ,y[1]) +fa(a,y[1]) + lambda*fc(x[1] ,c)+ lambda*omega[1][2]+ omega[2][1])/mu;
        if( fabs (omega[1][1] - z)>Norma)
            Norma = fabs (omega[1][1] - z);//recalcula a norma
        omega[1][1]= z;


        /*Na primeira coluna da matriz, computa seus elementos do centro
        de  omega[2][1] ate omega[N -2][1]*/
        for(i=2;i<=n-2;i++)
        {
            z=(-(h*h)*g(x[i],y [1]) + lambda*fc(x[i],c) + omega[i -1][1] + lambda*omega[i ][2]+ omega[i +1][1])/mu;

        if( fabs (omega[i][1] -z)>Norma )
            Norma = fabs (omega[i][1] -z);//recalcula a norma
        omega[i ][1]= z;
        }

        /* Na primeira coluna, efetua  seu ultimo elemento*/
        z=(- (h*h)*g(x[n -1] ,y [1]) +fb(b,y[1]) + lambda*fc(x[n -1] ,c)+omega[n -2][1] + lambda*omega[n -1][2])/mu;
        if( fabs (omega[n -1][1] - z)>Norma ) 
            Norma = fabs(omega[n -1][1] - z);//recalcula a norma
        omega[n -1][1]= z;

        //conta as iteracoes
        itera++;

    }while(Norma>=tolerance);
    
    printf ("Iteracoes: %d \n\n", itera);

    //
    /*mostra os resultados em um arquivo*/
    FILE *fpt; 

    fpt = fopen("MMedpPoisson630.dat","w");

    //joga os dados em edpPoisson.dat
    for(i=1;i<=n-1;i++)
    {
        for(j=1;j<=m-1;j++)
        {
            fprintf(fpt ,"%.9lf\t %.9lf\t %.9lf \n", x[i],y[j],omega[i][j]);
        }
    }
  
    //fecha o arquivo
    fclose (fpt);

    //desalocar
    free(x);
    free(y);
    free(omega);


}



