#include <iostream>
#include <cmath>
#include <numeric>
#include <gmp.h>
#include <gmpxx.h>
#include <stdio.h>
using namespace std;

void buildup1 (size_t k, size_t r, mpz_t n, mpz_t ***M); //10.3 step 1: building up subproduct tree
void buildup2 (mpz_t *u,size_t k, size_t r, mpz_t n, mpz_t ***M);// 10.3 building up subproduct tree for step 2
void multi (mpz_t *a, mpz_t * b, mpz_t * ans, size_t deg, mpz_t n);//naive way for polynomial multiplication
void rem(mpz_t *a,mpz_t *b,mpz_t *&rem,size_t dA,size_t dB,size_t &dR, mpz_t n);//remainder of two polynomials a/b
void fastMultiEval(mpz_t *f,size_t k, size_t r,mpz_t * fu,mpz_t ***M, mpz_t n);//10.7 step 2: fast multipoint evaluation
void godown(mpz_t *f,mpz_t *u,mpz_t ***M,size_t deg,mpz_t *fu,size_t k,mpz_t n,size_t start); //10.5 going down subproduct tree

  int main()
{
  mpz_t n;
  size_t b,c,k;
    //set n
    mpz_init_set_str(n,"9797",10);
    //upperbound b
    b=mpz_get_ui(n);
    b=ceil(sqrt(b));
    //choose c to be the smallest power of 2 that is greater than sqrt(b). And k=log2(c)
    k = ceil(log2(ceil(sqrt(b))));
    c=pow(2,k);
    gmp_printf("n: %Zd",n);
    cout<<" b: "<<b<<" c: "<<c<<" k: "<<k<<endl;

    //3D array
    mpz_t ***M1;
    mpz_t ***M2;
    size_t  pow2_i,pow2_ki; 
    M1 = new mpz_t  **[k+1]; //for step 1
    M2 = new mpz_t **[k+1]; //for step 2
  
    for(size_t i = 0; i<=k; i++)
    {
      pow2_i =pow(2,i);
      pow2_ki=c/pow2_i;
      M1[i]=new mpz_t*[pow2_ki];
      M2[i]=new mpz_t*[pow2_ki];
      for(size_t j = 0; j<pow2_ki;j++)
	{
       
	  M1[i][j] = new mpz_t[pow2_i+1];
	  M2[i][j] = new mpz_t[pow2_i+1];
	}
    }
  

    //step1: 10.3
    buildup1(k,c,n,M1);
     cout<<"buildup1  result: f"<<endl;
  for(int i = 0; i<=c;i++)
    { cout<<"f["<<i<<"]:";
      gmp_printf("%Zd ",M1[k][0][i]);}
  cout<<endl;


  //step2: 10.7
  mpz_t *g;
  g = new mpz_t [c];
  fastMultiEval(M1[k][0],k,c,g,M2,n);

       cout<<"fastMultiEval result: g"<<endl;
  for(int i = 0; i<c;i++)
    { cout<<"g["<<i<<"]:";
      gmp_printf("%Zd ",g[i]);}
  cout<<endl;



  //final step: find the least prime factor of n
  mpz_t *gcd;
  gcd = new mpz_t [c];
  mpz_init(gcd[0]);
  mpz_gcd(gcd[0],g[0],n);
  mpz_t min;
  mpz_init(min);
  mpz_add_ui(min,n,1);

  mpz_t temp;
  mpz_init(temp);
  cout<<"GCD:"<<endl;
  for(size_t i=0;i<c;i++)
    {
      mpz_abs(temp,g[i]);
      mpz_init(gcd[i]);
      mpz_gcd(gcd[i],temp,n);
      cout<<"gcd["<<i<<"]:";
      gmp_printf("%Zd ",gcd[i]);
      if(mpz_cmp_ui(gcd[i],1)>0&&mpz_cmp(gcd[i],min)<0)
	{
	  mpz_set(min,gcd[i]);
	}
    }
  cout<<endl;
mpz_add_ui(temp,n,1);
if(mpz_cmp(min,temp)==0){
  gmp_printf("%Zd is a prime number \n",n);
  return 0;}
gmp_printf("The least prime factor of %Zd is %Zd \n",n,min);


  //delete arrays
  for(size_t i=0; i<=k;i++) 
    {
      pow2_i =pow(2,i);
      pow2_ki = c/pow2_i;
      for(size_t j=0;j<pow2_ki;j++)
	{
	  delete[] M1[i][j];
	  delete[] M2[i][j];
	}
      delete[] M1[i];
      delete[] M2[i];
    }
    delete[] M1;
    delete[] M2;
    return 0;
}

void buildup1 (size_t k, size_t r, mpz_t n, mpz_t ***M)
{
  for(size_t j=0; j<r;j++)
    {
       mpz_init_set_ui(M[0][j][0],j+1);
      mpz_init_set_ui(M[0][j][1],1);
    }

  for(size_t i=1; i<=k; i++)
    {
      for(size_t j=0; j<pow(2,k-i); j++)
	{
	  size_t pow2i1 = pow(2,i-1);
	  multi(M[i-1][2*j],M[i-1][2*j+1],M[i][j], pow2i1,n);
	}
    }
  return;
}




void buildup2 (mpz_t *u,size_t k, size_t r, mpz_t n, mpz_t ***M)
{
  for(size_t j=0; j<r;j++)
    {
      mpz_t temp;
      mpz_init_set(temp,u[j]);
      mpz_neg(temp,temp);
      mpz_mod(temp,temp,n);
       mpz_init_set(M[0][j][0],temp);
      mpz_init_set_ui(M[0][j][1],1);
    }

  for(size_t i=1; i<=k; i++)
    {
      for(size_t j=0; j<pow(2,k-i); j++)
	{
	  size_t pow2i1 = pow(2,i-1);
	  multi(M[i-1][2*j],M[i-1][2*j+1],M[i][j], pow2i1,n);
	}
    }
  return;
}

void multi (mpz_t *a, mpz_t * b, mpz_t * ans, size_t deg, mpz_t n) //naive 
{
  mpz_t A,B,ANS;
  mpz_init(A);
  mpz_init(B);
  mpz_init(ANS);
  size_t temp = 0;
  for(size_t k=0; k<=2*deg; k++)
    {
      mpz_init_set_ui(ans[k],0);
      size_t j=0;
      if(k>deg){j=k-deg;}
      for(j;j<=min(deg,k);j++)
	{
	  mpz_init_set(ANS,ans[k]);
	  mpz_init_set(A,a[j]);
	  mpz_init_set(B,b[k-j]);
	  mpz_mul(ans[k],A,B);
	  mpz_add(ans[k],ans[k],ANS);
	  mpz_mod(ans[k],ans[k],n);
	}
	}
  return;
}
  

void rem(mpz_t *a,mpz_t *b,mpz_t *&rem,size_t dA,size_t dB,size_t &dr, mpz_t n)
{

  size_t dR=dA;
  mpz_t* r;
  r = new mpz_t [dR];
  mpz_t q;
  for(size_t i = 0; i<=dR; i++)
    {
      mpz_init_set(r[i],a[i]);
    }

   size_t i=dA-dB;
   int stay=1;
   while(i>=0 &&stay!=0)
    {
      if (dR==dB+i && mpz_cmp_ui(r[dR],0)!=0)
	{

	  mpz_init_set(q,r[dR]);
	  int stay2=1;
	  size_t j=dR;
	  while(j>=i && stay2!=0)
	      {
		mpz_t temp;
		mpz_init(temp);
		//r[j]=(r[j]-q*b[j-i])%n;
		mpz_mul(temp,q,b[j-i]);
		mpz_sub(temp,r[j],temp);
		mpz_mod(temp,temp,n);
		mpz_set(r[j],temp);
		if(j>0){j--;}else{stay2=0;}
	      }
	    mpz_set_ui(r[dR],0);
	    dR=dR-1;
	}
      else
	{
	  dR=dR-1;
	}
      if(i>0){ i--;}else{stay=0;}
    }
    dr=dR;
  rem = new mpz_t [dr];
  for(int i =dr; i>=0;i--)
    {
      mpz_init_set(rem[i],r[i]);
      mpz_mod(rem[i],rem[i],n);
     
    }
  return;
}


void fastMultiEval(mpz_t *f,size_t k, size_t r,mpz_t * fu,mpz_t ***M, mpz_t n)
{
  mpz_t *u;
  u = new mpz_t [r];
  mpz_t temp;
  cout<<"u:"<<endl;
  for(int i=0;i<r;i++)
    {
      mpz_init_set_ui(temp,i*r);
      mpz_mod(temp,temp,n);
      mpz_init_set(u[i],temp);
      if(mpz_cmp_ui(u[i],0)<0){mpz_add(u[i],u[i],n);}
    }
  //10.3
  buildup2(u,k,r,n,M);
  //10.5
  godown(f,u,M,r,fu,k,n,0);
  return;
}


  
void godown(mpz_t *f,mpz_t *u,mpz_t ***M,size_t deg,mpz_t *fu,size_t k,mpz_t n,size_t start)
  {
    mpz_t temp;
    mpz_init(temp);
  if(k==0)
    {
      mpz_mod(temp, f[0],n);
      mpz_init_set(fu[start],temp);
      return;
    }
  mpz_t *r0;
   mpz_t * r1;
  size_t dR0=floor(deg/2),dR1=floor(deg/2);
  r0 = new mpz_t [dR0];
  r1 = new mpz_t [dR1];
  rem(f,M[k-1][2*start],r0,deg,pow(2,k-1),dR0,n);
  rem(f,M[k-1][2*start+1],r1,deg,pow(2,k-1),dR1,n);
  godown(r0,u,M,dR0,fu,k-1,n,2*start);
  godown(r1,u,M,dR1,fu,k-1,n,2*start+1);
   return;
}
