//Filename: PollardStrassen_multiprec.cpp
//Created by Lillian Zeng
//baby steps giant steps algorithm for integer factorization
/*This program should be link against gmp library and mpfr library, for example
 g++ PollardStrassen_multiprec.cpp -lgmp -lgmpxx -lmpfr */
#include <iostream>
#include <cmath>
#include <numeric>
#include <gmp.h>
#include <gmpxx.h>
#include <stdio.h>
#include<mpfr.h>
using namespace std;
mp_rnd_t rnd = GMP_RNDN;
void buildup1 (size_t k, size_t r, mpz_t n, mpz_t ***M); //10.3 step 1: building up subproduct tree
void buildup2 (mpz_t *u,size_t k, size_t r, mpz_t n, mpz_t ***M);// 10.3 building up subproduct tree for step 2
void multi (mpz_t *a, mpz_t * b, mpz_t * ans, size_t deg, mpz_t n);//naive way for polynomial multiplication
void FastMultiplication(mpz_t n,mpz_t*a,mpz_t* b,size_t deg, mpz_t *c);
void FFT(size_t num,const size_t length,mpfr_t **f, mpfr_t *w,size_t deg,mpfr_t **ans, size_t index,size_t multiple);
void IFFT(size_t num,const size_t length, mpfr_t **alpha, mpfr_t *w,mpfr_t **ans);
void rem(mpz_t *a,mpz_t *b,mpz_t *rem,size_t dA,size_t dB,size_t &dR, mpz_t n);//remainder of two polynomials a/b
void fastMultiEval(mpz_t *f,size_t k, size_t r,mpz_t * fu,mpz_t ***M, mpz_t n);//10.7 step 2: fast multipoint evaluation
void godown(mpz_t *f,mpz_t *u,mpz_t ***M,size_t deg,mpz_t *fu,size_t k,mpz_t n,size_t start); //10.5 going down subproduct tree

  int main()
{
  size_t precision=200; //may be set differently based on how big the integer to be factored is
  mpfr_set_default_prec(precision); 
  cout<<"precision:"<<precision<<endl;
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
    /*cout<<"buildup1  result: f"<<endl;
	 for(int i = 0; i<=c;i++)
    { cout<<"f["<<i<<"]:";
      gmp_printf("%Zd ",M1[k][0][i]);}
  cout<<endl;
	*/

  //step2: 10.7
  mpz_t *g;
  g = new mpz_t [c];
  fastMultiEval(M1[k][0],k,c,g,M2,n);

  /*     cout<<"fastMultiEval result: g"<<endl;
  for(size_t i = 0; i<c;i++)
    { cout<<"g["<<i<<"]:";
      gmp_printf("%Zd ",g[i]);}
  cout<<endl;
  */

  for(size_t i=0;i<c;i++)
    {if(mpz_cmp_ui(g[i],0)==0)//use trial division when the two prime factors are very close to each other
	{
	  
	  for(size_t j=i*c;j<=c*(i+1);j++)
	    {
	      if(mpz_divisible_ui_p(n,j)!=0)
		{
		  cout<<"The least prime factor is "<<j<<endl;return 0;
		}
	    }
	}
    }
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
  //cout<<"GCD:"<<endl;
  for(size_t i=0;i<c;i++)
    {
      
      mpz_abs(temp,g[i]);
      mpz_init(gcd[i]);
      mpz_gcd(gcd[i],temp,n);
      //cout<<"gcd["<<i<<"]:";
      // gmp_printf("%Zd ",gcd[i]);
      if(mpz_cmp_ui(gcd[i],1)>0&&mpz_cmp(gcd[i],min)<0)
	{mpz_set(min,gcd[i]);	}
      
    }
 
mpz_add_ui(temp,n,1);
if(mpz_cmp(min,temp)==0){
  gmp_printf("%Zd is a prime number \n",n);
  return 0;}
gmp_printf("The least prime factor of %Zd is %Zd \n",n,min);


 
  //clear
 mpz_clear(n);
 mpz_clear(temp);
 mpz_clear(min);
 for(size_t i=0;i<c;i++)
   {
     mpz_clear(g[i]);
     mpz_clear(gcd[i]);
   }
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
	  //  multi(M[i-1][2*j],M[i-1][2*j+1],M[i][j], pow2i1,n);
	  FastMultiplication(n,M[i-1][2*j], M[i-1][2*j+1],pow2i1,M[i][j]);
	   
	}
    }
  return;
}




void buildup2 (mpz_t *u,size_t k, size_t r, mpz_t n, mpz_t ***M)
{
  mpz_t temp;
  mpz_init(temp);
  for(size_t j=0; j<r;j++)
    {
      
      mpz_set(temp,u[j]);
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
	  // multi(M[i-1][2*j],M[i-1][2*j+1],M[i][j], pow2i1,n);
	  
	     FastMultiplication(n,M[i-1][2*j], M[i-1][2*j+1],pow2i1,M[i][j]);
	}
    }
  mpz_clear(temp);
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
  mpz_clear(A);mpz_clear(B);mpz_clear(ANS);
  return;
}



void FastMultiplication(mpz_t n,mpz_t*a,mpz_t* b, size_t deg, mpz_t *c)
{
  mpfr_t PI;
mpfr_init(PI);
mpfr_const_pi(PI,rnd);
  mpfr_t temp1;
  mpfr_init(temp1);
  // step1: compute primitive dCth root of unity, w
  size_t dC=2*deg; //length of c-1
   mpfr_t theta;
  mpfr_init(theta);
    mpfr_t *w;
  w = new mpfr_t [2];
  mpfr_init(w[0]);
  mpfr_init(w[1]);
  mpfr_mul_ui(theta,PI,2,rnd);
  mpfr_div_ui(theta,theta,dC,rnd);
  mpfr_sin_cos(w[1],w[0],theta,rnd);
 
  mpfr_t **a1;
  mpfr_t **b1;
  mpfr_t **c1;
  a1 = new mpfr_t *[dC+1];
  b1 = new mpfr_t *[dC+1];
  c1 = new mpfr_t *[dC+1];
 
  for(size_t i=0;i<=dC;i++)
    {  a1[i] = new mpfr_t [2];b1[i] = new mpfr_t [2];c1[i]=new mpfr_t [2];
      mpfr_init(a1[i][0]); mpfr_init(a1[i][1]); mpfr_init(b1[i][0]); mpfr_init(b1[i][1]);
      mpfr_init(c1[i][0]); mpfr_init(c1[i][1]);
      mpfr_set_ui(a1[i][0],0,rnd);
      mpfr_set_ui(b1[i][0],0,rnd);
      mpfr_set_ui(a1[i][1],0,rnd);
      mpfr_set_ui(b1[i][1],0,rnd);
      if(i<=deg)
	{
	  mpfr_set_z(a1[i][0],a[i],rnd);
	  mpfr_set_z(b1[i][0],b[i],rnd);
	}

    }
 
  //step2: find the DFT of a and b
  mpfr_t **A;
  mpfr_t **B;
  mpfr_t **C;
  A =new mpfr_t *[dC+1];
  B =new mpfr_t *[dC+1];
  C =new mpfr_t *[dC+1];

  for(size_t i=0; i<=dC;i++)
    {
     
      A[i] = new mpfr_t [2];
      B[i] = new mpfr_t [2];
      C[i] = new mpfr_t [2];
      mpfr_init(A[i][0]);mpfr_init(A[i][1]);
      mpfr_init(B[i][0]);mpfr_init(B[i][1]);
      mpfr_init(C[i][0]);mpfr_init(C[i][1]);
    }
  FFT(dC,dC,a1,w,dC,A,0,1);
  FFT(dC,dC,b1,w,dC,B,0,1);
  //step 3: pointwise multiplication
  for(size_t i=0; i<=dC;i++)
    {
      //C=A*B
      mpfr_mul(C[i][0],A[i][0],B[i][0],rnd);
      mpfr_mul(temp1,A[i][1],B[i][1],rnd);
      mpfr_sub(C[i][0],C[i][0],temp1,rnd);
      
      mpfr_mul(C[i][1],A[i][0],B[i][1],rnd);
      mpfr_mul(temp1,A[i][1],B[i][0],rnd);
      mpfr_add(C[i][1],C[i][1],temp1,rnd);
     
    }
   //step4: inverse of the result
   IFFT(dC, dC, C,w,c1);

   mpfr_set_ui(temp1,1,rnd);
   mpfr_div_ui(temp1,temp1,2,rnd);
  
   for(size_t i=0;i<dC;i++)
     {
       mpz_init(c[i]);
       if(mpfr_cmp_ui(c1[i][0],0)<0)
	 {
	   mpfr_sub(c1[i][0],c1[i][0],temp1,rnd);
	   mpfr_rint_trunc(c1[i][0],c1[i][0],rnd);
	   mpfr_get_z(c[i],c1[i][0],rnd);
	   mpz_mod(c[i],c[i],n);}
       else{
	 mpfr_add(c1[i][0],c1[i][0],temp1,rnd);
	 mpfr_rint_trunc(c1[i][0],c1[i][0],rnd);
	 mpfr_get_z(c[i],c1[i][0],rnd);	 
	 mpz_mod(c[i],c[i],n);}
      
     }
   mpz_init_set_ui(c[dC],1);    
   mpz_sub_ui(c[0],c[0],1);
   mpfr_clear(PI); mpfr_clear(temp1);mpfr_clear(theta);mpfr_clear(w[0]);mpfr_clear(w[1]);
   for(size_t i=0;i<dC;i++)
     {
        mpfr_clear(a1[i][0]); mpfr_clear(b1[i][0]); mpfr_clear(c1[i][0]);
       mpfr_clear(a1[i][1]); mpfr_clear(b1[i][1]); mpfr_clear(c1[i][1]);
       mpfr_clear(A[i][0]); mpfr_clear(B[i][0]); mpfr_clear(C[i][0]);
       mpfr_clear(A[i][1]); mpfr_clear(B[i][1]); mpfr_clear(C[i][1]);
        delete[] a1[i]; delete[] b1[i]; delete[] c1[i];
         delete[] A[i]; delete[] B[i];delete[] C[i];
     }
   
     delete[] w;
    delete[] a1; delete[] b1;delete[] c1;
   delete[] A; delete[] B; delete[] C;
    
   return;
}

   void FFT(size_t num,const size_t length,mpfr_t **f, mpfr_t *w,size_t deg,mpfr_t **ans, size_t index,size_t multiple)
{
 
  mpfr_t PI;
mpfr_init(PI);
mpfr_const_pi(PI,rnd);
  if(deg==0)
    {
      mpfr_set(ans[index][0],f[0][0],rnd);
      mpfr_set(ans[index][1],f[0][1],rnd);
     
       return;
    }
  size_t dR=num/2-1;
  size_t dR_1 = dR+1;
  mpfr_t **rr0;
  mpfr_t **rr1;
  mpfr_t theta;
  mpfr_t thetaj;
  mpfr_t wj0;
  mpfr_t wj1;
  mpfr_t *w2; //w^2
  mpfr_t temp0;
  mpfr_t temp1;
  mpfr_t temp2;

  mpfr_init(theta);
  mpfr_mul_ui(theta,PI,2,rnd);
  mpfr_div_ui(theta,theta,num,rnd);
  mpfr_init(thetaj);
  mpfr_init(wj0);
  mpfr_init(wj1);
  w2 = new mpfr_t [2];
  mpfr_init(w2[0]);
  mpfr_init(w2[1]);
  mpfr_init(temp0);
  mpfr_init(temp1);
  mpfr_init(temp2);
  rr0 = new mpfr_t *[dR_1];
  rr1 = new mpfr_t *[dR_1];
  for (size_t j=0;j<=dR;j++)
    {
      rr0[j] = new mpfr_t [2];
      rr1[j] = new mpfr_t [2];
      mpfr_init(rr0[j][0]);
      mpfr_add(rr0[j][0],f[j][0],f[j+dR_1][0],rnd);
      mpfr_init(rr0[j][1]);
      mpfr_add(rr0[j][1],f[j][1],f[j+dR_1][1],rnd);
      
      mpfr_mul_ui(thetaj,theta,j,rnd);
      mpfr_init(rr1[j][0]); mpfr_init(rr1[j][1]);
      mpfr_sub(temp0,f[j][0],f[j+dR_1][0],rnd);
      mpfr_sub(temp1,f[j][1],f[j+dR_1][1],rnd);
      mpfr_sin_cos(wj1,wj0,thetaj,rnd);
      mpfr_mul(rr1[j][0],temp0,wj0,rnd);
      mpfr_mul(temp2,temp1,wj1,rnd);
      mpfr_sub(rr1[j][0],rr1[j][0],temp2,rnd);
      mpfr_mul(rr1[j][1],temp1,wj0,rnd);
      mpfr_mul(temp2,temp0,wj1,rnd);
      mpfr_add(rr1[j][1],rr1[j][1],temp2,rnd);
    }

      mpfr_mul_ui(thetaj,theta,2,rnd);
       mpfr_sin_cos(w2[1],w2[0],thetaj,rnd);
       
  FFT(num/2,length,rr0,w2,dR,ans,index,multiple*2); 
  FFT(num/2,length,rr1,w2,dR,ans,index+length/num,multiple*2);
  
   for(size_t i=0;i<dR;i++)
    {
      mpfr_clear(rr0[i][0]);mpfr_clear(rr0[i][1]);mpfr_clear(rr1[i][0]);mpfr_clear(rr1[i][1]);
      delete[] rr0[i];
      delete[] rr1[i];
    }
  mpfr_clear(PI); mpfr_clear(theta); mpfr_clear(thetaj); mpfr_clear(wj0); mpfr_clear(wj1); mpfr_clear(temp0); mpfr_clear(temp1); mpfr_clear(temp2); mpfr_clear(w2[0]); mpfr_clear(w2[1]);
  delete[]  w2;
  delete[] rr0;
  delete[] rr1;
  
  return;
}



void IFFT(size_t num,const size_t length, mpfr_t **alpha, mpfr_t *w,mpfr_t **ans)
{
 
  for(size_t i=0;i<num;i++)
    {
      mpfr_neg(alpha[i][1],alpha[i][1],rnd); 
    }
  FFT(num,length, alpha,w,num-1,ans,0,1);
  for(size_t i=0;i<num;i++)
    {
      mpfr_neg(ans[i][1],ans[i][1],rnd);
      mpfr_div_ui(ans[i][0],ans[i][0],num,rnd);
      mpfr_div_ui(ans[i][1],ans[i][1],num,rnd);
    }
  return;
}




void rem(mpz_t *a,mpz_t *b,mpz_t *rem,size_t dA,size_t dB,size_t &dr, mpz_t n)
{

 
  size_t dR=dA;
  mpz_t* r;
 
  r = new mpz_t [dR+1]; //
 
  mpz_t q;
  mpz_init(q);

  for(size_t i = 0; i<=dR; i++)
    {
      mpz_init_set(r[i],a[i]);
    }
 
   size_t i=dA-dB;
   int stay=1;
   mpz_t temp;
   mpz_init(temp);
 
   while(i>=0 &&stay!=0)
    {
      
      if (dR==dB+i && mpz_cmp_ui(r[dR],0)!=0)
	{
	 
	  mpz_set(q,r[dR]);
	  int stay2=1;
	  size_t j=dR;
	  while(j>=i && stay2!=0)
	      {
 	
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
  
   
 
  for(size_t i =0; i<=dr;i++)
    {
      mpz_init_set(rem[i],r[i]);
      mpz_mod(rem[i],rem[i],n);
      
    }
 
  mpz_clear(temp);
 
  return;
}


void fastMultiEval(mpz_t *f,size_t k, size_t r,mpz_t * fu,mpz_t ***M, mpz_t n)
{
  mpz_t *u;
  u = new mpz_t [r];
  mpz_t temp;
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

  for(size_t i=0;i<r;i++)
    {mpz_clear(u[i]);}
  delete[] u;
  mpz_clear(temp);
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
   size_t dR0=floor(deg/2);
    size_t dR1=floor(deg/2);

   
    mpz_t* r0= new mpz_t [dR0+1]; 
    mpz_t* r1=new mpz_t [dR1+1];
   
  rem(f,M[k-1][2*start],r0,deg,pow(2,k-1),dR0,n);
  
  rem(f,M[k-1][2*start+1],r1,deg,pow(2,k-1),dR1,n);
 
  godown(r0,u,M,dR0,fu,k-1,n,2*start);
  godown(r1,u,M,dR1,fu,k-1,n,2*start+1);

  
  mpz_clear(temp);

 
   delete[] r0;
  delete[] r1;
 
   return;
}
