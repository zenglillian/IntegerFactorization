#include <iostream>
#include <cmath>
#include <numeric>
#include <complex>
using namespace std;

void buildup1(int k,int r,int n,int ***M);//10.3 step 1: building up subproduct tree
void fastMultiEval(int *f,int k,int r,int * fu,int ***M, int n);//10.7 step 2: fast multipoint evaluation
void godown(int *f,int *u,int ***M,int deg,int *fu,int k,int n,int start);//10.5 going down subproduct tree
void rem(int*a,int*b,int*&rem,int dA,int dB,int &dr, int n); //remainder of two polynomials a/b
void buildup2(int *u,int k,int r,int n,int ***M); // 10.3 building up subproduct tree for step 2
void multi(int *a,int *b,int *ans,int size,int n);//naive way for polynomial multiplication
int GCD( int a, int b,int*x, int *y); //extended Euclidean algorithm for gcd
void fastMulti(int n,int *a, int *b,int deg,int *c); //fast polynomial multiplication
void FFT(int num,const int length,complex<double> *f,complex <double>w,int deg, complex<double> *ans,int index,int multiple);//Fast fourier transformation
void IFFT(int num,const int length, complex<double> *alpha,complex<double> w, complex<double> *ans); //inverser FFT

const double pi = 4*atan(1);
int main()
{
  int n,b,c,k,pow2_ki,pow2_i;
  //set n
  n=31*29;
  //upper bound:b
  b = ceil(sqrt(n));
  //choose c to be the power of 2 that is greater than sqrt(b). And 2^k= c with integer k
  
  k=ceil(log2(ceil(sqrt(b))));
  c = pow(2,k);
  cout<<"n :"<<n<<" b: "<<b<<" c: "<<c<<" k: "<<k<<endl;

  //3D array
  int ***M1;
  int ***M2;
  M1 = new int **[k+1]; //for step 1
  M2 = new int **[k+1]; //for step 2
  
  for(int i=0; i<=k;i++) 
    {
      pow2_i =pow(2,i);
      pow2_ki = c/pow2_i;
      M1[i]=new int*[pow2_ki];
      M2[i]=new int*[pow2_ki];
      for(int j=0;j<pow2_ki;j++)
	{M1[i][j] = new int[pow2_i+1];M2[i][j] = new int[pow2_i+1]; }
    }

  
  //step 1: 10.3
  buildup1(k,c,n,M1); //coefficients are M[k][0]
  cout<<"step1: buildup1 result: f"<<endl;
  for(int i = 0; i<=c;i++)
    { cout<<"f["<<i<<"]:"<<M1[k][0][i]<<" ";}
  cout<<endl;


  //step 2: 10.7
  int * g;
  g = new int [c];
  
  fastMultiEval(M1[k][0],k,c,g,M2,n);
  
  cout<<"step2: fastMultiEval result: gi" <<endl;
  for(int i = 0; i<c; i++)
    {
      cout<<"g["<<i<<"]:"<<g[i]<<" ";
    }
  cout<<endl;

  //final step: find the least prime factor of n
  int *gcd;
  int x,y;
  gcd = new int [c];
  gcd[0] = GCD(g[0],n,&x,&y);
  int min=n+1;

 
  cout<<"GCD:"<<endl;
  for(int i=0;i<c;i++)
    {
      
      gcd[i]=GCD(abs(g[i]),n,&x,&y);
      cout<<"gcd["<<i<<"]:"<<gcd[i]<<" ";
      if((gcd[i]>1)&&(gcd[i]<min))
	{
	  min=gcd[i];
	}
    }
  cout<<endl;
  if(min==n+1){cout<<n<<" is a prime number"<<endl;return 0;}
  cout<<"The least prime factor of "<<n<<" is: "<<min<<endl;




  //delete arrays
    for(int i=0; i<=k;i++) 
    {
      pow2_i =pow(2,i);
      pow2_ki = c/pow2_i;
      for(int j=0;j<pow2_ki;j++)
	{delete[] M1[i][j];delete[] M2[i][j]; }
      delete[] M1[i];
      delete[] M2[i];
    }
    delete[] M1;
    delete[] M2;
  return 0;
  
}





//Building up subproduct tree for step 1
void buildup1(int k,int r,int n,int ***M)
{
for(int j = 0; j<r;j++)
  {

    //put values in M0
    M[0][j][0]=j+1;
    M[0][j][1]=1;
  }

  //put values in Mij
 for(int i=1;i<=k;i++)
   {
     for(int j=0;j<pow(2,k-i);j++)
       {
	 int pow2i1 = pow(2,i-1);
	 fastMulti(n,M[i-1][2*j],M[i-1][2*j+1],pow2i1,M[i][j]);
	 //  multi(M[i-1][2*j],M[i-1][2*j+1],M[i][j],pow2i1,n);//multi
       
       }
   }			        
 return;
}





void fastMultiEval(int *f,int k,int r,int *fu,int ***M, int n)
{
  int *u;
  u = new int [r];
  cout<<"u:"<<endl;
  for(int i=0;i<r;i++)
    {
      u[i] = (i*r)% n;
      if(u[i]<0){u[i]=u[i]+n;}
    }
  //10.3
  buildup2(u,k,r,n,M);
  //10.5
  godown(f,u,M,r,fu,k,n,0);
  return;
}


//Building up subproduct tree for step 2
void buildup2(int *u,int k,int r,int n,int ***M)
{
for(int j = 0; j<r;j++)
  {

    //put values in M0
    M[0][j][0]=(-u[j])%n;
    if(M[0][j][0]<0){M[0][j][0]=M[0][j][0]+n;}
    M[0][j][1]=1;
  }

  //put values in Mij
 for(int i=1;i<=k;i++)
   {
     for(int j=0;j<pow(2,k-i);j++)
       {
	 int pow2i1 = pow(2,i-1);
	  fastMulti(n,M[i-1][2*j],M[i-1][2*j+1],pow2i1,M[i][j]);
	 //	  multi(M[i-1][2*j],M[i-1][2*j+1],M[i][j],pow2i1,n);//multi

       }
     
   }			        
 return;
}


void godown(int *f,int *u,int ***M,int deg,int *fu,int k,int n,int start) 
{

  if(k==0)
    {
      fu[start]=f[0]%n;
      return;
    }
  int *r0;
  int dR0=floor(deg/2),dR1=floor(deg/2);
  int* r1;
  r0 = new int [dR0];
  r1 = new int [dR1];
  rem(f,M[k-1][2*start],r0,deg,pow(2,k-1),dR0,n);
  rem(f,M[k-1][2*start+1],r1,deg,pow(2,k-1),dR1,n);
  godown(r0,u,M,dR0,fu,k-1,n,2*start);
  godown(r1,u,M,dR1,fu,k-1,n,2*start+1);
   return;
}







void  multi(int *a,int *b,int *ans,int size,int n) //naive way of polynomial multiplication
{
  int A,B,ANS;
  for(int k=0;k<=2*size;k++)
    {
      ans[k]=0;
      for(int i=max(0,k-size);i<=min(size,k);i++)
	{
	  ANS=ans[k];
	  A=a[i];
	  B=b[k-i];
	  ans[k] = (ANS+A*B)%n;
	}
    }
  
  return;
}



int GCD( int a, int b,int*x, int *y)
{
  if (a==0)
    {
      *x=0;
      *y=1;
      return b;
    }
  int x1,y1;
  int gcd = GCD (b%a,a,&x1,&y1);

  *x=y1-(b-a)*x1;
  *y=x1;
  return gcd;
}
    


void rem(int*a,int*b,int*&rem,int dA,int dB,int &dr, int n)
{
  int dR=dA;
  int* r;
  r = new int [dR];
  int q;
  for(int i = 0; i<=dR; i++) {r[i] = a[i];}
  for(int i=dA-dB; i>=0; i--)
    {
      if (dR==dB+i && r[dR]!=0)	{
	  q = r[dR];
	    for(int j=dR; j>=i; j--)
	      {
		r[j]=(r[j]-q*b[j-i])%n;
		if(r[j]<0){r[j]=r[j]+n;}  }
	    r[dR]=0;
	    dR=dR-1;	}
      else{dR=dR-1;}
    }
    dr=dR;
  rem = new int[dr];
  for(int i =dr; i>=0;i--) {
      rem[i]=r[i]%n;
      if(rem[i]<0){rem[i]=rem[i]+n;}  }
  return;
} 






void fastMulti(int n,int *a, int *b,int deg,int *c)
{
   //step1: compute primitive dCth  root of unity w & w^2 ...w^(dC-1)
  if(deg<1){cout<<"error!"<<endl;return;}
  int dC = 2*(deg); //length of c -1 but this is the length of a and b to be compute in fft and ifft
   complex<double> w = polar(1.0,2*pi/dC);
  complex<double> *a1;
  complex<double>*b1;
  complex<double>*c1;
  a1= new complex<double> [dC];
  b1= new complex<double> [dC];
  c1= new complex<double> [dC];
  for(int i=0;i<dC;i++) //transfer from int to complex<double
    {
      if(i<=deg){a1[i]=a[i];b1[i]=b[i];}
      else{a1[i]=0;b1[i]=0;}
    }
  //step2: find the DFT of a and b
  complex<double> *A;
  complex<double> *B;
  complex<double> *C;
  A= new complex<double> [dC];
  B= new complex<double> [dC];
  C= new complex<double> [dC];
  FFT(dC,dC,a1,w,dC-1,A,0,1);
  FFT(dC,dC,b1,w,dC-1,B,0,1);
  //step3: pointwise multiplication
  for(int i=0;i<dC;i++)
    {
      C[i]=A[i]*B[i];
    }
  //step4: inverse of the result
  IFFT(dC,dC, C, w, c1);
  for(int i=0;i<dC;i++)
    {if(real(c1[i])<0){c[i]=int(real(c1[i])-0.5)%n+n;}
      else{c[i]=int(real(c1[i]+0.5))%n+n;}
      c[i]=c[i]%n;
    }
      
      //because the double is not represent exactly in binary, so add or subtract 0.5 to avoid mistake
  c[dC]=1;
  //the highest term in polynomial multiplication in pollard-strassen method of integer factorization is always 1

  c[0]=c[0]-1; 
 // coeff(first term)+coeff(last term) in the true result = coeff(first term) in this computation
 
  delete []a1;
  delete []b1;
  delete []c1;
  delete []A;
  delete []B;
  delete []C;
  return;
}
  






void FFT(int num,const int length,complex<double> *f,complex <double> w,int deg, complex<double> *ans,int index,int multiple)
{
  if(deg==0)
    {
      ans[index] = f[0];
      return;
    }
 
  int dR =num/2-1;
  int dR_1= dR+1;
  complex<double> *r0;
  r0 = new complex<double> [dR_1];
  complex <double> *r1;
  r1 = new complex<double> [dR_1];
  for(int j=0; j<=dR; j++)
    {
      r0[j] = f[j]+f[j+dR_1];
      r1[j] = pow(w,j)*(f[j]-f[j+dR_1]);
    }

  FFT(num/2,length,r0,pow(w,2),dR,ans,index,multiple*2);
  FFT(num/2,length,r1,pow(w,2),dR,ans,index+length/num,multiple*2);
    return;
}
      
      

  void IFFT(int num,const int length, complex<double> *alpha, complex<double> w,complex<double> *ans)
{
  for(int i=0;i<num;i++)
    {
      alpha[i] = conj(alpha[i]);
    }
  FFT(num,length,alpha,w,num-1,ans,0,1);
  for(int i=0;i<num;i++)
    {
      ans[i] = conj(ans[i])/double(num);
    }
  return;
}
