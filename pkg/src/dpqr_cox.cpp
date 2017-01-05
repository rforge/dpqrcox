#include <Rcpp.h>
using namespace Rcpp;
// [[Rcpp::export]]

NumericVector dpqr_cox(NumericVector xa){
double mybase(NumericVector, NumericVector);
NumericVector rety(xa[1]);
mybase(xa, rety);
return rety;
}

//ȫ�ֱ���
double *pbv,*lmbd,*srt,*xv;//��������probvec����������lambdvec����������,�Ա�������
double *betav;//�޸ĵĲ�������lambdvec
double *weightv;//���㺯��ֵ�õ�Ȩ������
double *outv;//�������
int pn,xn;//pn�ǲ����������ȣ�xn���Ա�����������
int ord,lower_tail,log_p,log_ord;//���ú�����,lower.tail,log.p,log
double tolerance=0.001;//����
int suggested_precision=20;//ȫ�ֳ��þ���
int hpnflag=1;//ʹ�ø߾��������ȱʡ���á�
//ȫ�ֱ�������


#include<math.h>
#include<assert.h>
#include<stdio.h>
//#include<alloc.h>
#include<stdlib.h>
//#include<limits.h>

//!�±���뿪ʼ
/*
#hc ������һ����ֵ����inx
#hc inx[0] Ϊ���ú����� 1 ��Ӧ dcox��2 ��Ӧ pcox��3 ��Ӧ qcox
#hc inx[1] Ϊ�Ա�����������
#hc inx[2] Ϊ������������
#hc inx[3] Ϊlower.tail��0 ��Ӧ FALSE��1 ��Ӧ TRUE
#hc inx[4] Ϊlog.p��0 ��Ӧ FALSE��1 ��Ӧ TRUE
#hc inx[5] Ϊlog��0 ��Ӧ FALSE��1 ��Ӧ TRUE
#hc inx[6]-inx[5+inx[1]] Ϊ�Ա�������
#hc inx[6+inx[1]]-inx[5+inx[1]+inx[2]] Ϊ��������probvec
#hc inx[6+inx[1]+inx[2]]-inx[5+inx[1]+2*inx[2]] Ϊ��������lambdvec
#hc inx[6+inx[1]+2*inx[2]]-inx[5+inx[1]+3*inx[2]] Ϊ��������lambdvec�Ĵ���˵��
#hc �����������lambdvecΪc(0.5,0.2,0.3)��
#hc ����˵��Ϊ 2.0 3.0   1.0
#hc ����С��һ��Ϊ�ڶ�����0.2������С��Ϊ��������0.3��������Ϊ��һ����0.5��
#hc ע��inx����Ϊ6+inx[1]+3*inx[2]
*/

  void basic_coef_mod(double * vec, int * ord, double * res, double bt_er, int n){
    int i;
    double bdt, bdts, bdms, x;
    res[ord[0]] = vec[ord[0]];
    for(i=1;i<n;i++){
      bdt = vec[ord[i]];
      bdts = vec[ord[i-1]];
      bdms = res[ord[i-1]];
      if(bdt<=1.0) x = bdt; else x =1.0;
      x /= (bdt-bdts+1.0);
      x *= bt_er;
      res[ord[i]] = x+bdt-bdts+bdms;
    }
  }

  void coef_mod(double * vec, int * ord, double * res, double error, int n){
    int i;
    double d_er=0.0, p_er=0.0, bt_er=1.0, err;
    basic_coef_mod(vec, ord, res, bt_er, n);
    if(n==1)return;
    for(i=0;i<n;i++){
      err = fabs(res[i]-vec[i]);
      p_er += err/vec[i];
      d_er += err;
    }
    if(p_er<d_er) bt_er = error/(1.5*d_er); else bt_er = error/(0.5*p_er);
    basic_coef_mod(vec, ord, res, bt_er, n);
  }



//���������Ŀռ�
double * create_vec(int n){
 double *v;
 v = (double *)malloc((n+1)*sizeof(double));
 if(!v){
   Rprintf("Error: cannot allocate memory of %d*sizeof(double)\n",n+1);
   throw(0);
 }
 return v;
}

//�ͷ�����
void free_vec(double *v){
  free(v);
}

//�ͷ�ϵ������
void free_c_matrix(double **cm, int n){
 int i;
 for(i=1; i<=n; i++)free(cm[i]);
 free(cm);
 return; 
}

//����ϵ������Ŀռ䲢����ϵ�������ֵcm[i][j]��X_1+...X_i��X_j�ķֲ�������ϵ��
//Ȼ�����Ȩ��weightv
void weightv_ini(void){
 int i,j,n=pn;
 double **cm;
 cm = (double **)malloc((n+1)*sizeof(double *));
 for(i=1; i<=n; i++)cm[i] = (double *)malloc((n+1)*sizeof(double));
 for(i=1; i<=n; i++){
  cm[i][i]=1.0;
  for(j=1;j<i;j++)cm[i][i] *= betav[j]/(betav[j]-betav[i]);
  if(i>=2){
    for(j=1;j<i;j++)cm[i][j] = cm[i-1][j]*betav[i]/(betav[i]-betav[j]);
  }
 }
 for(j=1;j<=n;j++){
   weightv[j]=0.0;
   for(i=j;i<=n;i++){
     weightv[j] += pbv[i]*cm[i][j];
   }
 }
 free(cm);
 return;
}


void coef_mod_from_1(double error, int n){
    int ord[n], i;
    double *vec, *res;
    double d_er=0.0, p_er=0.0, bt_er=1.0, err;//

    for(i=1; i<=n; i++)ord[i-1]=(int)srt[i]-1;
    vec=&(lmbd[1]);
    res=&(betav[1]);

    basic_coef_mod(vec, ord, res, bt_er, n);
    if(n==1)return;
    for(i=0;i<n;i++){
      err = fabs(res[i]-vec[i]);
      p_er += err/vec[i];
      d_er += err;
    }
    if(p_er<d_er) bt_er = error/(1.5*d_er); else bt_er = error/(0.5*p_er);
    basic_coef_mod(vec, ord, res, bt_er, n);
}

double prec_low_bound(void){
double pb=0,db=0,x,bb;//���ʺ������ܶȺ���������Ϊ��֤tolerance������Ҫ����С����
int i;
bb=fabs(log10(tolerance));
for(i=1;i<=pn;i++){
  x = log10(fabs(weightv[i])); 
  if(x>=pb)pb=x;
  x += log10(betav[i]);
  if(x>=db)db=x;
}
if(pb>db)return pb+bb;
else return db+bb;
}

void compute_dcox(double * xv, double * outv, int xn){
 int i,j;
 double y;
 if(hpnflag){
   Rprintf("Error: Precision is too small for dcox.\n");
   throw(0);
 }
 for(i=1;i<=xn;i++){
   y=0.0;
   for(j=1;j<=pn;j++){
     y += weightv[j]*betav[j]*expl(-betav[j]*xv[i]);
   }
   outv[i]=y;
 }
 return;
}

void compute_pcox(double * xv, double * outv, int xn){
 int i,j;
 double y;
 if(hpnflag){
   Rprintf("Error: Precision is too small for pcox.\n");
   throw(0);
 }
 for(i=1;i<=xn;i++){
   y=0.0;
   for(j=1;j<=pn;j++){
  // Rprintf("weightv[j]=%f.\n",weightv[j]);//!
  // Rprintf("betav[j]=%f.\n",betav[j]);//!
  // Rprintf("xv[i]=%f.\n",xv[i]);//!
     y += weightv[j]*(1.0-expl(-betav[j]*xv[i]));
   }
  // Rprintf("y=%f.\n",y);//!
   outv[i]=y;
 }
 return;
}

double mybase(NumericVector xa, NumericVector rety){//��ʼ��
double mymain(double);//�ѱ������������ʱ������ɾ����һ��
double iny;
int i;
pn=xa[2];
xn=xa[1];
//Rprintf("pn=%d.\n",pn);//!
//Rprintf("xn=%d.\n",xn);//!
pbv = create_vec(pn);
betav = create_vec(pn);
lmbd = create_vec(pn);
srt = create_vec(pn);
weightv = create_vec(xn);
xv = create_vec(xn);
outv = create_vec(xn);
for(i=1;i<=pn;i++){
  pbv[i]=xa[5+xn+i];
//  Rprintf("   pbv[%d]=%f\n",i,pbv[i]);
  lmbd[i]=xa[5+xn+pn+i];
//  Rprintf("   lmbd[%d]=%f\n",i,lmbd[i]);
  srt[i]=xa[5+xn+2*pn+i];
//  Rprintf("   srt[%d]=%f\n",i,srt[i]);
}
for(i=1;i<=xn;i++){
  xv[i]=xa[5+i];
//  Rprintf("   xv[%d]=%f\n",i,xv[i]);
}
ord = (int)xa[0];
lower_tail = (int)xa[3];
log_p = (int)xa[4];
log_ord = (int)xa[5];
//Rprintf("   ord=%d, lower_tail=%d, log.p=%d, log=%d\n",ord,lower_tail,log_p,log_ord);

coef_mod_from_1(tolerance, pn);//coef_mod_from_1(double error, int n)error��������������tolerance
//betav�趨����

for(i=1;i<=pn;i++){
//  Rprintf("   betav[%d]=%f\n",i,betav[i]);
}

weightv_ini();//����Ȩ��weightv

for(i=1;i<=pn;i++){
//  Rprintf("   weightv[%d]=%50.24f\n",i,weightv[i]);
}


if(prec_low_bound()>suggested_precision*log10(256)){
        Rprintf("Error: Precision is too small.\n");
        throw(0);
}else if(prec_low_bound()<=17)hpnflag=0;

if(ord==1)compute_dcox(xv,outv,xn);
else if(ord==2)compute_pcox(xv,outv,xn);
else{
 Rprintf("Error: wrong order.\n");
 throw(0);
}

iny=mymain(1+xv[1]);//�ģ������ǽ��������ĳ���
//Rprintf("   y=%50.24f\n",iny);//��

for(i=1;i<=xa[1];i++){
  rety[i-1] = outv[i];
}

free(outv);
free(xv);
free(weightv);
free(srt);
free(lmbd);
free(betav);
free(pbv);

return 1.0;
}

//!�±�������

//================================================
//innerhpn.h ȥ����extern int d_inner_precision��extern���η�
//================================================

#ifndef INCLUDE_DLOAT
#define INCLUDE_DLOAT 1
#ifndef Precision
#define Precision d_inner_precision 
int d_inner_precision;  //���Ǽ�¼�ĵ�ǰȫ�־��ȣ�
#endif
struct dloat{
int sign;
long int exponent;
unsigned char *mantissa;
};
struct dloat *create_dloat(void);
void free_dloat(struct dloat *v);
#define dloatis0(v) (!(v)->mantissa[1])
void letitbe0(struct dloat *v);
#define MAX_LMT (LONG_MAX/2-1) 
#define MIN_LMT (LONG_MIN/2+1) 
#define LetItGreat(u)    (u->exponent = MAX_LMT,    u->sign = 1, u->mantissa[1]=(unsigned char)1 )
#define LetItTiny(u)    (u->exponent = MIN_LMT,    u->sign = 1, u->mantissa[1]=(unsigned char)1 )
#define BASE 256.0
void float2d(float fv, struct dloat *v);
float dloat2f(struct dloat *v);
void dinv(struct dloat *u, struct dloat *v);
void dmul(int *is, struct dloat *w, struct dloat *u, struct dloat *v);
int dabsbig(struct dloat *u, struct dloat *v);
void dcopy(struct dloat *a, struct dloat *b);
void dadd(int *is, struct dloat *w, struct dloat *u, struct dloat *v);
void dexp(int *is, struct dloat *w, struct dloat *u);
#define dneg(u) if(!dloatis0(u)) (u)->sign*=-1
struct dloat **create_dloat_vector(int n);
void free_dloat_vector(struct dloat **v, int n);
void sum(struct dloat *w, struct dloat **a, struct dloat **b, int n);
void int2d(int i, struct dloat *v);
#endif

//================================================
//hpn.hɾ��ͷ�ļ��Ͷ�dloat���ⲿ˵��
//================================================

struct hpn{
  int precision;
  struct dloat *value;
};

/*����һ������Ϊ precision �ĸ߾�����*/
struct hpn *create_hpn(int precision);

/*�ͷŸ߾����� p*/
void free_hpn(struct hpn *p);

/*�Ѹ߾����� p ��ֵд������i, i��Ϊ0*/
void int2hpn(int i, struct hpn *p);

/*��ʮ��������߾����� p�����mλ�����ʱp���ƻ�*/
void show_hpn(struct hpn *p, int m);

/*�Ը߾��ȸ����� u, exp(u) --> w
�����������־λ is ���� -1��
����������is ���� 0 */
void hpn_exp(int *is, struct hpn *w, struct hpn *u);


//================================================
//innerhpn.chɾ��ͷ�ļ���int d_inner_precision=24;
//================================================


// //=======================================
// //ԭ����mpopos�Ĵ�����
// //=======================================

#define LOBYTE(x) ((unsigned char)((x)& 0xff))
#define HIBYTE(x) ((unsigned char)((x)>> 8 & 0xff))

/*�޷����� u+v --> w 
ע�� u��vΪnλ�ģ�������u[1..n]��v[1..n]
�� w Ϊn+1λ�ģ�������w[1..n+1]
*/
void d_inner_add(unsigned char w[], unsigned char u[], unsigned char v[], int n)
{
    unsigned short i=0;
    int j;
    for (j=n;j;j--) {
        i=v[j]+u[j]+HIBYTE(i);
        w[j+1]=LOBYTE(i);
    }
    w[1]=HIBYTE(i);
}

/*�޷����� u-v --> w 
����Ϊ������־λ flag ���� -1�����Ϊ 0 ������flag ���� 0 
u,v,w Ϊ n λ��*/
void d_inner_sub(int *flag, unsigned char w[], unsigned char u[], unsigned char v[], int n)
{
    unsigned short i=256;
    int j;
    for (j=n;j;j--) {
        i=255-v[j]+u[j]+HIBYTE(i);
        w[j]=LOBYTE(i);
    }
    *flag=HIBYTE(i)-1;
}

/*ȡ�����*/
void ddfun1(unsigned char u[], int n)
{
    unsigned short i=256;
    int j;
    for (j=n;j;j--) {
        i=255-u[j]+HIBYTE(i);
        u[j]=LOBYTE(i);
    }
}

/*�� v[1..n] ���Ƶ� u[1..n]*/
void d_inner_mov(unsigned char u[], unsigned char v[], int n)
{
    int j;
    for (j=1;j<=n;j++) u[j]=v[j];
}

#undef LOBYTE
#undef HIBYTE

/*ʹ��FFT�˷�*/
static int ivar1,ivar2;
#define IMAX(a,b) (ivar1=(a),ivar2=(b),(ivar1) > (ivar2) ? (ivar1) : (ivar2))

double *ddfun2(int nl, int nh){
    double *v;
    v=(double*)malloc((unsigned) (nh-nl+1)*sizeof(double));
    assert(v);
    return v-nl;
}

void ddfun3(double *v, int nl){
    free((char*)(v+nl));
}

#define RX 256.0

#define SWAP(a,b) c=(a);(a)=(b);(b)=c

void ddfun4(double data[], unsigned long nn, int sign){
    double c,tempi,wtemp,r,pr,pi,wi,theta;
    unsigned long n,g,m,t,i,j=1;
    n=nn<<1;
    for (i=1;i<n;i+=2) { 
        if (i<j) {
            SWAP(data[j+1],data[i+1]);
            SWAP(data[j],data[i]);
        }
        m=n>>1;
        while (m>=2 && m<j) {
            j-=m;
            m>>=1;
        }
        j+=m;
    }
    g=2;
    while(g<n){
        wi=0.0;
        r=1.0;
        theta=sign*(6.28318530717959/g);
        wtemp=sin(0.5*theta);
        pr=-2.0*wtemp*wtemp;
        pi=sin(theta);
        t=g << 1;
        for(m=1;m<g;m+=2){
            for(i=m;i<=n;i+=t) {
                j=g+i;
                tempi=r*data[j+1]+wi*data[j];
                c=r*data[j]-wi*data[j+1];
                data[j]=data[i]-c;
                data[j+1]=data[i+1]-tempi;
                data[i]+=c;
                data[i+1]+=tempi;
            }
            r=(wtemp=r)*pr-wi*pi+r;
            wi=wi*pr+wtemp*pi+wi;
        }
        g=t;
    }
}

void ddfun5(double data[], unsigned long n, int sign){
    unsigned long n1,i,i1,i2,i3,i4;
    double c1=0.5,c2,h1i,h1r,h2i,h2r,r,wi,pr,pi,wtemp,theta;
    theta=3.141592653589793/(double) (n>>1);
    if (sign==1){
        c2=-0.5;
        ddfun4(data,n>>1,1);
    }else{
        c2=0.5;
        theta*=-1.0;
    }
    wi=sin(theta);
    pi=wi;
    wtemp=sin(0.5*theta);
    pr=-2.0*wtemp*wtemp;
    r=1.0+pr;
    n1=n+3;
    for(i=2;i<=(n>>2);i++){
        i4=1+(i3=n1-(i2=1+(i1=i+i-1)));
        h1i=(data[i2]-data[i4])*c1;
        h1r=(data[i3]+data[i1])*c1;
        h2i=(data[i1]-data[i3])*c2;
        h2r=-(data[i4]+data[i2])*c2;
        data[i1]=r*h2r+h1r-wi*h2i;
        data[i2]=r*h2i+h1i+wi*h2r;
        data[i3]=wi*h2i+h1r-r*h2r;
        data[i4]=wi*h2r-h1i+r*h2i;
        r=(wtemp=r)*pr-wi*pi+r;
        wi=wi*pr+wtemp*pi+wi;
    }
    if(sign==1) {
        data[1]=(h1r=data[1])+data[2];
        data[2]=h1r-data[2];
    }else{
        data[1]=c1*((h1r=data[1])+data[2]);
        data[2]=c1*(h1r-data[2]);
        ddfun4(data,n>>1,-1);
    }
}


/*��λ�˷�
�޷����� u[1..n] �� v[1..m] --> w[1..n+m].
ע��˷��Ľ�������� n+m λ�ģ�Ҳ������ n+m-1 λ����
�ں�������£�w[1] ��0
*/
/* Uses Fast Fourier Transform to multiply the unsigned radix 256 integers u[1..n] and v[1..m], yielding a product w[1..n+m].*/
void d_inner_mul(unsigned char w[], unsigned char u[], unsigned char v[], int n, int m)
{
    int j,m1,n1=1;
    double y,t,*a,*b;

    m1=IMAX(m,n);
    while (m1>n1)n1 <<= 1;
    n1 <<= 1;
    a=ddfun2(1,n1);
    b=ddfun2(1,n1);
    for(j=1;n>=j;j++)
        a[j]=(double)u[j];
    for(j=n+1;n1>=j;j++) a[j]=0.0;
    for(j=1;m>=j;j++)    b[j]=(double)v[j];
    for(j=m+1;n1>=j;j++) b[j]=0.0;
    ddfun5(a,n1,1);
    ddfun5(b,n1,1);
    b[2]*=a[2];
    b[1]*=a[1];
    for(j=3;n1>=j;j+=2){
        b[j]=(t=b[j])*a[j]-b[j+1]*a[j+1];
        b[j+1]=b[j+1]*a[j]+t*a[j+1];
    }
    ddfun5(b,n1,-1);
    y=0.0;
    for(j=n1;j;j--) {
        t=b[j]/(n1>>1)+0.5+y;
        y=(unsigned long)(t/RX);
        b[j]=t-y*RX;
    }
    assert(y<RX);
    w[1]=(unsigned char)y;
    for(j=2;n+m>=j;j++)
        w[j]=(unsigned char)b[j-1];
    ddfun3(b,1);
    ddfun3(a,1);
}
#undef RX

/*
v[1..m] �ĵ��� --> u[1..n]
v ��С������ v[1] ֮��ע�� v[1] ���� 
u ��С������ u[1] ֮��ע�� u[1] ���� 
*/

#define MF 4
#define BI (1.0/256)
void d_inner_inv(unsigned char u[], unsigned char v[], int n, int m)
{
    float f1,f2;
    int i,j,g,m1;
    unsigned char *r,*s; 
    if(m<n)g=n;
    else g=m;
    r=(unsigned char*)malloc((2*g+1)*sizeof(unsigned char)); /*���丨������*/
    s=(unsigned char*)malloc((g+1)*sizeof(unsigned char));
    if(m<MF)m1=m;
    else m1=MF;
    f2=(float)v[m1];  /* v�ĸ������-->f2 */
    for(j=m1-1;j;j--) {
        f2*=BI;
        f2+=v[j];
    }
    f1=1.0/f2; /* ������ֵ�趨 */
    for(j=1;j<=n;j++){
        i=(int)f1;
        u[j]=(unsigned char)i;
        f1=256.0*(f1-i);
    }
    for(;;){
        d_inner_mul(r,u,v,n,m);
        d_inner_mov(s,&r[1],n);
        ddfun1(s,n);
        s[1]-=254;
        d_inner_mul(r,s,u,n,n);
        d_inner_mov(u,&r[1],n);
        for(j=2;j<n;j++)if(s[j])break;
        if(j==n){
            free(s);
            free(r);
            return;
        }
    }
}
#undef MF
#undef BI


// //=======================================
// //ԭ����dloat�Ĵ�����
// //=======================================
struct dloat *create_dloat(void){
    struct dloat *v;
    v = (struct dloat *)malloc(sizeof(struct dloat));
    v->mantissa = (unsigned char*) malloc((Precision+1)*sizeof(unsigned char));
    return v; 
}

void free_dloat(struct dloat *v){
    free(v->mantissa);
    free(v);
}

/*�ѷǸ�����i��ֵ��v*/
void int2d(int i, struct dloat *v){
    int k,j,base=(int)BASE;
    if(i<=0){
        Rprintf("Error: Input i of function `void int2d(int i, struct dloat *v)' must be not less than 0, Now, i= %d\n",i);
        throw(0);
    }
    v->sign = 1;
    v->exponent = 0;
    k=i;
    for(j=Precision; j>=1; j--){
        v->mantissa[j]=(unsigned char)k%base;
        k -= (int)v->mantissa[j];
        k/=base;
    }
    if(k){
        Rprintf("Error: Precision is too small, cannot transform i= %d\n",i);
        throw(0);
    }
    for(j=1; j<=Precision; j++){
        if(v->mantissa[j]) break;
    }
    if(j>Precision||v->mantissa[j]==0){
        letitbe0(v);
    }else{
    d_inner_mov(v->mantissa, &(v->mantissa[j-1]), Precision-j+1);
    v->exponent = Precision-j;
  }
}

/*�߾��ȸ����� v �� 0������ 1�� v ���� 0������ 0
ע����� 0 ���� 1<β��<256�� ��β���Ǹ�������λ��С��
*/
#define dloatis0(v) (!(v)->mantissa[1])

/*���߾��ȸ����� v ��ֵΪ 0*/
void letitbe0(struct dloat *v){
    int i;
    v->exponent = 0;
    v->sign = 1;
    for(i=1;i<=Precision;i++)v->mantissa[i]=(unsigned char)0;
    return;
}

/*�߾��ȸ�����ָ����������Сֵ*/
#define MAX_LMT (LONG_MAX/2-1) 
#define MIN_LMT (LONG_MIN/2+1) 

/*���߾��ȸ����� u ��һ���ܴ��ֵ*/
#define LetItGreat(u)    (u->exponent = MAX_LMT,    u->sign = 1, u->mantissa[1]=(unsigned char)1 )
/*���߾��ȸ����� u ��һ���ӽ�0��ֵ*/
#define LetItTiny(u)    (u->exponent = MIN_LMT,    u->sign = 1, u->mantissa[1]=(unsigned char)1 )


/*��ʾ�߾��ȸ�����
*/
/*#define showd(v) show_a_dloat((v), #v)
void show_a_dloat(struct dloat *v, char *vstrr){
    int i;
    printf("%s �ķ���:%d ָ��:%d β��: %d . ", vstrr, v->sign, v->exponent, v->mantissa[1]);    
    if(!dloatis0(v)){ 
        for(i=2; i<=Precision; i++)printf("%d ", v->mantissa[i]);
    }
    puts("");
}
*/

/* ������ fv ת�� Ϊ�߾��ȸ����� v */
#define BASE 256.0
void float2d(float fv, struct dloat *v){
    double tmp, fe;
    int i,j;
    if ( fv < 0 ){
        fv = -fv;
        v->sign = -1;
    } else v->sign = 1;
    if ( fv == 0 ){
        letitbe0(v);
        return;
    }
    tmp = log((double) fv)/log( BASE );
    fe = floor (tmp);
    tmp = exp((tmp-fe)*log( BASE ));
    v->exponent = (long int) fe;
    for (j=1 ; j<= Precision ; j++) {
        i=(int) tmp;
        v->mantissa[j]=(unsigned char) i;
        tmp = BASE *(tmp-i);
    }
}


#define MF 4
#define BI (1.0/256)
/* �߾��ȸ����� v ת��Ϊ������ */
float dloat2f(struct dloat *v){
    int j,mm;
    float fu,fv;
    if(dloatis0(v)) return 0.0;
    mm = (MF>Precision) ? Precision : MF;
    fv =(float) v->mantissa[mm];  
    for (j=mm-1;j>=1;j--) {
        fv *= BI;
        fv += v->mantissa[j];
    }
    fu = exp( log( BASE )*((double) v->exponent) );
    return ((float) fu*fv*(v->sign));
}

/*���Ը������͸߾��ȸ��������໥ת��*/
/*
#include <stdio.h>
void main(){
    struct dloat *v;
    float x=0;//1.E-16;
    int i;
    v = create_dloat();
    float2d(x, v);
    showd(v);
    printf("x:%f f:=%f\n", x, dloat2f(v));
    free_dloat(v);
}
*/

/*�߾��ȸ����� v �ĵ��� --> u*/
void dinv(struct dloat *u, struct dloat *v)
{
    unsigned char *r;
    assert(!dloatis0(v));//��0��������
    u->sign = v->sign;
    u->exponent = -v->exponent-1;
    if(u->exponent<=MIN_LMT){
        letitbe0(u);
        return;
    }
    r = (unsigned char *)malloc((Precision+2)*sizeof(unsigned char));
    d_inner_inv(r, v->mantissa, Precision+1, Precision);
    if(r[1])
    {//r=1
         d_inner_mov(u->mantissa, r, Precision); 
         u->exponent++;
    } else d_inner_mov(u->mantissa, &(r[1]), Precision); //r<1
  free(r);
    assert(u->exponent<=MAX_LMT);//���־���ֵ̫��
}


/*���� dinv*/
/*
#include <stdio.h>
void main(){
    struct dloat *v, *u;
    float x=-1.0/32000;
    int i;
    v = create_dloat();
    u = create_dloat();
    float2d(x, v);
    dinv(u, v);
    printf("x:%f v:=%f\n", x, dloat2f(v));
    showd(v);
    showd(u);
    printf("1/x:%f u:=%f\n", 1/x, dloat2f(u));
    free_dloat(u);
    free_dloat(v);
}
*/

//�����λ
//��λ����
//�޽�λ����
#define carry(r,w)    if((r)[1] != 0){\
  w->exponent += 1;\
  d_inner_mov(w->mantissa, (r), Precision);\
} else d_inner_mov(w->mantissa, &((r)[1]), Precision)


/*�߾��ȸ����� u*v --> w
�����������־λ is ���� -1��
����������is ���� 0
*/
void dmul(int *is, struct dloat *w, struct dloat *u, struct dloat *v){
    unsigned char *r;
    *is = 0;
    w->sign = u->sign * v->sign;
    w->exponent = u->exponent + v->exponent;
    r = (unsigned char *)malloc((2*Precision+1)*sizeof(unsigned char));
    d_inner_mul(r, u->mantissa, v->mantissa, Precision, Precision);
    carry(r,w);//if(r[1] != 0){//��λ
        //w->exponent += 1;
        //d_inner_mov(w->mantissa, r, Precision);
    //} else d_inner_mov(w->mantissa, &(r[1]), Precision);  
    free(r);
    if( dloatis0(w) ) letitbe0(w);
    if( w->exponent <= MIN_LMT ) letitbe0(w);
    if( w->exponent >= MAX_LMT ) *is = -1;
    return;
}

/*���� dmul*/
/*
void main(){
    struct dloat *v, *u, *w;
    float x=255, y=2550, z;
    int i, is;
    z=x*y;
    v = create_dloat();
    u = create_dloat();
    w = create_dloat();
    float2d(x, v);
    float2d(y, u);
    printf("x:%f v:=%f\n", x, dloat2f(v));
    showd(v);
    printf("y:%f u:=%f\n", y, dloat2f(u));
    showd(u);
    dmul(&is, w, u, v);
    assert(is!=-1);
    printf("x*y:%f w:=%f\n", z, dloat2f(w));
    showd(w);
    free_dloat(w);
    free_dloat(u);
    free_dloat(v);
}
//*/


/*�Ƚϸ߾��ȸ����� u v �ľ���ֵ
��|u|>=|v|����־λ is ���� 1; |u|<|v|������ 0
*/

int dabsbig(struct dloat *u, struct dloat *v){
    unsigned char *r;
    int is;
    if(u->exponent > v->exponent){
        is=1;
  } else if (u->exponent < v->exponent){
      is=0;
    } else {//u ��ָ�� = v ��ָ��
      r = (unsigned char *)malloc((Precision+1)*sizeof(unsigned char));
    d_inner_sub(&is, r, u->mantissa, v->mantissa, Precision);
    free(r);
    is += 1;
    }
    return is;
}



/*�߾��ȸ����� u+v --> w
�����������־λ is ���� -1��
����������is ���� 0
*/
void dadd(int *is, struct dloat *w, struct dloat *u, struct dloat *v){
    struct dloat *b, *s;
    long dif;
    int i,k;
    unsigned char *r1, *r2;

  if(dloatis0(v)){
      *is = 0 ;
      dcopy(w,u);
      return;
  }else if(dloatis0(u)){
       *is = 0 ;
      dcopy(w,v);
      return;
  }
   
    if(dabsbig(u,v)){
        b=u;
        s=v;
    } else{
        b=v;
        s=u;
    }//���ˣ�b,s ���� u,v������֤�� b �ľ���ֵ >= s �ľ���ֵ

     w->sign = b->sign;
     w->exponent = b->exponent;
    dif =    b->exponent - s->exponent;
    if( dif >= Precision ){ //���ھ������ޣ��ӷ�������� b
        *is = 0;
        dcopy(w,b);
        return;
    }
    
    r1 = (unsigned char *)malloc((Precision+1)*sizeof(unsigned char));
    for(i=1; i<=dif; i++) r1[i]=(unsigned char)0;
    d_inner_mov(&(r1[dif]), s->mantissa, Precision-dif);//����β��������ʹ�ö���ӷ�
    
    r2 = (unsigned char *)malloc((Precision+2)*sizeof(unsigned char));
    if( (s->sign)*(b->sign) == 1){ //b s ͬ��
      d_inner_add(r2, r1, b->mantissa, Precision);
      carry(r2,w);
      d_inner_mov(r2, w->mantissa, Precision);//һ���Ѿ�����������0+0.2�����λ�Ҫ�����빤��
    }else{//b s ���
          d_inner_sub(is, r2, b->mantissa, r1, Precision);
          assert(!*is);//��Ϊ b �ľ���ֵ >= s �ľ���ֵ
        }
        for(i=1;i<=Precision;i++)if(r2[i])break;//Ѱ�ҷ���λ
        k = i-1;// ǰ k λ����0,Ҫ��β����ǰ k λ��ָ��ҲҪ����
        if( (k>=Precision-1) && !r1[k+1]){//�����0
        letitbe0(w);
        }else{
          w->exponent -= k;
          d_inner_mov(w->mantissa, &(r2[k]), (Precision-k));
          for(i=Precision-k+1;i<=Precision;i++)w->mantissa[i]=(unsigned char)0;
        }
        free(r2);
        free(r1);
    *is = 0;
    if( w->exponent <= MIN_LMT ) letitbe0(w);
    if( w->exponent >= MAX_LMT ) *is = -1;
  return;
}



/*���� dadd*/
/*
void main(){
    struct dloat *v, *u, *w;
    float x=0.2,y=0,z;
        //x=2, y=-1,z;
    // x=1,y=-1,z;
    //y=1111111111, x=35670510543, z;
    int i, is;
    z=x+y;
    v = create_dloat();
    u = create_dloat();
    w = create_dloat();
    float2d(x, v);
    float2d(y, u);
    printf("x:%f v:=%f\n", x, dloat2f(v));
    showd(v);
    printf("y:%f u:=%f\n", y, dloat2f(u));
    showd(u);
    dadd(&is, w, u, v);
    assert(is!=-1);
    printf("x+y:%f w:=%f\n", z, dloat2f(w));
    showd(w);
    free_dloat(w);
    free_dloat(u);
    free_dloat(v);
}
*/


/*�Ը߾��ȸ����� u��|u|ԽСԽ�죩, exp(u) --> w */
void d_exp(struct dloat *w, struct dloat *u){
    struct dloat *tmp0,*tmp1,*tmp2,*tmp3;
    long i, dif;
    int is;
    tmp0 = create_dloat();
    tmp1 = create_dloat();
    tmp2 = create_dloat();
    tmp3 = create_dloat();
    int2d(1, tmp0);
    dcopy(tmp2, u);
    dadd(&is, tmp1, tmp2, tmp0);
    assert(is==0);
    dcopy(tmp0,tmp1);//u+1 -> tmp0
    dcopy(tmp1,u);
    for(i=2;;i++)
    {
        int2d(i, tmp2);
        dinv(tmp3, tmp2);
        dcopy(tmp2, tmp3);//tmp2=1/i
        dmul(&is, tmp3, tmp1, tmp2);//tmp3=tmp1/i
        assert(is==0);
        dmul(&is, tmp1, tmp3, u);// tmp1*(u/i) -> tmp1
        assert(is==0);
        dif =    tmp0->exponent - tmp1->exponent;
        if( dif >= Precision )  break; //�ﵽ������
        dadd(&is, w, tmp0, tmp1);//tmp1+tmp0 -> w
        dcopy(tmp0, w);
    }            
    dcopy(w,tmp0);
    free_dloat(tmp3);
    free_dloat(tmp2);
    free_dloat(tmp1);
    free_dloat(tmp0);
}


void dexp1(int *is, struct dloat *w, struct dloat *u){
    struct dloat *tmp1,*tmp2;
    long expt,j;
    int i; 
    *is = 0;
    if( dloatis0(u) ){
        int2d(1.0, w);
        return;
    }
    if(u->exponent < 0){ //��Ӧu<1
      d_exp(w, u);
      return;
    }
    expt =  u->exponent+1;
    tmp1 = create_dloat();
    dcopy(tmp1,u);
    tmp1->exponent=-1;
    tmp2 = create_dloat();
    d_exp(tmp2, tmp1);
    ///*
    for(i=1; i<=8; i++)
    {
//        printf("expt=%d\n",expt);//!
        for(j=1; j<=expt; j++)//ѭ�������⣡
        {
            dmul(is, tmp1, tmp2, tmp2);
            if(*is) break;
            if( tmp1->exponent <= MIN_LMT )
            {
                letitbe0(tmp1);
                break;
            }
            if( tmp1->exponent >= MAX_LMT )
            {
                *is = -1;
                break;
            }
            dcopy(tmp2, tmp1);
        }
        if(*is||dloatis0(tmp1)) break;
    }
    dcopy(w, tmp1);//*/
    free_dloat(tmp2);
    free_dloat(tmp1);
    return;
}


/*����dexp���Ȳ������⡣��Ϊ dexp �����ݴμ��㣬��
a(1+\delta)^n=a^n(1+n*\delta+0(\delta))
Ϊ��֤ dexp ����ľ��ȣ�Ҫ�Ӵ�������������
����struct dloat����ƣ����ǿ�������dloat.c ֻ��Ҫ����Precision
*/

/*�Ը߾��ȸ����� u, exp(u) --> w
�����������־λ is ���� -1��
����������is ���� 0
*/

void dexp(int *is, struct dloat *w, struct dloat *u){
    struct dloat **tmp;
    int i, default_precision,inc;
    default_precision=Precision;//��¼ȱʡ����
    if(Precision>=200)inc=log(Precision);
  else inc=2;
    inc+=(int)((8*fabs(u->exponent)+1)/256)+2;//���㾫��Ҫ���Ӷ��٣������㷨δ��׼ȷ
    Precision+=inc;
  tmp=create_dloat_vector(2);
    d_inner_mov(tmp[1]->mantissa, u->mantissa, default_precision);
    for(i=default_precision+1;i<=Precision;i++)tmp[1]->mantissa[i]=(unsigned char)0;
    tmp[1]->sign = u->sign;
    tmp[1]->exponent = u->exponent;
    dexp1(is, tmp[2], tmp[1]);
    w->sign = tmp[2]->sign;
    w->exponent = tmp[2]->exponent;
    d_inner_mov(w->mantissa, tmp[2]->mantissa, default_precision);
    free_dloat_vector(tmp,2);//free_dloat(inu);
    Precision=default_precision;
}


/*���� dexp*/
/*
void test(void){
    struct dloat *u, *w;
    float //x=-5e9;//0;//0;//00000;//..����
        //x=5e9; //�ٴ�������...
     x=-3.4;
    int is;
    u = create_dloat();
    w = create_dloat();
    float2d(x, u);
    printf("x:%f u:=%f\n", x, dloat2f(u));
    showd(u);
    dexp(&is, w, u);
  //d_exp(w, u);
  printf("is=%d\n",is);//!
    assert(is!=-1);
    printf("exp(x):%f w:=%f\n", exp(x), dloat2f(w));
    showd(w);
    free_dloat(w);
    free_dloat(u);
}
*/

/*
void main(){
test();
}
*/


#define show_dloat_as_float(s,u)      printf(s #u ":%f\n", dloat2f((u)))

//*���ɸ߾��ȸ������飬��n�����ݣ�����ʹ���±���1...n
struct dloat **create_dloat_vector(int n)
{
    struct dloat **v;
    int i;
    v = (struct dloat **)malloc((n+1)*sizeof(struct dloat*));
    for(i=1; i<=n; i++)    v[i]=create_dloat();
    return v;
}

void free_dloat_vector(struct dloat **v, int n)
{
    int i;
    for(i=n; i>=1; i--)    free_dloat(v[i]);
    free(v);
}

/* ͬά�� n �߾��ȸ������� a �� b �ڻ� --> w*/
void sum(struct dloat *w, struct dloat **a, struct dloat **b, int n)
{
    struct dloat **tmp;
    int i,is;
    tmp = create_dloat_vector(3);
    int2d(0, tmp[1]);
    for(i=1; i<=n;i ++)
    {
        dmul(&is, tmp[2], a[i], b[i]);
        assert(is==0);
        dadd(&is, tmp[3], tmp[2], tmp[1]);
        assert(is==0);
        dcopy(tmp[1], tmp[3]);
    }
    dcopy(w, tmp[1]);
    free_dloat_vector(tmp,3);
}


/*�Ѹ߾��ȸ����� b ���Ƶ� a , b --> a */
void dcopy(struct dloat *a, struct dloat *b)
{
    a->sign = b->sign;
    a->exponent = b->exponent;
  d_inner_mov(a->mantissa, b->mantissa, Precision);
}








//================================================
//hpn.c
//================================================

struct hpn *create_hpn(int prec){
    struct hpn *p;
    int default_precision;
    assert(prec>2);
    p = (struct hpn *)malloc(sizeof(struct hpn));
    default_precision=Precision;//��¼ȱʡ����
    Precision=prec;
    p->precision=Precision;
    p->value=create_dloat();
    Precision=default_precision;
    return p;
}

void free_hpn(struct hpn *p){
    free_dloat(p->value);
    free(p);
}

void int2hpn(int i, struct hpn *p){
    int default_precision;
    default_precision=Precision;//��¼ȱʡ����
    Precision=p->precision;
    int2d(i, p->value);
    Precision=default_precision;
}

/*��һ���ڲ��߾�������ʮ������������mλ*/
void myshow(struct dloat *v, int m){
    int i,j,is,shift,flag=0;
    char *out;
    long int exponent10=0;
    struct dloat **tmp;
  
  if(dloatis0(v)){
      Rprintf("0\n");
      return;
  }

    out = (char *)malloc((m+1)*sizeof(char));
  tmp = create_dloat_vector(9);
   
  if(v->sign<0){
    Rprintf("%s","-");
    v->sign=1;
  }  else   Rprintf("%s","+");/*������ţ���ȡ����ֵ*/
  
  int2d(10,tmp[3]);//tmp[3]=10
  dinv(tmp[4],tmp[3]);//tmp[4]=0.1

  if(v->exponent<0){//����С��1
      while(v->exponent <0){//���ϳ�10��ֱ�����ִ��ڻ����1
        dcopy(tmp[5],v);
        dmul(&is, tmp[8], tmp[5], tmp[3]);
        assert(is==0);
        exponent10 -= 1;//ͬʱ��¼10�ݴεı仯
        dcopy(v,tmp[8]);
      }
  }
    
  if(v->exponent>=0){//���ִ��ڻ����1
      while(v->exponent >= 0){//���ϳ�0.1��ֱ������С��1
        dcopy(tmp[5],v);
        dmul(&is, tmp[8], tmp[5], tmp[4]);
        assert(is==0);
        exponent10 += 1;//ͬʱ��¼10�ݴεı仯
        dcopy(v,tmp[8]);
      }
  }
  
  //��ʱ�������� [0.1, 1)
  //ע���п�����0.999999����Ŀǰ��λ�µĲ��ϵ�255
  
  for(i=0;i<=m;i++)out[i]=(char) 0; //����ĳ�ʼ�� 
    for(i=1;i<=m;i++){
      dmul(&is,tmp[1],v,tmp[3]); /*�����ֲ��ϳ���10*/
      assert(is==0);
    dcopy(v,tmp[1]);
        if(!flag && v->exponent==0){
            if((int)v->mantissa[1]==10){
              flag=1;
              out[i-1]=(char) (1+(int)out[i-1]);
            }else{
              out[i]=v->mantissa[1];
            }
            dcopy(tmp[6],v);
            for(j=2;j<=Precision;j++)tmp[6]->mantissa[j]=(unsigned char)0;/*tmp[6]�൱�ڽ���������һλ��v*/
            dneg(tmp[6]);
            dadd(&is,tmp[1],tmp[6],v);
            assert(is==0);
            dcopy(v,tmp[1]);
        }
    }
    
    if(out[0])shift=1; else shift=0;
  Rprintf("%1d.",(int)out[1-shift]);
  for(i=2;i<=m;i++) Rprintf("%1d",(int)out[i-shift]);
  Rprintf("E%+d",(int)exponent10-1+shift);

    free_dloat_vector(tmp,9);
    free(out);
}

/*��һ���߾�������ʮ������������mλ*/
void show_hpn(struct hpn *p, int m){
    int default_precision;
  default_precision=Precision;
  Precision=p->precision;
    assert(m>1);
    myshow(p->value, m);
    Precision=default_precision;
}

/*�Ը߾��ȸ����� u, exp(u) --> w
�����������־λ is ���� -1��
����������is ���� 0
*/
void hpn_exp(int *is, struct hpn *w, struct hpn *u){
    int default_precision;
    assert(w->precision==u->precision);
  default_precision=Precision;
  Precision=w->precision;
    dexp(is, w->value, u->value);
    Precision=default_precision;
}


//==================================================
//test.cpp
//==================================================
double mymain(double x){
    double y;
    long double lx,ly;
    struct hpn *tmp1,*tmp2;
    lx=(long double)x;
    tmp1=(struct hpn *)create_hpn(42);
  tmp2=(struct hpn *)create_hpn(42);
  int2hpn((int)(x),tmp1);//�ڲ������exp(x��ȡ��)
    //hpn_exp(&is, tmp2, tmp1);
    //assert(is==0);
    ly=expl(lx);
    y=exp(x);
//    Rprintf("  ly=%50.24Lf\n",ly);
    //Rprintf("   y=%50.24f\n",y);
  //show_hpn(tmp2,100);
    //Rprintf("\n");
    free_hpn(tmp2);
    free_hpn(tmp1);
        return y;
}
