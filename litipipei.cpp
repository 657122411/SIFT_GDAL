/*==========================================*/   
/*             最小二乘匹配                  */   
/*                             1995。11。16 */   
/*==========================================*/   
#include <stdlib.h>   
#include <malloc.h>   
#include <stdio.h>   
#include <math.h>   
   
   
/*==================================================*/   
/*               高斯解方程子程序  
a[n1][n1]*x[n1][1]=b[n1][1]   AX=B   */   
/*==================================================*/   
int gs(float *a,float *b,int n1,float *x)   
{int *js,l,k,i,j,is,p,q;   
float d,t;   
js=(int *)malloc(n1*sizeof(int));   
l=1;   
for(k=0;k<=n1-2;k++)   
{d=0.0;   
for(i=k;i<=n1-1;i++)   
for(j=k;j<=n1-1;j++)   
{t=(float)fabs(*(a+i*n1+j));   
if(t>d) {d=t;js[k]=j;is=i;}   
}     
if(d+1.0==1.0) l=0;   
else   
{if(js[k]!=k)   
for(i=0;i<=n1-1;i++)   
{p=i*n1+k;q=i*n1+js[k];   
t=*(a+p);*(a+p)=*(a+q);*(a+q)=t;   
}   
if(is!=k)   
{for(j=k;j<=n1-1;j++)   
{p=k*n1+j;q=is*n1+j;   
t=*(a+p);*(a+p)=*(a+q);*(a+q)=t;   
}   
t=*(b+k);*(b+k)=*(b+is);*(b+is)=t;   
}   
}   
if(l==0)   
return(0);   
d=*(a+k*n1+k);   
for(j=k+1;j<=n1-1;j++)   
{p=k*n1+j;*(a+p)=*(a+p)/d;}   
*(b+k)=*(b+k)/d;   
for(i=k+1;i<=n1-1;i++)   
{for(j=k+1;j<=n1-1;j++)   
{p=i*n1+j;   
*(a+p)=*(a+p)-*(a+i*n1+k)*(*(a+k*n1+j));   
}   
*(b+i)=*(b+i)-*(a+i*n1+k)*(*(b+k));   
}   
}   
d=*(a+(n1-1)*n1+n1-1);   
if(fabs(d)+1.0==1.0)   
{free(js);   
return(0);   
}   
*(x+n1-1)=*(b+n1-1)/d;   
for(i=n1-2;i>=0;i--)   
{t=0.0;   
for(j=i+1;j<=n1-1;j++)   
t=t+*(a+i*n1+j)*(*(x+j));   
*(x+i)=*(b+i)-t;   
}   
js[n1-1]=n1-1;   
for(k=n1-1;k>=0;k--)   
if(js[k]!=k)   
{t=*(x+k);*(x+k)=*(x+js[k]);*(x+js[k])=t;}   
free(js);   
return(1);          //有解   
}   
   
   
/*++++++++++++++++++++++++++++++++++++++++++++++++++*/   
/*              最小二乘匹配子程序                    */   
/*                                                  */   
/*++++++++++++++++++++++++++++++++++++++++++++++++++*/   
   
void LSM_Match(unsigned char *dstim,unsigned char *srcim,int XO1,int YO1,int XO2,int YO2,int lWidth,float&xmatch,float&ymatch,int m)   
{   
   
   
    int ix = lWidth;            //传进变量   
    int i,j,p,q,solut;          //solut 是否有解   
    int Gx,Gy;   
    int k=0,k1,k2,k3,k4;   
    int i1,i2,ii2,j1,j2/*,i3,j3*/;   
    static float a[8]={0,1,0,0,0,1,0,1};        //a0,a1,a2,b0,b1,b2,h0,h1      
    float B20,B10;   
    float *Xr,*Yr /*,*L*/; //At[8][N*N];                //Xr>X 坐标改正，Yr>Y坐标改正   
    float A[9],f1[8][8],g1[8][1],x[8][1];           //f1,C'C   g1 C'L      x方程的解   
       
    double min1=1.0;/*,min2=1.50;*/   
    //float X,Y;   
       
    k1= YO1-m/2;      k3= XO1-m/2;   //i,j;         //匹配起点(k3,k1),匹配中心(X01,Y01)   
    k2= YO2- 51/2;    k4= XO2- 51/2; //i,j;         //搜索起点(k4,k2),搜索中心(X02,Y02)   
       
    unsigned char *F = new unsigned char[m*m];   //m是虾米东西？匹配窗口的大小 比如5*5 或 3*3   
    unsigned char *G = new unsigned char[m*m];      //变化，需要重采样   
    unsigned char *GO= new unsigned char[51*51]; // allocate the search image buffer. 分配搜索图像的缓冲区,保持不变   
    if( F == NULL || G == NULL || GO == NULL)   
        return;   
    Xr = new float[m*m];      Yr = new float[m*m];   
    //  L = new float[m*m];   
    if(Xr == NULL || Yr == NULL /*|| L == NULL || At == NULL*/)   
        return;   
       
    /* 读入匹配窗口的数据 */   
    for(i=0; i<m; i++)   
        for(j=0; j<m; j++)   
        {   
            Xr[i*m+j]= 0.;       Yr[i*m+j]= 0.;   
            F[i*m+j]= dstim[(unsigned long)(k1+i)*ix+k3+j]; //目标图像          //gfs ix未定义，猜测是图像化4整后的宽   
        }   
    for(i=0; i<51; i++)   
        for(j=0; j<51; j++)   
            GO[i*51+j] = srcim[(unsigned long)(k2+i)*ix+k4+j];  //搜索图像,保持不变，变的是G            //gfs ix也为搜索图像的宽   
           
    Xr[m*m/2] = 0;    Yr[m*m/2] = 0;                //多余 ，中心位置也为0   
    k2 = 25-m/2;      k4 = 25-m/2;                  //GO中心的m*m矩左上   
       
    for(i=0;i<8;i++)   
    {   
        a[i]=0;   
        if(i==1 ||i==5 || i==7)  a[i]=1;            //a0~a2,b0~b2,h0,h1   
    }   
   
    for(ii2=0; ii2 < 5; ii2 ++)   //  max loop number =5   
    {   
        //临近差值，重新计算搜索起点位置k4,k2   
        if(Xr[m*m/2]>= 0.5)  k4=k4+(int)(Xr[m*m/2]+0.5);   
        if(Xr[m*m/2]<=-0.5)  k4=k4+(int)(Xr[m*m/2]-0.5);   
        if(Yr[m*m/2]>= 0.5)  k2=k2+(int)(Yr[m*m/2]+0.5);   
        if(Yr[m*m/2]<=-0.5)  k2=k2+(int)(Yr[m*m/2]-0.5);   
           
        for(i=0;i<m;i++)   
            for(j=0;j<m;j++)   
            {   
                i1 = int( i+ Xr[i*m+j]+ 0.5);           //四舍五入取整用   
                j1 = int( j+ Yr[i*m+j]+ 0.5);   
                G[i*m+j]=GO[(i1+k2)*51+j1+k4];      //完成重采样（2）   
            }   
               
            //初值为0   
            for(i=0;i<8;i++)   
            {   
                g1[i][0]=0.;   
                for(j=0;j<8;j++)   
                    f1[i][j]=0.;   
            }   
               
            B20=a[7];           //h1   
            B10=a[6];           //h0   
            k=0;   
               
            /*   法化组成法方程*/   
            for(i=0;i<m;i++)   
                for(j=0;j<m;j++)   
                {   
                    Xr[k]= i-m/2.0f;                //像素中心化   
                    Yr[k]= j-m/2.0f;   
                       
                    i1=(i-1<0)?i:i-1;            //防止差分越界，以后无用   
                    i2=(i+1>m-1)?i:i+1;   
                    j1=(j-1<0)?j:j-1;   
                    j2=(j+1>m-1)?j:j+1;   
                    Gx=(G[i*m+j2]-G[i*m+j1])/2;     //差分运算   
                    Gy=(G[i2*m+j]-G[i1*m+j])/2;   
                       
                    A[0]=-B20*Gx;                   //A[0]=-h1*g'x   
                    A[1]=-B20*Gx*Xr[k];             //A[1]=-h1*g'x*x   
                    A[2]=-B20*Gx*Yr[k];             //A[2]=-h1*g'x*y   
                    A[3]=-B20*Gy;                   //A[3]=-h1*g'y*x   
                    A[4]=-B20*Gy*Xr[k];             //A[4]=-h1*g'y*y   
                    A[5]=-B20*Gy*Yr[k];             //A[5]=-h1*g'y*y   
                    A[6]=-1.0;   
                    A[7]=-G[i*m+j]*1.0f;                    //A[7]=-g   
                    A[8]=B10+B20*G[i*m+j]-F[i*m+j];         //A[8] = h0+h1*g-f;   
   
                    k=k+1;   
                    for(p=0;p<8;p++)   
                    {   
                        g1[p][0]=g1[p][0]+A[p]*A[8];            //C'L   
                        for(q=p;q<8;q++)   
                            f1[p][q]=f1[p][q]+A[p]*A[q];   
                    }   
                    for(p=0;p<8;p++)   
                        for(q=p;q<8;q++)   
                            f1[q][p]=f1[p][q];          //对称阵 C'C   
                           
                }   
                solut=gs(*f1,*g1,8,*x);         /*  高斯求解改正数值（5）  */   
                if(solut)                       //有解   
                {   
                    for(i=0;i<8;i++)   
                    { a[i]=a[i]+x[i][0];}           //对a进行改正（6）   
                    min1=0;   
                    for(i=0;i<8;i++)   
                        if(fabs(x[i][0])>min1) min1=fabs(x[i][0]);  //求x最大值   
                }   
                   
                k=0;   
                for(i=0;i<m;i++)          /*    坐标改正 （1）    */   
                    for(j=0;j<m;j++)   
                    {   
                        Xr[k]=a[0]+a[1]*Xr[k]+a[2]*Yr[k];   
                        Yr[k]=a[3]+a[4]*Xr[k]+a[5]*Yr[k];   
                        k=k+1;   
                    }   
                    xmatch= XO2 + Xr[m*m/2];           /*    重取样    */  //计算改正后中心匹配位置（7）   
                    ymatch= YO2 + Yr[m*m/2];                                   
                       
                    if(Xr[m*m/2]>= 0.5) XO2=XO2+(int)(Xr[m*m/2]+0.5);   
                    if(Xr[m*m/2]<=-0.5) XO2=XO2+(int)(Xr[m*m/2]-0.5);   
                    if(Yr[m*m/2]>= 0.5) YO2=YO2+(int)(Yr[m*m/2]+0.5);   
                    if(Yr[m*m/2]<=-0.5) YO2=YO2+(int)(Yr[m*m/2]-0.5);   
                    if(min1<0.0005)                  /*   控制循环域值设定    */   
                        break;   
    }   
       
    delete []F;  delete []G;  delete []GO;   
    delete []Xr; delete []Yr; /*delete At;*/   
}   
   
   