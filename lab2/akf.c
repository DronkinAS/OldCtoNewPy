#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define M 4000


int main ()
{
double step, *t, *akf, temp_sum=0, srednee, koeff, dz, mod_akf;  
double *x, *ogib, *time_ogib, temp, temp1, temp2, temp3, temp4,sum1=0,sum2=0;
double Integral;
int i=0,k=0, N=0, m, ind_ogib=0;
FILE *in_date;
FILE *out_date;
FILE *out_ogib_date;
FILE *in_ogib_date;
//РѕС‚РєСЂС‹РІР°РµРј С„Р°РёР» РґР»СЏ С‚РѕРіРѕ, С‡С‚РѕР±С‹ СѓР·РЅР°С‚СЊ РєРѕР»РёС‡РµСЃС‚РІРѕ С‚РѕС‡РµРє РІ СЂРµР°Р»РёР·Р°С†РёСЏС…
  in_date = fopen ("realization.dat","r");
//  printf ("РЎС‡РёС‚С‹РІР°СЋ РєРѕР»РёС‡РµСЃС‚РІРѕ С‚РѕС‡РµРє РІ СЂРµР°Р»РёР·Р°С†РёРё \n");
  while (!feof(in_date))
  {
     fscanf (in_date, "%lf %lf", &temp, &temp2);
     N=N+1 ;
  }
  fclose (in_date);
  N=N-1;
//  printf ("РљРѕР»РёС‡РµСЃС‚РІРѕ С‚РѕС‡РµРє РІ СЂРµР°Р»РёР·Р°С†РёРё СЂР°РІРЅРѕ %d \n", N);
//  printf ("Р’С‹РґРµР»СЏСЋ РїР°РјСЏС‚СЊ \n");
//РІС‹РґРµР»СЏРµРј РїР°РјСЏС‚СЊ РґР»СЏ РєР°Р¶РґРѕРіРѕ РјР°СЃСЃРёРІР°
  x = (double *)calloc(N, sizeof(double));
  t = (double *)calloc(2, sizeof(double));
  akf = (double *)calloc(M+1, sizeof(double));
//  printf ("РЈСЃРїРµС€РЅРѕ \n");

//РѕС‚РєСЂС‹РІР°РµРј С„Р°РёР» РґР»СЏ СЃС‡РёС‚С‹РІР°РЅРёСЏ РІС…РѕРґРЅС‹С… РґР°РЅРЅС‹С… Рё СЃРѕР·РґР°РЅРёСЏ РјР°СЃСЃРёРІР° РґР»СЏ РѕР±СЂР°Р±РѕС‚РєРё
  in_date = fopen ("realization.dat","r");

//РѕРїСЂРµРґРµР»СЏРµРј С€Р°Рі РІСЂРµРјРµРЅРЅРѕР№ СЂРµР°Р»РёР·Р°С†РёРё
//  printf ("РћРїСЂРµРґРµР»СЏСЋ С€Р°Рі РІСЂРµРјРµРЅРЅРѕР№ СЂРµР°Р»РёР·Р°С†РёРё \n");  
  for (i=0; i<2; i++)
  {
   fscanf (in_date, "%lf %lf", &t[i], &temp);
  }
  step=t[1]-t[0];
//  printf ("РЁР°Рі РІСЂРµРјРµРЅРЅРѕР№ СЂРµР°Р»РёР·Р°С†РёРё СЂР°РІРµРЅ %f\n", step);
  fclose (in_date);
  
//РЎРѕР·РґР°РµРј РјР°СЃСЃРёРІ, СЃРѕРґРµСЂР¶Р°С‰РёР№ РІ СЃРµР±Рµ РІСЂРµРјРµРЅРЅСѓСЋ СЂРµР°Р»РёР·Р°С†РёСЋ, РїРѕРґСЃС‡РёС‚С‹РІР°РµРј СЃСЂРµРґРЅРµРµ Р·РЅР°С‡РµРЅРёРµ
  in_date = fopen ("realization.dat","r");
  while (!feof(in_date))
  {
     fscanf (in_date, "%lf %lf", &temp, &x[k]);
     temp_sum=temp_sum+x[k];
     k=k+1 ;
  }
  fclose (in_date);
  
  srednee=temp_sum/N;
//  printf ("РЎСЂРµРґРЅРµРµ Р·РЅР°С‡РµРЅРёРµ <x(t)> = %f\n", srednee);
  
//РџСЂРёСЃС‚СѓРїР°РµРј Рє СЂР°СЃС‡РµС‚Сѓ Р°РІС‚РѕРєРѕСЂСЂРµР»СЏС†РёРѕРЅРЅРѕР№ С„СѓРЅРєС†РёРё Рё РµРµ РјРѕРґСѓР»СЏ
//РђР»РіРѕСЂРёС‚Рј Рё С„РѕСЂРјСѓР»Сѓ СЃРј. РІ РєРЅРёРіРµ Р”.Р­. РџРѕСЃС‚РЅРѕРІ, Рђ.Рќ. РџР°РІР»РѕРІ, РЎ.Р’. РђСЃС‚Р°С…РѕРІ "РњРµС‚РѕРґС‹ РЅРµР»РёРЅРµР№РЅРѕР№ РґРёРЅР°РјРёРєРё", СЃС‚СЂ. 108  
  out_date = fopen ("akf.dat","w");
  for (m=0; m<=M; m++)
  { 
    temp_sum=0;
    k=N-m;
    for(i=0; i<=k; i++)
        {
          koeff = (x[i]-srednee)*(x[i+m]-srednee);
          temp_sum=temp_sum+koeff;
        }
    akf[m]=temp_sum/k;
    
    fprintf (out_date,"%f\t%f\n", m*step,akf[m]/akf[0]);
  }
  
  fclose (out_date);
  
//РџСЂРёСЃС‚СѓРїР°РµРј Рє СЂР°СЃС‡РµС‚Сѓ РёРЅС‚РµРіСЂР°Р»Р° РѕС‚ РјРѕРґСѓР»СЏ Р°РІС‚РѕРєРѕСЂСЂРµР»СЏС†РёРѕРЅРЅРѕР№ С„СѓРЅРєС†РёРё 
//Рё РѕРєРѕРЅС‡Р°С‚РµР»СЊРЅРѕРіРѕ РїРѕРґСЃС‡РµС‚Р° РІСЂРµРјРµРЅРё РєРѕСЂСЂРµР»СЏС†РёРё  
//Р?СЃРїРѕР»СЊР·СѓРµРј РјРµС‚РѕРґ СЃРёРјРїСЃРѕРЅР°  
  for (i=2; i<=M; i=i+2)
  {  
    sum1=sum1+fabs(akf[i-1]);
  }
  
  for (i=2; i<=M-2; i=i+2)
  {  
    sum2=sum2+fabs(akf[i]);
  }
  
  Integral=step*(fabs(akf[0])+4*sum1+2*sum2+fabs(akf[M]))/3;
  printf("correlation time =%f\n", Integral/akf[0]);

  
  free(x);
  free(t);
  free(akf);

return(0);
} 
