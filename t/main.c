#include <iostream> 
#include <math.h> 
#include <fstream> 
#include <iomanip> 
#include <stdio.h> 
#include <stdlib.h> 

using namespace std; 

double k1[3], k2[3], k3[3], k4[3], x=1, y=2, t=200, z=3, h=0.001; 
double f1 (double y, double z, double x) 
{ 
double m=2.0; 
return m*x+y-x*z; 
} 
double f2 (double x) 
{ 
return -x; 
} 
double f3 (double x, double z) 
{ double g=0.2; 

return -g*z+g*x*x;} 

int main(int argc, char** argv) 

{setlocale(LC_ALL,"Russian"); 
ofstream G; 
G.open("C:\\sites\\file_1.txt", ios::out); 

for (t=0; t<=200; h=t+h) 
{ 
k1[0]=h*f1(y,z,x); 
k1[1]=h*f2(x); 
k1[2]=h*f3(x,z); 

k2[0]=h*f1(y+0.5*k1[1],z+0.5*k1[2],x+0.5*k1[0]); 
k2[1]=h*f2(x+0.5*k1[0]); 
k2[2]=h*f3(x+0.5*k1[0],z+0.5*k1[2]); 

k3[0]=h*f1(y+0.5*k2[1],z+0.5*k2[2],x+0.5*k2[0]); 
k3[1]=h*f2(x+0.5*k2[0]); 
k3[2]=h*f3(x+0.5*k2[0],z+0.5*k2[2]); 

k4[0]=h*f1(y+k3[1],z+k3[2],x+k3[0]); 
k4[1]=h*f2(x+k3[0]); 
k4[2]=h*f3(x+k3[0],z+k3[2]); 

x=x+(k1[0]+2*k2[0]+2*k3[0]+k4[0])/6; 
y=y+(k1[1]+2*k2[1]+2*k3[1]+k4[1])/6; 
z=z+(k1[2]+2*k2[2]+2*k3[2]+k4[2])/6; 
if (t>100) { 
G«x«" "«y«" "«z«endl;} 
} 

G.close(); 

system("pause"); 
return 0; 
}