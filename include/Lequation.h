#include <iostream>
#include <cmath>
#include <iostream>
#include <stdio.h>
#include <conio.h>
#include <locale.h>
#ifndef LEQUATION_H
#define LEQUATION_H
using namespace std;


class Lequation{
private:
    double root;
    double e;
    double pi;
public:
    Lequation(){
        pi = 3.14159265358979323846;
        e = 2.71828182846;
    }

    /* Створення допоміжної ф-ії sign */

    template <typename T> double sgn(T val) {
    return (T(0) < val) - (val < T(0));
    }

    /* Створення допоміжної функії arsh */

    double arsh(double x){
        return log(x+sqrt(x*x+1));
    }

    /* Створення допоміжної функії arch */

    double arch(double x){
        return log(x+sqrt((x+1)*(x-1)));
    }

    /* Створення допоміжної функії ch */

    double ch(double x){
        return (pow(e, x)-pow(e, -x))*1./2;
    }

    /* Створення допоміжної функії sh */

    double sh(double x){
        return (pow(e, x)+pow(e, -x))*1./2;
    }

    /* Функція для розв*язання лінійного рівняння  */

    double linear_equation(double a, double b){
        if (a==-b==0){
            cout<<"This equation has an infinite number of solutions."<<endl;
        }
        else if (a==0 & b!=0){
            cout<<"This equation has not have any solutions."<<endl;
        }
        else if (a!=0){
            root = (-b)/a;
            cout<<"Solutions of this equation is "<<root<<endl;
        }
        return 0;
    }

    /* Фунція для розв*язання квадратичного  рівнянння */

    double quadratic_equation(double a, double b, double c){
        float x1, x2, discriminant, realPart, imaginaryPart;
        discriminant = b*b - 4*a*c;
        if (discriminant > 0) {
            x1 = (-b + sqrt(discriminant)) / (2*a);
            x2 = (-b - sqrt(discriminant)) / (2*a);
            cout << "Roots are real and different." << endl;
            cout << "x1 = " << x1 << endl;
            cout << "x2 = " << x2 << endl;
        }
        else if (discriminant == 0) {
            cout << "Roots are real and same." << endl;
            x1 = (-b + sqrt(discriminant)) / (2*a);
            cout << "x1 = x2 =" << x1 << endl;
        }
        else {
            realPart = -b/(2*a);
            imaginaryPart =sqrt(-discriminant)/(2*a);
            cout << "Roots are complex and different."  << endl;
            cout << "x1 = " << realPart << "+" << imaginaryPart << "i" << endl;
            cout << "x2 = " << realPart << "-" << imaginaryPart << "i" << endl;
    }
    }

    /* Функція для ров*язання кубічного рівняння */

    double cub(double a, double b, double c, double d){
        double p,q,s,f,x1,x2,x3,x2i,x3i;
        p=((3.*a*c-b*b)/(3.*a*a));
        q=((2.*b*b-9.*a*b*c+27.*a*a*d));
        s=(((q*q)/4.)+(p*p*p)/27.);
        if (q<0){
            f=(atan(pow(-s,0.5)/(-q/2)));
            }
        else if (q>0){
            f=(atan(pow(-s,0.5)/(-q/2))+pi);
            }
        else{
            f=(pi/2);
            }
        if (s<0){
            x1=(2.*pow((-p/3.),0.5)*cos(f/3.)-b/3.*a);
            x2=(2.*pow((-p/3.),0.5)*cos(f/3.+(2.*pi)/3.)-b/3.*a);
            x3=(2.*pow((-p/3.),0.5)*cos(f/3.+(2.*pi)/3.)-b/3.*a);
        }
        else if (s>0){
            x1=(pow(-q/2.+pow(s,0.5),1./3.)+pow(-q/2.-(pow(s,0.5)),1./3.)-b/(3.*a));
            x2=(-0.5*(pow(-q/2.+pow(s,0.5),1./3.)+pow(-q/2.-(pow(s,0.5)),1./3.)-b/(3.*a)));
            x2i=((pow(3.,0.5)/2.)*(pow(-q/2.+(pow(s,0.5)),1./3.)-pow(-q/2.-(pow(s,0.5)),1./3.)));
            x3=(-0.5*(pow(-q/2+pow(s,0.5),1/3)+pow(-q/2-(pow(s,0.5)),1/3)-b/(3*a)));
            x3i=((pow(3.,0.5)/2.)*(pow(-q/2.+(pow(s,0.5)),1./3.)-pow(-q/2.-(pow(s,0.5)),1./3.)));
        }
        else{
            x1=(2.*pow(-q/2.,1./3.)-b/(3.*a));
            x2=(-1.*pow(-q/2.,1./3.)-b/(3.*a));
            x3=(-1.*pow(-q/2.,1./3.)-b/(3.*a));
        }
        cout<<"x1= "<<x1<<endl;
        if (s>0){
            cout<<"x2= "<<x2<<"+"<<x2i<<"i"<<endl;
            cout<<"x3= "<<x3<<"-"<<x3i<<"i"<<endl;
        }
        else{
            cout<<"X2= "<<x2<<endl;
            cout<<"X3= "<<x3<<endl;
        }
    }

       /* Функція для ров*язання рівняння 4-ого степеня  */


    double forth_equation(double a, double b, double c, double d, double e){
        double x1,x3,x2,x4,p,q, del0, del1,S, Q;
        if (a==0){
            cub(b,c,d,e);
        }
        else{
            del0 = c*c-3*b*d+12*a*e;
            del1 = 2*c*c*c-9*b*c*d+27*b*b*e+27*a*d*d-72*a*c*e;
            p = (8*a*c-3*b*b)/(8*a*a);
            q = (b*b*b-4*a*b*c+8*a*a*d)/(8*a*a*a);
            Q = pow(del1+sqrt(del1*del1-4*del0*del0*del0)/2, 1/3);
            S = 1/2*sqrt(-2/3*p+1/(3*a)*(Q+del0/Q));
            x1 = -b/(4*a)-S+1/2*sqrt(-4*S*S-2*p+q/S);
            x2 = -b/(4*a)-S-1/2*sqrt(-4*S*S-2*p+q/S);
            x3 = -b/(4*a)+S-1/2*sqrt(-4*S*S-2*p+q/S);
            x4 = -b/(4*a)+S+1/2*sqrt(-4*S*S-2*p+q/S);
            cout<<"x1 = "<< x1 <<endl;
            cout<<"x2 = "<< x2 <<endl;
            cout<<"x3 = "<< x3 <<endl;
            cout<<"x4 = "<< x4 <<endl;
        }
    }
};

#endif // LEQUATION_H
