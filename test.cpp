#include "Lequation.h"

int main(){
 Lequation le;
 cout<< "Equation of linear: "<<endl;
 le.linear_equation(1,-1);
 cout<< "Equation of quadratic: "<<endl;
 le.quadratic_equation(1,1,1);
 cout<< "Equation of cubes: "<<endl;
 le.cub(2,-11,-12,-9);
 cout<< "Equation of forth: "<<endl;
 le.forth_equation(0, 212, -1121, -123, -100);
 return 0;
}
