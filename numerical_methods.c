#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <stdbool.h>
#include <string.h>
#define N 100

int factorial(int n){
	int i, sum = 1;
	if(n < 0) return -1;
	else if(n == 1) return 1;
	else{
		for(i = n; i > 0; i--){
			sum *= i;
		}
	}
}


double polynom(int degree, double coefs[], double x, double constant){
    int i;
    int power;
    double sum = 0;
    for(i = 0; i < degree; i++){
        power = degree - i;
        sum += pow(x, power) * coefs[i];
    }
    sum += constant;
    return sum;
}

double absolute(double val){
	if(val < 0) val = -val;
	return val;
}

double merkeziturev(int degree, double coefs[], double x, double constant, double h){
  double y = (polynom(degree, coefs, x+h, constant) - polynom(degree, coefs, x-h, constant))/(2*h);
  return y;
}

double geriturev(int degree, double coefs[], double x, double constant, double h){
  double y = (polynom(degree, coefs, x, constant) - polynom(degree, coefs, x-h, constant))/(h);
  return y;
}

double ileriturev(int degree, double coefs[], double x, double constant, double h){
  double y = (polynom(degree, coefs, x+h, constant) - polynom(degree, coefs, x, constant))/(h);
  return y;
}

void bisection(int degree, double coefs[], double constant, double error){
	int i;
	double y;
	bool isFound = false;
	bool isRootThere = false;
	double leftPoint;
    double rightPoint;
    double midPoint;
	while(!isRootThere){
    	printf("araligin sol tarafini giriniz: \n");
	    scanf("%lf", &leftPoint);
	    printf("araligin sag tarafini giriniz: \n");
	    scanf("%lf", &rightPoint);
	    if(polynom(degree, coefs, leftPoint, constant) * polynom(degree, coefs, rightPoint, constant) < 0) isRootThere = true;
	    else{
	    	printf("aralikta kok yok, lutfen bir daha giriniz.\n");
		}
	}
	i = 1;
	while(!isFound){
		midPoint = (leftPoint + rightPoint)/2;
		y = polynom(degree, coefs, midPoint, constant);
		if((rightPoint - leftPoint)/pow(2,i) < error) isFound = true;
		else{
			if(polynom(degree, coefs, midPoint, constant) * polynom(degree, coefs, leftPoint, constant) < 0) rightPoint = midPoint;
			else leftPoint = midPoint;
		}
		i++;
	}
	if(isFound) printf("kok: %lf", midPoint);
}

void RegulaFalsi(int degree, double coefs[], double constant, double error ){
	int i;
	double leftPoint;
	double rightPoint;
	double c;
	bool isFound = false;
	bool isRootThere = false;
	while(!isRootThere){
    	printf("araligin sol tarafini giriniz: \n");
	    scanf("%lf", &leftPoint);
	    printf("araligin sag tarafini giriniz: \n");
	    scanf("%lf", &rightPoint);
	    if(polynom(degree, coefs, leftPoint, constant) * polynom(degree, coefs, rightPoint, constant) < 0) isRootThere = true;
	    else{
	    	printf("aralikta kok yok, lutfen bir daha giriniz.\n");
		}
	}
	c = (polynom(degree, coefs, leftPoint, constant)*rightPoint - polynom(degree, coefs, rightPoint, constant)*leftPoint)/(polynom(degree, coefs, leftPoint, constant)-polynom(degree, coefs, rightPoint, constant));

	i = 1;
	while(!isFound){
		if(absolute(rightPoint-leftPoint)/pow(2,i) < error){
			printf("kok: %lf\n", c);
			isFound = true;
		}
		if(polynom(degree, coefs, leftPoint, constant)*c < 0){
			rightPoint = c;
		}
		else if(polynom(degree, coefs, leftPoint, constant)*c > 0){
			leftPoint = c;
		}
		else{
			printf("%lf ya da %lf kok", leftPoint, rightPoint);
			isFound = true;
		}
		i++;		
		
		c = (polynom(degree, coefs, leftPoint, constant)*rightPoint - polynom(degree, coefs, rightPoint, constant)*leftPoint)/(polynom(degree, coefs, leftPoint, constant)-polynom(degree, coefs, rightPoint, constant));


	}
}


void newtonRalphson(int degree, double constant, double coefs[], double error){
	double leftPoint;
	double rightPoint;
	
	double X0;
	double X1;
	bool isFound = false;
	printf("Baslangic noktasini giriniz: \n");
	scanf("%lf", &leftPoint);
	
	X0 = leftPoint;
	
	while(!isFound){
		X1 = X0 - (polynom(degree, coefs, X0, constant)/merkeziturev(degree, coefs, X0, constant, 0.001));
		if(X1 - X0 < error){
			isFound = true;
			printf("kok : %lf", X1);
		}
		else{
			X0 = X1;
		}
	}
}


void getCofactor(double mat[N][N], double temp[N][N], int p,
                 int q, int n)
{
    int i = 0, j = 0;
	int row,col;
    for (row = 0; row < n; row++)
    {
        for (col = 0; col < n; col++)
        {
            if (row != p && col != q)
            {
                temp[i][j++] = mat[row][col];
                if (j == n - 1)
                {
                    j = 0;
                    i++;
                }
            }
        }
    }
}
double determinantOfMatrix(double mat[N][N], int n)
{
    double D = 0; // Initialize result
 	double temp[N][N];
	int i;
    if (n == 1)
        return mat[0][0];
    int sign = 1;
    for (i = 0; i < n; i++)
    {
        getCofactor(mat, temp, 0, i, n);
        D += sign * mat[0][i] * determinantOfMatrix(temp, n - 1);
        sign = -sign;
    }
 
    return D;
}

void adjoint(double matris[N][N], double adj[N][N], int size){
	if(size == 1){
		adj[0][0] = 1;
		return;
	}
	int i ,j;
	int sign = 1;
	double temp[N][N];
    for (i=0; i<size; i++)
    {
        for (j=0; j<size; j++)
        {
            getCofactor(matris, temp, i, j, size);
            sign = ((i+j)%2==0) ? 1: -1;
            adj[j][i] = (sign)*(determinantOfMatrix(temp, size-1));
        }
    }
}


void inverse(double cof[N][N], int size, double determinant){
	if(determinant == 0.0){
		printf("TERSI YOK");
		return;
	}
	int i,j;
	for (i = 0; i < size; i++)
	{
		for(j = 0; j < size; j++){
			cof[i][j] /= determinant;
		}
	}
	
}


void nxninv(double matris[N][N],double cofactors[N][N], int unknowns){
	int i,j;
	adjoint(matris, cofactors, unknowns);
	double determinant = determinantOfMatrix(matris, unknowns);
	inverse(cofactors, unknowns, determinant);
	for (i = 0; i < unknowns; i++){
			for (j = 0; j < unknowns; j++)
			{
				printf("%lf ", cofactors[i][j]);
			}
			printf("\n");
	}
}

void dividerow(double row[], double num, int unknowns){
	int i;
	for(i = 0; i <= unknowns; i++){
		row[i] = row[i]/num;
	}
}

void addrow(double coef, double row1[], double row2[], int unknowns){
	
	int i;
	for(i = 0; i <= unknowns; i++) row2[i]-= coef*row1[i]; 
}

void printmatris(double matris[N][N], double vector[], int unknowns){
	int i,j;
	for (i = 0; i <= unknowns; i++){
			for (j = 0; j <= unknowns; j++)
			{
				printf("%lf ", matris[i][j]);
			}
			
	}
}


bool isApplicable(double matris[N][N], int size){
	int i,j;
	for(i = 0; i < size; i++){
		for(j = 0; j < size; j++){
			if(matris[i][i] < matris[i][j]){
				printf("bu matris icin gauss seidel yontemi uygulanamaz");
				return false;
			}
		}
	}
	return true;
}

void gaussElimination(double matris[N][N], int unknowns){
	double c;
	double vector[N];
	double sum;
	int i,j;
	
		for(i = 0; i < unknowns; i++){
		for(j = i+1; j < unknowns; j++){
			c = matris[j][i]/matris[i][i];
			addrow(c,matris[i], matris[j], unknowns);
		}
	}
	for(i = 0; i < unknowns; i++){
		dividerow(matris[i], matris[i][i], unknowns);
	}
	vector[unknowns-1] = matris[unknowns-1][unknowns];
	for(i = unknowns - 2; i >= 0; i--){
		sum = 0;
		for(j = i+1; j < unknowns; j++){
			sum += matris[i][j]*vector[j];
		}
		vector[i] = matris[i][unknowns] - sum;
	}
	for (i = 0; i < unknowns; i++){
		printf("%d. degisken: %lf\n", i+1, vector[i]);
	}
	
	
}

void gaussSiedel(double matris[N][N], int unknowns){
	int i,j, counter = 0;
	double error = 0.001;
	bool isFound = false; 
	double X0;
	double vector[N];
	if(isApplicable(matris, unknowns)){
		printf("Baslangic degerlerini giriniz: \n");
		for(i = 0; i < unknowns; i++){
			printf("%d. bilinmeyen: ", i+1);
			scanf("%d", &vector[i]);
			printf("\n");
		}
		while(!isFound){
			counter = 0;
			for(i = 0; i < unknowns; i++){
				X0 = vector[i];
				vector[i] = matris[i][unknowns];
				for(j = 0; j < unknowns; j++){
					if(i != j){
						vector[i] -= matris[i][j] * vector[j];
					}
				}
				vector[i] /= matris[i][i];
				if(absolute(vector[i] - X0) < error){
					//printf("za %lf\n", vector[i] - X0);
					counter++;
				} 
			}
			if(counter == unknowns) isFound = true;
		}
		for(i = 0; i < unknowns; i++) printf("%d. bilinmeyen: %lf\n",i+1 ,vector[i]);
		}
		else{
			printf("bu denklem takimina gauss siedal uygulanamaz.");
		}
	}


void trapez(int degree, double coefs[], double constant){
	double leftPoint, rightPoint, prev, h, X0, temp, area = 0;
	int i,n;
	printf("integralin sol tarafini giriniz: \n");
	scanf("%lf", &leftPoint);
	printf("integralin sag tarafini giriniz: \n");
	scanf("%lf", &rightPoint);
	printf("kac araliga bolmek istiyorsunuz?\n");
	scanf("%d", &n);
	prev = polynom(degree, coefs, leftPoint,constant);
	h = (rightPoint - leftPoint)/n;
	X0 = leftPoint + h;
	for(i = 0; i < n; i++){
		temp = polynom(degree, coefs, X0, constant);
		area += ((temp + prev)*h)/2;
		prev = temp;
		X0 += h;
	}

	printf("Girilen polinomun %lf ile %lf arasindaki belirli integral: %lf", leftPoint, rightPoint ,area);
}
void printPoly(int degree, double coefs[], double constant){
	int i;
	for(i = 0; i < degree; i++){
		printf("%lf*x^%d + ", coefs[i], degree-i);
	}
	printf("%lf\n", constant);
}
void simpson1over3(int degree, double coefs[], double constant){
	double leftPoint, rightPoint, h, area = 0;
	int i,n;
	printf("integralin sol tarafini giriniz: \n");
	scanf("%lf", &leftPoint);
	printf("integralin sag tarafini giriniz: \n");
	scanf("%lf", &rightPoint);
	printf("kac araliga bolmek istiyorsunuz?\n");
	scanf("%d", &n);

	h = (rightPoint - leftPoint)/n;
	area += polynom(degree, coefs, leftPoint,constant) + polynom(degree, coefs, rightPoint,constant);
	for(i = 1; i < n; i++){
		if(i % 2 == 1){
			area+= 4*polynom(degree, coefs, leftPoint + h*i,constant);
		}
		else{
			area+= 2*polynom(degree, coefs, leftPoint + h*i,constant);
		}
	}
	area *= h/3;
	printf("Girilen polinomun %lf ile %lf arasindaki belirli integral: %lf", leftPoint, rightPoint ,area);
}

void forwardDif(double values[], double h, double dif[N][N], int n){
	bool isZero = false;
	int i,j;
	for(i = 1; i < n; i++){
		dif[0][i-1] = values[i] - values[i-1];
	}
	n--;
	i = 1;
	while(!isZero){
		for(j = 1; j < n; j++){
			dif[i][j-1] = dif[i-1][j] - dif[i-1][j-1];
		}
		if(dif[i][0] == 0) isZero = true;
		i++;
	}
	i--;
	for(j = 0; j < i; j++){
		dif[i][j] = dif[i-1][1] - dif[i-1][0];
	}
	
}

void simpson3over8(int degree, double coefs[], double constant){
	double leftPoint, rightPoint, X0, h, temp, area = 0, print1, print2;
	int i,n;
	printf("integralin sol tarafini giriniz: \n");
	scanf("%lf", &leftPoint);
	printf("integralin sag tarafini giriniz: \n");
	scanf("%lf", &rightPoint);
	print1 = leftPoint;
	print2 = rightPoint;
	printf("kac araliga bolmek istiyorsunuz?\n");
	scanf("%d", &n);
	X0 = (rightPoint-leftPoint)/n;
	h = (X0)/3;
	for(i = 0; i < n; i++){
		temp = leftPoint + (X0);
		area += (polynom(degree, coefs, leftPoint,constant) + 3*polynom(degree, coefs, leftPoint+h,constant) + 3*polynom(degree, coefs, leftPoint+2*h,constant) + polynom(degree, coefs, leftPoint+3*h,constant))*((temp-leftPoint)/8);
		leftPoint = temp;
	}
	printf("Girilen polinomun %lf ile %lf arasindaki belirli integral: %lf", print1,print2,area);
}

void interpolation(double matris[N][N], double vector[]){
	int unknowns,i,j;
	double X0,h,x,temp, sum = 0;
	printf("kac deger var ?\n");
	scanf("%d", &unknowns);
	printf("baslangic degerini giriniz: \n");
	scanf("%lf", &X0);
	printf("aradaki farki giriniz: \n");
	scanf("%lf", &h);
	printf("x: \n");
	scanf("%lf", &x);
	for(i = 0; i < unknowns; i++){
		printf("f(%lf): ", X0+(i*h));
		scanf("%lf", &vector[i]);
	}
	forwardDif(vector, h, matris, unknowns);
	printf("Geri fark tablosu:\n");
	for(i = 0; i < unknowns; i++){
		for(j = 0; j < unknowns - i -1; j++){
			if(matris[i][j] != 0)
				printf("%lf ", matris[i][j]);
		}
		printf("\n");
	}
	sum += vector[0];
	for(i = 1; i < unknowns; i++)
	{
		temp = matris[i-1][0]/factorial(i);
		for(j = i; j > 0; j--){
			temp *= x-j+1;
		}
		sum += temp;
	}
	printf("%lf noktasindaki ongorulen deger: %lf", x,sum);
}



int main(){
	bool isRootThere = false;
	bool isLineer;
	bool isFound = false;
	char ch;
    int i,j,k,n;
	int counter, itN = 0;
    double y;
	double prev;
	double temp;
	double area = 0;
    double leftPoint;
    double rightPoint;
    double midPoint;
    double error;
    double x;
    double c;
	double X0;
	double h;
	double X1;
    double zero = 1;
    double deltaZero = 1;
    double constant, sum = 0;
	double vector[N];
	double matris[N][N];
	double cofactors[N][N];
    double coefs[200];
    int degree;
	int size = 0;
	int unknowns;

		
		
	/*else{
		
		
	}*/
	int input;

    do
    {
		printf("\n0. Cikis\n");
        printf( "1. Bisection yontemi\n" );
        printf( "2. Regula-Falsi yontemi\n" );
        printf( "3. Newton-Rapshon yontemi\n" );
        printf( "4. NxN'lik bir matrisin tersi\n" );
        printf( "5. Gauss Eleminasyon\n" );
        printf( "6. Gauss Seidal yontemleri\n" );
        printf( "7. Sayisal Turev\n" );
        printf( "8. Simpson yontemleri ile integral\n" );
        printf( "9. Trapez yontemi ile integral\n" );
        printf( "10. Degisken donusumsuz Gregory newton Enterpolasyonu\n" );
        printf( "Selection: " );
        scanf( "%d", &input );

        switch ( input ) 
        {
        case 1:
			printf("Kacinci dereceden bir polinom ?\n");
			scanf("%d", &degree);
			for(i = 0; i < degree; i++){
				printf("%d dereceli terimin katsayisini giriniz: ", degree-i);
				scanf("%lf", &coefs[i]);
			}
			
			printf("sabit degeri giriniz: \n");
			scanf("%lf", &constant);
			printPoly(degree, coefs, constant);
			printf("hatayi belirleyiniz: \n");
			scanf("%lf", &error);
            bisection(degree, coefs, constant, error);
            break;
        case 2:          
            printf("Kacinci dereceden bir polinom ?\n");
			scanf("%d", &degree);
			for(i = 0; i < degree; i++){
				printf("%d dereceli terimin katsayisini giriniz: ", degree-i);
				scanf("%lf", &coefs[i]);
			}
			printf("sabit degeri giriniz: \n");
			scanf("%lf", &constant);
			printPoly(degree, coefs, constant);
			printf("hatayi belirleyiniz: \n");
			scanf("%lf", &error);
			RegulaFalsi(degree, coefs, constant, error);
            break;
        case 3:         
            printf("Kacinci dereceden bir polinom ?\n");
			scanf("%d", &degree);
			for(i = 0; i < degree; i++){
				printf("%d dereceli terimin katsayisini giriniz: ", degree-i);
				scanf("%lf", &coefs[i]);
			}
			printf("sabit degeri giriniz: \n");
			scanf("%lf", &constant);
			printPoly(degree, coefs, constant);
			printf("hatayi belirleyiniz: \n");
			scanf("%lf", &error);
			newtonRalphson(degree, constant,coefs, error );
            break;
        case 4:        
            printf("Matris Boyutu?\n");
			scanf("%d" ,&unknowns);
			for (i = 0; i < unknowns; i++)
			{
				for (j = 0; j < unknowns; j++)
				{
					printf("tersi alinacak matrisin %d %d elemanini giriniz:\n", i, j);
					scanf("%lf", &matris[i][j]);
				}
				
			}
			printf("Matris: \n");
			for(i = 0; i < unknowns; i++){
				for(j = 0; j < unknowns; j++){
					printf("%lf ",matris[i][j]);
				}
				printf("\n");
			}
			printf("Matris Tersi: \n");
			nxninv(matris, cofactors, unknowns);
            break;
		case 5:
			printf("kac bilinmeyenli bir denklem?\n");
			scanf("%d" ,&unknowns);
			for (i = 0; i < unknowns; i++)
			{
				for (j = 0; j < unknowns; j++)
				{
					printf("matrisin %d %d elemanini giriniz:\n", i, j);
					scanf("%lf", &matris[i][j]);
				}
				printf("%d. denklemin sonucunu giriniz.\n", i+1);
				scanf("%lf", &matris[i][unknowns]);
			}
			printf("Matris: \n");
			for(i = 0; i < unknowns; i++){
				for(j = 0; j <= unknowns; j++){
					printf("%lf ",matris[i][j]);
				}
				printf("\n");
			}
			gaussElimination(matris, unknowns);
			break;	
		case 6:
			printf("kac bilinmeyenli bir denklem?\n");
			scanf("%d" ,&unknowns);
			for (i = 0; i < unknowns; i++)
			{
				for (j = 0; j < unknowns; j++)
				{
					printf("matrisin %d %d elemanini giriniz:\n", i, j);
					scanf("%lf", &matris[i][j]);
				}
				printf("%d. denklemin sonucunu giriniz.\n", i+1);
				scanf("%lf", &matris[i][unknowns]);
			}
			gaussSiedel(matris, unknowns);
			break;	
		case 7:
			printf("1. Geri Fark Yontemi\n");
			printf("2. Ileri Fark Yontemi\n");
			printf("3. Merkezi Fark Yontemi\n");
			printf("Istediginiz yontemi seciniz\n");
			scanf("%d", &k);
			printf("Kacinci dereceden bir polinom ?\n");
			scanf("%d", &degree);
			for(i = 0; i < degree; i++){
				printf("%d dereceli terimin katsayisini giriniz: ", degree-i);
				scanf("%lf", &coefs[i]);
			}
			printf("sabit degeri giriniz: \n");
			scanf("%lf", &constant);
			printPoly(degree, coefs, constant);
			printf("Hangi noktanin turevini bulmak istiyorsunuz?\n");
			scanf("%lf", &x);
			printf("Hata miktarini giriniz: \n");
			scanf("%lf", &h);
			if(k == 1) printf("Girilen polinomun %lf noktasindaki turevi: %lf\n", x,geriturev(degree, coefs, x, constant,h));
			else if(k == 2) printf("Girilen polinomun %lf noktasindaki turevi: %lf\n",x, ileriturev(degree, coefs, x, constant,h));
			else if(k == 3)	printf("Girilen polinomun %lf noktasindaki turevi: %lf\n",x,merkeziturev(degree, coefs, x, constant,h));
			else printf("Yanlis girdi girdiniz.\n");
			break;	
		case 8:
			printf("1. Simpson 1/3 yontemi\n");
			printf("2. Simpson 3/8 yontemi\n");
			printf("Istediginiz yontemi seciniz.\n");
			scanf("%d", &k);
			printf("Kacinci dereceden bir polinom ?\n");
			scanf("%d", &degree);
			for(i = 0; i < degree; i++){
				printf("%d dereceli terimin katsayisini giriniz: ", degree-i);
				scanf("%lf", &coefs[i]);
			}
			printf("sabit degeri giriniz: \n");
			scanf("%lf", &constant);
			printPoly(degree, coefs, constant);
			if(k == 1) simpson1over3(degree, coefs, constant);	
			else if(k == 2) simpson3over8(degree, coefs, constant);
			else printf("Yanlis girdi girdiniz.\n");
			break;	
		case 9:
			printf("Kacinci dereceden bir polinom ?\n");
			scanf("%d", &degree);
			for(i = 0; i < degree; i++){
				printf("%d dereceli terimin katsayisini giriniz: ", degree-i);
				scanf("%lf", &coefs[i]);
			}
			printf("sabit degeri giriniz: \n");
			scanf("%lf", &constant);
			printPoly(degree, coefs, constant);
			trapez(degree, coefs, constant);
			break;	
		case 10:
			interpolation(matris, vector);
			break;
		case 0: 
			isFound = true;
			break;
		default:            
            printf( "Yanlis girdi girdiniz.\n" );
            break;
        }
    }  while (!isFound);
	
    return 0;
    
}
