#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <time.h>

/** define some convergence threshold. This is a user defined parameter which can be changed. */
#define CONVERGE_THRESHOLD 0.000001

//#define M 4058
//#define N 1400


//double array_my[161464]={0};
struct mtx
{
    int row;
    int col;
	int fre;
};
struct mtx bbcnews[161462];
//struct mtx bbcnews[161462];
//double mtx[161462][3]={0};
//double mtx2[4058][1400]={1};

//void a(double *array);


void MatrixMultiply(double* a,double* b,double* c, int m, int n, int p);
double Minus(double* m1, double* m2, int row, int col);
void RandMtx(double* Arr, int row, int col);
void MatrixTranspose(double *a, double *b, int row, int col);
void bubble_sort(double a[],int b[],int n);


int run_nmf( int k )
{
	printf("** RUNNING NMF FOR k=%d ...\n", k );
    char tmp[512];
    char terms[4058][15];
    int i = 0;
	int j = 0;
	long int l = 0;
//	long int k = 0;
	int no_rows = 4058;
	int no_cols = 1400;
	long int N = 4058*1400;
	//double (*mtx2)[1400]=malloc(sizeof(double)*4058*1400);
	FILE *fp1,*fp2;
	
	int loop=0;
	int c = 0;

	double sum[1400]={0};
	double df[4058]={0};
	double w0[4058]={0};
	double w1[4058]={0};
	int p[4058]={0};
	double v[4058]={0};
//	double w[4058][k];
//	double wt[k][4058];
//	double h[k][1400];
//	double ht[1400][k];

	double *mtx2 = (double *) malloc(sizeof(double)*N);
	double *a = (double *) malloc(sizeof(double)*N);
	double *w = (double *) malloc(sizeof(double)*4058*k);
	double *wt = (double *) malloc(sizeof(double)*k*4058);
	double *h = (double *) malloc(sizeof(double)*k*1400);
	double *ht = (double *) malloc(sizeof(double)*1400*k);
	double *wta = (double *) malloc(sizeof(double)*k*1400);
	double *a0 = (double *) malloc(sizeof(double)*N);
	double *wta0 = (double *) malloc(sizeof(double)*k*1400);
	double *aht = (double *) malloc(sizeof(double)*4058*k);
	double *a0ht = (double *) malloc(sizeof(double)*4058*k);
	memset(mtx2, 0, sizeof(double)*N);
	memset(a, 0, sizeof(double)*N);
	memset(w, 0, sizeof(double)*4058*k);
	memset(wt, 0, sizeof(double)*k*4058);
	memset(h, 0, sizeof(double)*k*1400);
	memset(ht, 0, sizeof(double)*1400*k);
	memset(wta, 0, sizeof(double)*k*1400);
	memset(a0, 0, sizeof(double)*N);
	memset(wta0, 0, sizeof(double)*k*1400);
	memset(aht, 0, sizeof(double)*4058*k);
	memset(a0ht, 0, sizeof(double)*4058*k);

	double *a1 = (double *) malloc(sizeof(double)*N);
	memset(a1, 0, sizeof(double)*N);

	if (NULL == mtx2) {
		printf("memory allocation failed\n");
		return -1;
	}

	if (NULL == a) {
		printf("memory allocation failed\n");
		return -1;
	}
	
	if (NULL == w) {
		printf("memory allocation failed\n");
		return -1;
	}

	if (NULL == wt) {
		printf("memory allocation failed\n");
		return -1;
	}

	if (NULL == h) {
		printf("memory allocation failed\n");
		return -1;
	}

	if (NULL == ht) {
		printf("memory allocation failed\n");
		return -1;
	}

	if (NULL == wta) {
		printf("memory allocation failed\n");
		return -1;
	}

	if (NULL == a0) {
		printf("memory allocation failed\n");
		return -1;
	}

	if (NULL == wta0) {
		printf("memory allocation failed\n");
		return -1;
	}

	if (NULL == aht) {
		printf("memory allocation failed\n");
		return -1;
	}
	
	if (NULL == a0ht) {
		printf("memory allocation failed\n");
		return -1;
	}

	printf("Value of N: %ld\n", N);



	//double *mtx2=(double *)malloc(sizeof(double)*4058*1400);
	//long double a[]={0};
	//double b[161464];
    //struct mtx bbcnews[8000]; 
	//double array_my[5681200]={0};

//	int (*p)[1400];
//	p = ( int (*)[1400])malloc(4058);


    /* read into bbcnews.terms */
	fp1=fopen("bbcnews.terms","r");
	if(fp1==NULL)
    {
        printf("\nFailed to open the file");
        exit(0);
    }
    while( fgets( tmp,sizeof(tmp),fp1 )!=NULL ) 
    {
        sscanf( tmp,"%s",terms[i] ); 
//        printf("%s\n",terms[i]);
        i++;
//		printf("%s\n",tmp);
    }


	/* read into bbcnews.mtx */
    fp2=fopen("bbcnews.mtx","r");
    if(fp2==NULL)
    {
        printf("\nFailed to open the file");
        exit(0);
    }
    
	for(i=0;i<2;i++){
    	fgets( tmp, sizeof(tmp), fp2 );
	//	fgets( tmp, sizeof(tmp), fp2 );
    }
	i=0;

    while( fgets( tmp,sizeof(tmp),fp2 )!=NULL ) /*¨º?¨¨?¨ºy?Y*/
    {
        
		sscanf( tmp,"%d %d %d",&(bbcnews[i].row),&(bbcnews[i].col),&(bbcnews[i].fre)); 
//		printf("%s\n",tmp);
		//sscanf( tmp,"%d %d %d",&(mtx[i][0] ),&(mtx[i][1]),&(mtx[i][2])); 
	        //for(j=0;j<161664;j++) {
			//fscanf(fp,"%1d",array_my[j]);
		    //printf("Number is: %d\n", array_my[j]);
	//}
//    	printf("%d %d %d\n",bbcnews[i].row,bbcnews[i].col,bbcnews[i].fre);
        i++;   
		//printf("%d\n",i);
    }
	
	
// put the data into new matrix 
	for(l=0; l<161462; l++)
	{
		//mtx2[bbcnews[i].row][bbcnews[i].col]=bbcnews[i].fre;
		i = bbcnews[l].row - 1;		
		j = bbcnews[l].col - 1;
		mtx2[i*no_cols+j] = bbcnews[l].fre;
	}
		
//    for(i=0;i<4058;i++)
//		for(j=0;j<1400;j++)
//			if(j!=(bbcnews[i].col-1))
//			mtx2[i][j];
//	printf("Printing new matrix...\n");
	for(i=0;i<4058;i++)
	{
		for(j=0;j<1400;j++)
		{
			if (0 != mtx2[i*no_cols+j]) {
//				printf("i: %d, j: %d, value: %f\n", i, j, mtx2[i*no_cols+j]);
			//	printf("%d\n",mtx2[i*no_cols+j]);
			}
			//printf("%d ",mtx2[i*no_cols+j]);
		}
	}
 


// Apply TF-IDF normalization to the term-document matrix.
//void a(double *array)
//{
//	int i,j;
//	long int N = 4058*1400;

//	double *tf = (double *) malloc(sizeof(double)*N);
//	double *idf = (double *) malloc(sizeof(double)*N);
    
	for(j=0;j<1400;j++) {
		sum[j]=0;
		for(i=0;i<4058;i++){
			//sum[j]=0;
			sum[j] += mtx2[i*no_cols+j];
		}
		//printf("%d\n",sum[j]);
	}
	
	for(i=0;i<4058;i++)	{
		df[i]=0;
		for(j=0;j<1400;j++) {
		
			if(mtx2[i*no_cols+j] !=0 )
			df[i]++;
           
		}
	//	printf("%d\n",df[i]);
	}
	
	
	//for(i=0;i<4058;i++){
	//	if(df[i]!=0)
	//		printf("%d\n",df[i]);
	//}
	

	/*
	for(j=0;j<1400;j++)
		for(i=0;i<4058;i++)
			tf[i*no_cols+j]=mtx[i*no_cols+j]/sum[j];
	for(j=0;j<1400;j++)
		for(i=0;i<4058;i++)
			idf[i*no_cols+j]=log(1400/df[i]);
	for(j=0;j<1400;j++)
		for(i=0;i<4058;i++)
			a[i*no_cols+j]=tf[i][j]*idf[i][j];
	*/



	for(i=0;i<4058;i++){
		for(j=0;j<1400;j++){
			a[i*no_cols+j] = (mtx2[i*no_cols+j]/sum[j])*log(1400/df[i]+1);
			//if (mtx2[i*no_cols+j] != 0)
			//{

				//printf("i: %d, j: %d, mtx2: %.08f\n", i, j, mtx2[i*no_cols+j]);
			//   }
//			if (0 != mtx2[i*no_cols+j]) 
//			printf("i: %d, j: %d, a: %.08f, mtx2: %.08f \n", i, j, a[i*no_cols+j], mtx2[i*no_cols+j]);			
			//printf("i: %d, j: %d, value: %.08f\n", i, j, a[i][j]);
			//printf("%f \n",a[i*no_cols+j]);
		}
	}

		RandMtx(w,4058,k);
		RandMtx(h,k,1400);
	
/* Update factor H, W */
	for(loop=0;loop<100;loop++)
	{
/* Update factor H */
		MatrixMultiply(w,h,a0,4058,1400,k);
		MatrixTranspose(w,wt,4058,k);
		MatrixMultiply(wt,a,wta,k,1400,4058);
		MatrixMultiply(wt,a0,wta0,k,1400,4058);
		for(c=0;c<k;c++)
		{
			for(j=0;j<1400;j++)
			{
				h[c*1400+j] = h[c*1400+j]*(wta[c*1400+j]/wta0[c*1400+j]);
			}
		}
//		MatrixTranspose(h,ht,k,1400);

/* Update factor W */
		MatrixTranspose(h,ht,k,1400);
		MatrixMultiply(w,h,a0,4058,1400,k);
		MatrixMultiply(a,ht,aht,4058,k,1400);
		MatrixMultiply(a0,ht,a0ht,4058,k,1400);
		for(i=0;i<4058;i++)
		{
			for(c=0;c<k;c++)
			{
				w[i*k+c] = w[i*k+c]*(aht[i*k+c]/a0ht[i*k+c]);
			}
		}
//		MatrixMultiply(w,h,a0,4058,1400,k);
//		MatrixTranspose(w,wt,4058,k);

		/** Test convergence, by computing error of approximation to the original matrix A */
		int pos;
		double error = 0.0;
		MatrixMultiply(w,h,a1,4058,1400,k);
		for( pos = 0; pos < N; pos++)
		{
			error += ( (a1[pos]-a[pos]) * (a1[pos]-a[pos]) );
		}
		/* normalize by matrix size */
		error /= N;
		printf( "Iteration %d: Error = %.7f\n", (loop+1), error );
		/* have we reached convergence ? */
		if( error < CONVERGE_THRESHOLD )
		{
			printf("Algorithm converged after %d iterations.", (loop+1) );
			break;
		}
	}

/* Print the matrix w */
	for(i=0;i<4058;i++)
	{
		for(j=0;j<k;j++)
		{
			if (0 != w[i*k+j]) {
//			printf("i: %d, j: %d, w: %.08f \n", i, j, w[i*k+j]);
			}
		}
	}

/* Print the matrix h */
	for(i=0;i<k;i++)
	{
		for(j=0;j<1400;j++)
		{
			if (0 != h[i*1400+j]) {
//			printf("i: %d, j: %d, h: %.08f \n", i, j, h[i*1400+j]);
			}
		}
	}

/*	for(i=0;i<4058;i++)
	{
		for(j=0;j<1400;j++)
		{
			if (0 != a0[i*1400+j]) {
			printf("i: %d, j: %d, a0: %.08f \n", i, j, a0[i*1400+j]);
			}
		}
	}
*/


//	for(i=0;i<4058;i++)
//	{
//		w0[i]=w[i*k+0];
//		w1[i]=w[i*k+1];
//		p[i]=i;
//		printf("%d\n",p[i]);
//	}
	
	/* sort the array */ 
	for(j=0;j<k;j++) 
	{
		for(i=0;i<4058;i++) 
		{
			v[i]=w[i*k+j];
			p[i]=i;
		}                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                            
		bubble_sort(v,p,4058);
		printf("CLUSTER %d\n", (j+1) );
		for(i=0;i<10;i++)
		{
			printf("%s\n",terms[p[i]]);
		}
		printf("\n");
	}


//	for(j=0;j<k;j++)
//	bubble_sort(v[j],p,4058);
/*	for(i=0;i<10;i++)
	{
		printf("%s\n",terms[p[i]]);
	}
	
//	bubble_sort(w1,p,4058);
	for(i=0;i<10;i++)
	{
		printf("%s\n",terms[p[i]]);
	}
*/

/*	for(i=0;i<10;i++)
	{
		printf("%lf\n",w0[i]);
	}

	for(i=0;i<10;i++)
	{
		printf("%lf\n",w1[i]);
	} */
	
	free(mtx2);	
	free(a);
	free(w);
	free(wt);
	free(h);
	free(ht);     
	free(wta);
	free(a0);
	free(wta0);
	free(aht);
	fclose(fp1);
	free(a0ht);
	fclose(fp2);

	free(a1);

 	return 0;
}



/* multiplication of matrix */
void MatrixMultiply(double* a,double* b,double* c, int m, int n, int p)
{
  int i,j,k;
    for(i = 0; i < m; i++)
    {
		for(j = 0; j < n; j++) 
		{
			double sum = 0;
			for(k = 0; k < p; k++) {
                sum += a[i*p+k] * b[k*n+j];
//				c[i*n+j] += a[i*p+k] * b[k*n+j];
			}
			c[i*n+j] = sum;
        }
    }
}


/* Euclidean distance between two matrixes */
double Minus(double* m1, double* m2, int row, int col) {
		int i,j;	
 		double value = 0;
		for ( i = 0; i < row; i ++ )
		{
			for ( j = 0; j < col; j ++ ) 
			{
				value += pow((m1[i*col+j] - m2[i*col+j]), 2);
			}
		}
	return sqrt(value);
}


/* Randomly initialize factors for NMF. ---> W and H */               
void RandMtx(double* Arr, int row, int col) {
	int i,j;                                    
   	srand(time(NULL));                                                                
    for (i = 0; i < row; i++)            
    {                                  
        for (j = 0; j < col; j++)        
        {                              
            Arr[i*col+j] = (double)rand()/(double)RAND_MAX;
        }                              
    }
}                                  
                                                      
/* Update factor H, W */
/* transpose(w) */
/*double wt[k][4058];
for(i=0;i<k;i++)
{
	for(j=0;j<4058;j++)
	{
		wt[i][j]=w[j][i];
	}
}

double h.transpose()
double ht[1400][k];
for(i=0;i<1400;i++)
{
	for(j=0;j<k;j++)
	{
		ht[i][j]=h[j][i];
	}
}
*/

void MatrixTranspose(double *a, double *b, int row, int col)
{
	int i,j;
 	for(i = 0;i < row; i++)
 	{
 		for(j = 0; j < col; j++)
 		{
 			b[j*row+i] = a[i*col+j];
 		}
 	} 
}



/* rank words and rank the index of words at sametime */
void bubble_sort(double a[],int b[],int n)//n?a¨ºy¡Á¨¦a¦Ì??a????¨ºy
{
    int i,j,temp2;
    double temp1;
    for(j=0;j<n-1;j++)
        for(i=0;i<n-1-j;i++)
        {
            if(a[i]<a[i+1])//¨ºy¡Á¨¦?a??¡ä¨®D?¡ã¡ä¨¦yD¨°??¨¢D
            {
                temp1=a[i];
                a[i]=a[i+1];
                a[i+1]=temp1;
				temp2=b[i];
                b[i]=b[i+1];
                b[i+1]=temp2;
            }                             
        }
}


int main()
{	
	int k = 0;
	for(k=2;k<=6;k++)
	{
		run_nmf( k );
	}	
	return 0;
}
