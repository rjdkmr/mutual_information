/*
 * This file is part of mutual_information program
 *
 *              Author: Rajendra Kumar
 *
 *        Copyright (C) 2014  Rajendra Kumar
 *
 * WARNING: This program is not completely tested for the calculation accuracy.
 *        Use it at your own RISK.
 *
 * This program calculates mutual information by using the method developed
 * in the following article:
 * Alexander Kraskov, Harald Stögbauer, and Peter Grassberger (2004)
 * Estimating mutual information
 * Physical Review E 69, 066138.
 * URL: http://link.aps.org/doi/10.1103/PhysRevE.69.066138
 *
 *
 * This program also calculates generalized correlation by using the method
 * discussed in the following article:
 * Oliver F. Lange andHelmut Grubmüller (2005)
 * Generalized correlation for biomolecular dynamics
 * Proteins: Structure, Function, and Bioinformatics 62:1053-1061.
 * URL: http://onlinelibrary.wiley.com/doi/10.1002/prot.20784/full
 *
 * mututal_information is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License along
 * with mututal_information. If not, see <http://www.gnu.org/licenses/>.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
 * "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
 * LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
 * A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
 * OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
 * SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED
 * TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
 * PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
 * LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
 * NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
 * SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 *
 */

#include <unistd.h>
#include <iostream>
#include <cmath>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <pthread.h>
#include "digamma.hpp"
#include "MI_depend.hpp"
#include "Mutual_information.hpp"

using namespace std;
#define NTHREADS sysconf( _SC_NPROCESSORS_ONLN )

void *status;
pthread_t *thread;
pthread_attr_t attr;
int ithread;

MI_DATA miData;

void mi_init (double *x, double *y, int n, int k)	{
	miData.n = n;
	miData.k = k;
	miData.e = new double [n+1];
	miData.nx = new int [n+1];
	miData.ny = new int [n+1];
	miData.x = new double [n+1];
	miData.y = new double [n+1];
	for(int i=0;i<n;i++){
		miData.x[i] = x[i];
		miData.y[i] = y[i];
	}
	miData.prog = 0;
}

void mi_finished ()	{
	delete [] miData.e;
	delete [] miData.nx;
	delete [] miData.ny;
	delete [] miData.x;
	delete [] miData.y;
}

void *mi_calculate (void *arg) {
	double ex,ey;
	long t = long(arg);
	int start = ceil(t*((miData.n/NTHREADS)));
	int end = ceil((t+1)*(miData.n/NTHREADS)+NTHREADS);
	int i = start;

	if(NTHREADS==1){
		start = 0;
		end = miData.n;
	}
	if(end>miData.n)
		end=miData.n;
	//cout<<miData.n<<"\t"<<t<<"\t"<<i<<"\t"<<start<<"\t"<<end<<"\n";

	for(i=start;i<end;i++)	{
		double *d, *arr;
		int *s;
		d = new double [miData.n+1];
		s = new int [miData.n+1];
		arr = new double [miData.n+1];

		for(int j=0; j<miData.n; j++)
			arr[j] = miData.x[j];

		for(int j=0;j<miData.n;j++)
			s[j] = j;
		for (int j=0;j<miData.n;j++) {
			if(i==j)	{
				d[j]= 0;
				arr[j]=0;
			}
			else	{
				d[j] = dist(miData.x[i],miData.y[i],miData.x[j],miData.y[j]);
				arr[j] = d[j];
			}
		}

		//k-nearest neighbor searching
		if((i+1)%100==0)	{
			int flush_output = (float(miData.prog)/float(miData.n) )* 100;
				cout<<"# Progress:  "<<flush_output<<" % \r"<<flush;
		}
		miData.e[i] =0;
		nbSearch(arr,s,miData.n,miData.k);
		miData.e[i] = d[s[miData.k]]*2;
		//d[s[k]];
		//for (j=0;j<=miData.k;j++)
			//cout<<arr[j]<<"\t"<<d[s[j]]<<"\n";
		//cout<<t<<"\t"<<miData.e[i]<<"\n";

		//To fullfill Maximum Norm
		ex = fabs(miData.x[i]-miData.x[s[miData.k]]);
		ey = fabs(miData.y[i]-miData.y[s[miData.k]]);
		//cout<<ex<<"\t"<<ey<<"\n";
		miData.nx[i]=0; miData.ny[i]=0;
		if(ex>=ey)	{
			NumPoint(miData.x,&miData.x[i],&ex,&miData.nx[i],miData.n);
			NumPoint(miData.y,&miData.y[i],&ex,&miData.ny[i],miData.n);
		}
		else	{
			NumPoint(miData.x,&miData.x[i],&ey,&miData.nx[i],miData.n);
			NumPoint(miData.y,&miData.y[i],&ey,&miData.ny[i],miData.n);
		}
		miData.prog++;

		delete [] d;
		delete [] s;
		delete [] arr;
	}

//	for(i=start;i<end;i++)	{
//		cout<<t<<"\t"<<i<<"\t"<<miData.nx[i]<<"\t"<<miData.ny[i]<<"\t"<<miData.e[i]<<"\n";
//	}
	if(NTHREADS>1)
		pthread_exit(NULL);
	return (void *)0;
}

double mi (double *x, double *y, int n, int k)
{
	//cout<<"\n\n\n\nhelllooo\n=========================================\n";
	//cout<<n<<"testing\n";
	long i;
	double dx,dy, dxMin, dxMax, dyMin, dyMax;
	double psiN, psi1, psiK, *psiNx, *psiNy, HXY, HX, HY, mi;
	psiNx = new double [n+1];
	psiNy = new double [n+1];

	//Minimum and Maximum of X and Y and dimension of X-data and Y-data
	min_max(x,&dxMin,&dxMax,n);
	dx = fabs(dxMax-dxMin);
	min_max(y,&dyMin,&dyMax,n);
	dy = fabs(dyMax-dyMin);

	//Initialization of MI data structure
	mi_init (x, y, n, k);

	/* Initialize and set thread detached attribute */
	if(NTHREADS>1)	{
		thread = new pthread_t [NTHREADS];
		pthread_attr_init(&attr);
		pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);

		for(i=0;i<NTHREADS;i++)
			ithread = pthread_create(&thread[i],&attr, mi_calculate, (void *) i);
			//mi_calculate((void *)&i);

		 //Free attribute and wait for the other threads
		pthread_attr_destroy(&attr);
		for(i=0; i<NTHREADS; i++) {
			ithread = pthread_join(thread[i], &status);
			if (ithread) {
				printf("ERROR; return code from pthread_join() is %d\n", ithread);
				exit(-1);
			 }
		 }
	}
	for(i=0;i<n;i++){
		if(NTHREADS==1)
			mi_calculate((void *) i);
		//cout<<i<<"\t"<<miData.nx[i]<<"\t"<<miData.ny[i]<<"\t"<<miData.e[i]<<"\n";
		psiNx[i] = digama(miData.nx[i]+1);
		psiNy[i] = digama(miData.ny[i]+1);
	}

	psiN = digama(double(miData.n));
	psi1 = digama(double(1.000));
	psiK = digama(double(miData.k));
	HXY = psiN - psiK + log(1) + (((dx+dy)/n)*sum(miData.e,n));
	HX = (-sum(psiNx,n)/n) + psiN + log(1) + ((dx/n)*sum(miData.e,n));
	HY = (-sum(psiNy,n)/n) + psiN + log(1) + ((dy/n)*sum(miData.e,n));
	//cout<<"H(X,Y) = "<<HXY<<",\t H(X) = "<<HX<<",\t H(Y) = "<<HY<<",\t MI = "<<HX+HY-HXY<<",\t Norm MI = "<<(HX+HY-HXY)/((HX+HY)/2)<<"\n";
	mi = (psiK + psiN) - ((sum(psiNx,n) + sum(psiNy,n))/n);

	mi_finished();
	return mi;
}
