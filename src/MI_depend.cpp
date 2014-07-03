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

#include <iostream>
#include <cmath>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "digamma.hpp"
#include "MI_depend.hpp"
using namespace std;

double dist(double x1, double y1, double x2, double y2){
	double x =0, y=0;
	x = (x1-x2)*(x1-x2);
	y = (y1-y2)*(y1-y2);
	return sqrt(x+y);
}

void min_max (double *x, double *min, double *max, int n)	{
	int i=0;
	*min = x[0];
	*max = x[0];
	for(i=0;i<n;i++)	{
		if(x[i]>*max)
			*max = x[i];
		if(x[i]<*min)
			*min = x[i];
	}
}

void NumPoint (double *a, double *q, double *r, int *num, int n)	{
	int i =0;
	double temp[1],left,right;
	temp[0] = *q + *r;	temp[1] = *q - *r;
	min_max(temp,&left,&right,2);
	for(i=0;i<n;i++)	{
		if((a[i]>left) && (a[i]<right))
			*num = *num + 1;
	}
}

void nbSearch (double *arr, int *s, int n, int k)	{
	int i = 0, j = 0;
	int tmp_ind, min_ind;
	double tmp_arr, min_arr;
	for(i=0;i<n;i++)	{
		min_arr=arr[i];
		min_ind=s[i];
		for(j=0;j<=k;j++)	{
			if(arr[j]>min_arr)		{
				min_arr = arr[j];
				tmp_arr = arr[i];
				arr[i] = arr[j];
				arr[j] = tmp_arr;

				min_ind = s[j];
				tmp_ind = s[i];
				s[i] = s[j];
				s[j] = tmp_ind;
			}
		}
	}
}

double sum(double *a, int n)	{
	int i=0;
	double value=0;
	for(i=0;i<n;i++)
		value = value + a[i];
	return value;
}
