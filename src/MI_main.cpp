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
#include <fstream>
#include <cmath>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <iomanip>
#include <sstream>
#include <vector>
#include "digamma.hpp"
#include "MI_depend.hpp"

using namespace std;

void CopyRightMsg() {

    const char *copyright[] = {
            "                                                                        ",
            "             :-)  mututal_information (-:                               ",
            "                                                                        ",
            "               Author: Rajendra Kumar                                   ",
            "                                                                        ",
            "         Copyright (C) 2014  Rajendra Kumar                             ",
            "                                                                        ",
            "WARNING: This program is not completely tested for the calculation accuracy.",
            "         Use it at your own RISK.                                       ",
            "                                                                        ",
            "This program calculates mutual information by using the method developed",
            "in the following article:                                              ",
            "Alexander Kraskov, Harald Stögbauer, and Peter Grassberger (2004)       ",
            "Estimating mutual information                                           ",
            "Physical Review E 69, 066138.                                           ",
            "URL: http://link.aps.org/doi/10.1103/PhysRevE.69.066138                 ",
            "                                                                        ",
            "                                                                        ",
            "This program also calculates generalized correlation by using the method",
            "discussed in the following article:                                     ",
            "Oliver F. Lange andHelmut Grubmüller (2005)                             ",
            "Generalized correlation for biomolecular dynamics                       ",
            "Proteins: Structure, Function, and Bioinformatics 62:1053-1061.         ",
            "URL: http://onlinelibrary.wiley.com/doi/10.1002/prot.20784/full         ",
            "                                                                        ",
            "mututal_information is distributed in the hope that it will be useful,  ",
            "but WITHOUT ANY WARRANTY; without even the implied warranty of          ",
            "MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the           ",
            "GNU General Public License for more details.                            ",
            "                                                                        ",
            "You should have received a copy of the GNU General Public License along ",
            "with mututal_information. If not, see <http://www.gnu.org/licenses/>.   ",
            "                                                                        ",
            "THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS     ",
            "\"AS IS\" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT     ",
            "LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR   ",
            "A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT    ",
            "OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,   ",
            "SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED",
            "TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR  ",
            "PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF  ",
            "LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING    ",
            "NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS      ",
            "SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.            ",
            "                                                                        ",
            "                        :-)  mututal_information (-:                    ",
            "                                                                        ",
            "                                                                        "
    };
    int i = 0;
    char *str;
    for(i=0; i<47; i++) {
        str = strdup(copyright[i]);
        fprintf(stderr,"%s\n", str);
    }
}


int main (int argc, char *argv[]){

	double temp;
	char str[20000];
	int i=0,j=0;
	int row=0, col=0;
	vector< vector <double> > x;

	if (argc != 6)	{
		fprintf(stderr,"\n===================== HELP ========================\n");
		fprintf(stderr,"Number of argument should be 5 ...\n");
		fprintf(stderr,"Usage: ./mutual_information <k> <nc> <INPUT FILE> <OUTPUT FILE MI> <OUTPUT FILE R>\n");
		fprintf(stderr,"Where \"k\" is number of nearest neighbor points\n");
		fprintf(stderr,"Where \"nc\" is no. of column. e.g. if total no. of column is 5 and nc is 3 then 3x2 matrix will be computed. If nc=0, 5x5 matrix will be computed.\n");
		fprintf(stderr,"\n====================================================\n");
		exit(1);
	}

	CopyRightMsg();

	int k = atoi(argv[1]);
	int vs = atoi(argv[2]);
	//opening of input file stream
	ifstream fin;
	fin.open(argv[3], ios::in|ios::ate);

	//opening of input file stream
	ofstream foutMinf, foutR;
	foutMinf.open(argv[4], ios::out);
	foutR.open(argv[5], ios::out);

	//Determining the total size of the data
	size_t size;
	size = fin.tellg();
	x.reserve(size);

	//Determining the number of row
	fin.seekg(0,ios::beg);
	fin.getline(str,20000);
	istringstream is(str);
	while(is)	{
		is>>temp;
		col++;
	}
	col = col-1;

	//Reading the data from input file
	fin.seekg(0,ios::beg);
	while(fin.getline(str,20000))	{
		istringstream is(str);
		for(i=0;i<col;i++){
			is>>temp;
			x[i].push_back(temp);
			//cout<<x[i][col]<<"\t";
		}
		//cout<<"\n";
		row++;
	}

	//Calculation of Mutual Information and Generalized Correlation Coefficiant
	//vector< vector <double> > minf (col,vector<double> (col,0));
	//vector< vector <double> > r  (col,vector<double> (col,0));
	double **minf=NULL, **r=NULL;

	if(vs==0)	{
		minf = new double * [col];
		r = new double * [col];
		for(i=0;i<col;i++)	{
			minf[i] = new double [col];
			r[i] = new double [col];
		}
		cout.precision(5);

		for(i=0;i<col;i++)	{
			for(j=0;j<=i;j++)	{
				if(i==j){
					minf[i][j] = 8.0;
					minf[j][i] = 8.0;
					r[i][j] = 1.0;
					r[i][j] = 1.0;
					continue;
				}
				cout<<"\n# Column: "<<i+1<<"\tvs\tColumn: "<<j+1<<"\n"<<flush;
				//Mututal Information
				minf[i][j] = mi(&x[i][0],&x[j][0],row,k);
				minf[j][i] = minf[i][j];
				//Generalized correlation coefficiant
				if(minf[i][j]>0)
					r[i][j] = sqrt(1 - exp(-(2*minf[i][j])/3));
				else
					r[i][j] = 0;
				r[j][i] = r[i][j];
				cout.clear();
				cout<<"\n MI = "<<minf[i][j]<<"\n r = "<<r[i][j]<<"\n"<<flush;
			}
		}
	}
	else	{
		minf = new double * [vs];
		r = new double * [vs];
		cout.precision(5);
		for(i=0;i<vs;i++)	{
			minf[i] = new double [col-vs];
			r[i] = new double [col-vs];
			for(j=0;j<col-vs;j++)	{
				cout<<"\n# Column: "<<i+1<<"\tvs\tColumn: "<<j+vs+1<<"\n"<<flush;
				//Mututal Information
				minf[i][j] = mi(&x[i][0],&x[j+vs][0],row,k);

				//Generalized correlation coefficiant
				if(minf[i][j]>0)
					r[i][j] = sqrt(1 - exp(-(2*minf[i][j])/3));
				else
					r[i][j] = 0;
				cout.clear();
				cout<<"\n MI = "<<minf[i][j]<<"\n r = "<<r[i][j]<<"\n"<<flush;
			}
		}
	}

	if(vs==0)	{
		foutMinf.precision(5);
		for(i=0;i<col;i++)	{
			for(j=0;j<col;j++)	{
				foutMinf<<fixed<<minf[i][j]<<"\t\t";
			}
			foutMinf<<"\n";
		}
		foutR.precision(5);
		for(i=0;i<col;i++)	{
			for(j=0;j<col;j++)	{
				foutR<<fixed<<r[i][j]<<"\t\t";
			}
			foutR<<"\n";
		}
	}
	else	{
		foutMinf.precision(5);
		for(i=0;i<col-vs;i++)	{
			for(j=0;j<vs;j++)	{
				foutMinf<<fixed<<minf[j][i]<<"\t\t";
			}
			foutMinf<<"\n";
		}
		foutR.precision(5);
		for(i=0;i<col-vs;i++)	{
			for(j=0;j<vs;j++)	{
				foutR<<fixed<<r[j][i]<<"\t\t";
			}
			foutR<<"\n";
		}
	}

	fprintf(stdout, "\nThanks for using mutual_information!!!\n");
	fprintf(stdout, "\n++++ PLEASE READ AND CITE THE FOLLOWING REFERENCE ++++\n");
	fprintf(stderr, "-------- -------- ------------------- -------- ----------\n");
	fprintf(stderr, "Alexander Kraskov, Harald Stögbauer, and Peter Grassberger (2004)\n");
	fprintf(stderr, "Estimating mutual information A\n");
	fprintf(stderr, "Physical Review E 69, 066138.\n");
	fprintf(stderr, "-------- -------- ------------------- -------- ----------\n");
	fprintf(stderr, "Oliver F. Lange andHelmut Grubmüller (2005)\n");
	fprintf(stderr, "Generalized correlation for biomolecular dynamics \n");
	fprintf(stderr, "Proteins: Structure, Function, and Bioinformatics 62:1053-1061.\n");
	fprintf(stderr, "-------- -------- ------------------- -------- ----------\n\n");

	return EXIT_SUCCESS;
}
