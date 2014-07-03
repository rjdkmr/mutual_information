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

#ifndef MUTUAL_INFORMATION_HPP_
#define MUTUAL_INFORMATION_HPP_

typedef struct {
	double *x, *y, *e;
	int *nx, *ny;
	int n, k;
	int prog;
} MI_DATA;

void mi_init (double *x, double *y, int n, int k);
void mi_finished ();
void *mi_calculate (void *arg);
double mi (double *x, double *y, int n, int k);



#endif /* MUTUAL_INFORMATION_HPP_ */
