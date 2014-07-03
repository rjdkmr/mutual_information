mutual_information
==================

This program calculates mutual information by using the method developed
in the following article:                                              
Alexander Kraskov, Harald Stögbauer, and Peter Grassberger (2004)       
Estimating mutual information                                           
Physical Review E 69, 066138.                                           
URL: http://link.aps.org/doi/10.1103/PhysRevE.69.066138                 
                                                                        
                                                                        
This program also calculates generalized correlation by using the method
discussed in the following article:                                     
Oliver F. Lange andHelmut Grubmüller (2005)                             
Generalized correlation for biomolecular dynamics                       
Proteins: Structure, Function, and Bioinformatics 62:1053-1061.         
URL: http://onlinelibrary.wiley.com/doi/10.1002/prot.20784/full         

WARNING: This program is not completely tested for the calculation accuracy.
         Use it at your own RISK.
         
###Download
<pre><code>git clone https://github.com/rjdkmr/mutual_information
</code></pre>
***

###Installation
<pre><code>cd mutual_information
mkdir build
cd build
cmake .. -DCMAKE_INSTALL_PREFIX=/opt/mutual_information
make
make install
</code></pre>
***

###Usage
<pre><code> mi [k] [nc] [INPUT FILE] [OUTPUT FILE MI] [OUTPUT FILE R]
</code></pre>

Where "k" is number of nearest neighbor points.
"nc" is no. of columns in input file.
e.g. if total no. of column is 5 and nc is 3 then 3x2 matrix will be computed. If nc=0, 5x5 matrix will be computed.

***
