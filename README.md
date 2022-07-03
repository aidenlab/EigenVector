The main function is **bestEigen3.R**. It takes uo to 5 arguments.  
x - the square contacts matrix in Matrix (sparse) format; REQUIRED  
a1,a2 - two vectors - approximations to the first and second eigenvectors; do not use unless you know what you are doing  
tol - relative error allowed in the eigenvector (default 1.0e-6)  
maxiter - maximum number of iteretions (default 100)  

The second function is **eigenVectorRscript.R** - it is a wrapper for bestEigen built to run via Rscript so that it can be run in parallel on several chromosomes.  
Put eigenVectorRscript.R and bestEigen3.R in same folder and while in this folder a typical run is  
**Rscript --vanilla eigenVectorRscript.R [options] fin fout binsize**  

fin - input file as produced by juicer dump function (in sparse form); REQUIRED  
fout - file to store the eigenvector; REQUIRED
binsize - binsize (in base-pairs) used with juicer dump command; REQUIRED  

Run  
Rscript eigenVectorRscript.R --help
to see optional parameters. They are:  

-t,--tolerance - precision (error in the eigenvector) - default is 1.0e-6  
-m,--maxiter - maximum iterations - default is 100  
-s,--size - chromosome length (in basepairs) - used to determine the number of bins; if not supplied the number is determined based on the highest position encountered and may be slightly smaller than than it should be. **Hence it is recommended to use this option.**  
-v,--verbose - whether to output information to stdout (TRUE or FALSE) - default is FALSE

**eigFromHicRscript.R** is very similar to eigenVectorRscript.R but instead reading file produced by juicer tools dump it resads the data directly from the hic file. So its first argument is hic file. Its second (additional) argument is the chromosome. Nore that chr1 is sometimes encoded as 1 and in such a case it needs to be 1 when calling eigFromHicRscript.R. It also has two additional options:  
-n, --norm - which normalization to use; NONE for no normalization; other possibilities are VC, VC_SQRT, KR, SCALE, etc.  Default is NONE  
-o, --obs - whether to use observed/expected (oe) matrix (recommended) or observed. Use -o observed for observed, otherwise o/e is used  
Other arguments and options are like in eigenVectorRscript.R. So a typical usage is:  
**Rscript eigFromHicRscript.R [options] hicfile chr fout binsize**  
Note that you will need strawr R package (in addition to the other two). To install strawr on linux do:  
**sudo Rscript -e 'remotes::install_github("aidenlab/straw/R")'**  
alternative open an R session and do:  
**remotes::install_github("aidenlab/straw/R")**

The run produces two files: fout and fout.report. The first one holds the eigenvector. The second one contains run information. The last line lists the following quantities:  
lam1 - the approximation to the first eigenvalue  
lam2 - the approximation to the second eigenvalue  
er - L2 norm of C*v - lam1*v where v is the eigenvector (normalised to have L2 norm of 1) and C is the correlation matrix  
niter - number of iterations (always a multiple of 10)  

Please note that a bound on the relative error in v is er/(lam1 - lam2) (if lam1 and lam2 were exact). We stop when this becomes < tol or we reach maxiter iterations. Since each time we do 10 iterations we may exceed maxiter to the nearest multiple of 10.

REMARK: we start with a random vector so the results are not identical between cosequtive runs; to change this behaviour uncomment line 3 of bestEigen3.R (#       set.seed(12345) ). You can choose other seed (not necessarily 12345).  

PLEASE NOTE:
You need to have Matrix and optparse R packages installed.  
We do not compute the correlation matrix since we only need the result of it's multiplication with a vector. We manage to keep only one matrix in sparse Matrix format so the dimension of the matrix almost does not affect the performance. What matters is the number of non-zero elements in the contact matrix produced by juicer dump command.  
Do not attempt to use tolerance (tol) smaller than 1.0e-15 since it is comparable to the roundoff error and the process may never converge.

Author: Moshe Olshansky;  e-mail: moshe.olshansky@gmail.com
 
