Source code used for 

Hickey, G and Blanchette M.  "A Probabilistic Model For Sequence Alignment with Context-Sensitive Indels" RECOMB 2011.  (and submitted to JCB)

Released under the GNU General Public License (see gpl.txt for more information).

This is source code for a project that is still under development.  The interface, in particular, is presently being reimplemented so as to be somewhat usable.  If you are interested in running the present version in the meantime, some notes are provided below.  Also, please feel free to contact me (Glenn Hickey: hickey@mcb.mcgill.ca) and I will be happy to help. 

INSTALLATION

configure --prefix="installation dir"; make; make install

Note: The source code is huge because of the generated code to solve the integrals required by the context F84 model.   There shouldn't be any dependencies though.

QUICK START

Some parameter sets that were pre-computed using the training program are found in context/example/params.  For each species (described in the paper), four parameter sets are given: 

_sf84_djc : Felsenstien 84 substitution model with Jukes Cantor context dependence model
_sf84: Felsenstein 84 subsitution model with *no* context dependence
_sjc_djc : Jukes Cantor substitution model with Jukes Cantor context dependence model
_sjc : Jukes Cantor substitution model with *no* context dependence

Suppose we have a file in MFA-like format called human_dog.mfa that we would like to align (or realign) using the Felsenstein 84 / Jukes Cantor parameters for human/dog and a context window of length 2. 

human1  ACCACAAATTTA
dog1 ACCACAT--TTTA
human2  TCCC
dog2 CCAC

etc.

We would run

src/align infile=human_dog.mfa outfile=output.mfa fpfile=example/params/dog _sf84_djc.txt win=2

and the output would be in output.mfa 

Test data can be simulated using the "genpair" tool (todo: describe).  It can also be sampled from a given alignment.  Some such sample alignments (indeed, those used for the paper and sample parameters) are given in context/data.  In order to sample 100 pairwise alignments from the human/dog alignment, each of length [50 - 100], one can run

src/axtsample infile=data/chr22.hg19.canFam2.net.axt outfile=sample.mfa minlen=50 maxlen=100 n=100

This command creates the sampled file in sample.mfa, and displays some summary statistics on the screen.  These are estimated parameters (computed from the input file using simple counting statistics, not the training) of the entire input file, and of the sample.  They are mostly useful to give an idea of how representative the sample is of the input. 

The obtained sub-alignment sample can then be realigned as above:

src/align infile=sample.mfa outfile=sample_realigned fpfile=example/params/dog _sf84_djc.txt win=7


ALIGNMENT

align infile=<sequence pairs to align> outfile=<output alignment> fpfile=<model paramers> win=<context window size> outfile2=<posterior probability annotation (optional)>

HOMOLGY 

homolg infile=<sequence pairs to test> outfile=<total likelihood of each pair> fpfile=<model paramers> win=<context window size>

TRAINING

train infile=<sequence pairs for training> fpfile=<seed model parameters (optional)> fpout=<output model parameters> win=<context window size>

Advanced optional training options:

optits=<number of gradient descent iterations>
optruns=<number of gradient descent seeds tried per iteration>
optthreshold=<gradient descent threshold>
offset=<gradient descent offset>

emits=<EM loop erations>
emthreshold=<EM loops threshold>
ecr=<number of times likelihood can be within ecrthreshold before breaking>
ecrthreshold=<see above>
seed=<random seed>

sym=<true or false> : insertion and deletion forced to be the same

FILE FORMATS

Sequence pairs accepted for now in following format which describes a list of sequence pairs. 

name1 sequence 
name2 sequence
name1 sequence 
name2 sequence

Parameters are given in a text file in the followinf format

example:

<*t= 1 *mu= 0 a= 0.017441 b= 0.00895469 p0= 0.262795 p1= 0.241077 p2= 0.237599 p3= 0.258528 *ga= 0 *c= 0 d= 0.0531192 *q0= 0.25 *q1= 0.25 *q2= 0.25 *q3= 0.25 RD= 0.00017009 RI= 0.00017009 RMD= 0.000533192 RMI= 0.000533192 PE= 0.016656 PCD= 0.358262 PCI= 0.358262 KD= 0.747603 KI= 0.747603 >

F84 substitution model:
a + b = transition rate
a = transversion rate
p0, p1, p2, p3 = stationary distribution of A, C, G, T, respectively.

F84 context model:
c + d = transition rate
c = transversion rate
q0, q1, q2, q3 = stationary distribution of A, C, G, T, respectively.

RD: deletion rate
RI: insertion rate
RMD= context deletion rate
RMI= context insertion rate
1-PE= geometric parameter for sequence length
1-PCD=geometric parameter for context deletion
1-PCI=geometric parameter for context insertion
KD=geometric parameter for deletion
KI=geometric parameter for insertion

* before any identifier means that parameter won't be modified by training program (when passed as seed).  

 








  
