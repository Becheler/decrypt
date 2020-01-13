  seed = 123

  seqfile = temp_seq.phyl
  outfile = out.txt
  Imapfile = temp_imap.txt
  mcmcfile = mcmc.txt

* speciesdelimitation = 0          * fixed species tree
* speciesdelimitation = 1 0 2      * speciesdelimitation algorithm0 and finetune(e)
  speciesdelimitation = 1 1 2 1    * speciesdelimitation algorithm1 finetune (a m)
  speciestree = 1  0.4 0.2 0.1     * speciestree pSlider ExpandRatio ShrinkRatio

  speciesmodelprior = 0          * 0: uniform LH; 1:uniform rooted trees; 2: uniformSLH; 3: uniformSRooted
  species&tree = 2  Pop1  Pop2
                    30  30
                   (Pop1,Pop2);
  diploid = 0  0                 * 0: phased sequences, 1: unphased diploid sequences
  usedata = 1                    * 0: no data (prior); 1:seq like
  nloci = 10                      * number of data sets in seqfile

  cleandata = 1                  * remove sites with ambiguity data (1:yes, 0:no)?

  thetaprior = 3 0.8        * invgamma(a, b) for theta, theta = 4.N.mu
  tauprior =  3 1.2         * invgamma(a, b) for root tau & Dirichlet(a) for other tau's

  print = 1 0 0 0                * MCMC samples, locusrate, heredityscalars Genetrees
  burnin = 8000
  sampfreq = 2
  nsample = 1000
  finetune = 1: 0.2 0.002 0.01 0.001 0.1 0.1 0.1  * finetune for GBtj, GBspr, theta, tau, mix, locusrate, seqerr
