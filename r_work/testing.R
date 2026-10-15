CEBPGK562_Experiment1_Rep1=import("ExampleData/CEBPGK562_Experiment1_Rep1.bed.gz",format="narrowPeak")
CEBPGK562_Experiment1_Rep2=import("ExampleData/CEBPGK562_Experiment1_Rep2.bed.gz",format="narrowPeak")
CEBPGK562_Experiment2_Rep1=import("ExampleData/CEBPGK562_Experiment2_Rep1.bed.gz",format="narrowPeak")
CEBPGK562_Experiment2_Rep2=import("ExampleData/CEBPGK562_Experiment2_Rep2.bed.gz",format="narrowPeak")

Bedlist=list(CEBPGK562_Experiment1_Rep1,CEBPGK562_Experiment1_Rep2,CEBPGK562_Experiment2_Rep1,CEBPGK562_Experiment2_Rep2)
UnifiedSequence=Unify(Bedlist)
Chrom=ChromGet(UnifiedSequence)

BiasHMM=hmm.discnp::hmm(Chrom[["ChromSplit"]],K=3,itmax=1000)
BiasViterbi=hmm.discnp::viterbi(Chrom[["ChromSplit"]],BiasHMM)

ComputedScores=ScoreRegions(UnifiedSequence,Viterbi = BiasViterbi,Bedlist,Scorename = "signalValue")
# The signalValue here is not fold change from control, results will differ from paper

ContiguousRegions=GetRegions(Chrom[["ChromSplit"]],BiasViterbi,UnifiedSequence)