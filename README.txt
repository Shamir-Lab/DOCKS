How to generate decycling set

java -jar decycling.jar <output file> <k> <Alphabet>

Example run:

java -jar decycling.jar decycling_5_ACGT.txt 5 ACGT

This command will genearte a decycling set for k=5 over ACGT alphabet. The output will be saved to decycling_5_ACGT.txt.

How to find additional k-mers to hit all L-long sequences

java -jar DOCKS.jar <output file> <input decycling file> <k> <L - sequence length> <Alphabet> <ILP time limit> <0/1/2/3 - random/greedyL/greedyAny/none> <x - optional with any>
Example DOCKS run:

java -jar DOCKS.jar res_5_20_ACGT_0_1.txt decycling_5_ACGT.txt 5 20 ACGT 0 1

Example DOCKSAny run:

java -jar DOCKS.jar res_5_20_ACGT_0_2.txt decycling_5_ACGT.txt 5 20 ACGT 0 2

Example DOCKSAnyX (X=125) run:

java -jar DOCKS.jar res_5_20_ACGT_0_2_125.txt decycling_5_ACGT.txt 5 20 ACGT 0 2 125

In all of the runs, the use needs to provide the decycling set computed by decycling.jar.

Example ILP run (limited to 1000s), with no DOCKS solution:

java -jar DOCKS.jar res_5_20_ACGT_1000_3.txt decycling_5_ACGT.txt 5 20 ACGT 1000 3

You may need to use more memory for higher values of k, i.e. k>=10, by adding -Xmx4096m option for example.
You may need to increase heapspace size for random mode for DFS, by adding -Xss515m option for example.

Interpreting the output

decycling's and DOCKS's outputs should look like this:
AAAACAAA
AAAAGAAA
AAAATAAA
AAACCAAA
...
ATAACGAA
TCACCGAA
GCCTACTA
TCCTCCTA
Each line is a k-mer in the set.
