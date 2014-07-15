DNA_sequence_alignment
======================

This is a simple version of DNA sequence alignment written in C++, which is the procedure of comparing two (pair-wise alignment) or more (multiple alignment) DNA sequences by searching for a series of characters that are in the same order in all sequences. Two sequences can be aligned by writing them across a page in two rows. Identical or similar characters are placed in the same column, and non identical ones can either be placed in the same column as a mismatch or against a gap (-) in the other sequence. The program is written by using dynamic programming method, with the following recurrence relation:

	   D[i,j] = max{ D[i-1,j-1] + pair(s[j],t[i]),
        	         D[i-1,j] + deletion,
                	 D[i,j-1] + deletion }

    In other words, if you have an optimal alignment up to D[i-1,j-1] there are only three possibilities of what could happen next:

	1.	the characters for s[i] and t[j] match (if not match, the pair() function would be a negative score for subtituting a character to match the other sequence)
	2.	a gap (in our program, deletion) is introduced in t and
	3.	a gap is introduced in s

By implementing the recurrence relation above, aligning DNA sequences and finding the best matching DNAs can be achieved.


== Instruction ==
* type "make" in the terminal to activate 'Makefile' which compiles the program
* use following command line prompts to run the program: (for -m -c -d, using double type, ex: 2.0 instead 2)

   -1 followed by a filename that stores the 1st DNA sequence
   -2 followed by a filename that stores the 2nd DNA sequence  
   -m the similarity measure for matching two identical bases
   -c the similarity measure for changing a base in one of the sequences to match with a base in the other
   -d the similarity measure for skipping a base in one of the sequences and not matching it

Examples:
	./prog -1 seq1.txt -2 seq2.txt
	./prog -1 seq1.txt -2 seq2.txt -c -1 -d -10	
	./prog -1 fragment1.txt -2 fragment2.txt -m 10 -c -5 -d -2
