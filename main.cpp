
//***** Freddy Yang **********************
//***** DNA Sequence Alignment ***********
//***** March 2014 ***********************

#include<iostream>
#include<string>
#include<iomanip>
#include <sstream>
#include <fstream>
#include <stdio.h>
#include <stdlib.h>
using namespace std;

const int k=5000; //maximum length of the comparing DNA sequences

double db[k][k]; 
double matching(2), changing(-0.5), deletion(-1);

typedef struct 
{
   string bp;
   int size;
} sequence;

// measure constants when matching and subtituting
double pairs(sequence &seq1,sequence &seq2,int i,int j)
{   
	double d = 0;
	if(seq1.bp[j-1]==seq2.bp[i-1])
		d = matching; // the similarity measure for matching two identical bases
	else
		d = changing; // the similarity measure for changing a base in one of the sequences to match with a base in the other
	return d;
}

// find max among a,b,c
double domax(double a,double b,double c)
{
	if( a>=b && a>=c )
		return a;
	else if(b>=a&&b>=c)
		return b;
	else
		return c;
}

// compare sequences
void alig(sequence &seq1,sequence &seq2,double db[][k])  //sequence alignment
{  
	double a=0,b=0,c=0;
	for(int i=0;i<=seq1.size;i++)
		db[0][i]=0-i*2;
	for(int j=0;j<=seq2.size;j++)
		db[j][0]=0-j*2;
	for(int m=1;m<=seq2.size;m++)
	{
		for(int n=1;n<=seq1.size;n++)
		{
		  // Recurrence Relation
			a=db[m][n-1] + deletion;   // the similarity measure for skipping a base in one of the sequences and not matching it
			b=db[m-1][n] + deletion;   // same as above
			c=db[m-1][n-1]+ pairs(seq1,seq2,m,n);
		db[m][n]=domax(a,b,c);
		}
	}

}

/*   ****************TESTING*****************
void oput(sequence seq1,sequence seq2,double db[][k])
{
	cout<<"alig result:"<<endl;
    cout<<setw(7)<<" ";
    for(int i=0;i<seq1.size;i++)
     cout<<setw(6)<<seq1.bp[i];
     cout<<endl;
	 cout<<endl;
   for(int m=0; m<=seq2.size;m++)
   {     
	   if(m>=1)
			   cout<<seq2.bp[m-1]<<setw(7);
       if(m==0)
			  cout<<setw(1)<<" ";
	   for(int n=0;n<=seq1.size;n++)
	   {
		  cout<<setw(6)<<db[m][n];
		  
	   }
	   cout<<"\n\n";
	   }
}
*/

// Read a .txt file
string get_element(char* filename)
{
        string line,outp;
	ifstream myfile(filename);
	while (getline(myfile, line))
	{
		outp = line;
	}
	myfile.close();
	return outp;
}

int main(int argc, char* argv[])
{
	sequence seq1,seq2;
	string seq1_name,seq2_name;

	// Command Line Prompt
	for (int i = 0; i < argc; i++)
	{
		if (std::string(argv[i]) == "-1")
		{
			seq1.bp = get_element(argv[i+1]);
			seq1.size=seq1.bp.length();
			seq1_name = argv[i+1];
		}
		else if (std::string(argv[i]) == "-2")
		{
			seq2.bp = get_element(argv[i+1]);
			seq2.size=seq2.bp.length();
			seq2_name = argv[i+1];
		}
		else if (std::string(argv[i]) == "-m")
			matching = atof(argv[i+1]);
		else if (std::string(argv[i]) == "-c")
			changing = atof(argv[i+1]);
		else if (std::string(argv[i]) == "-d")
		        deletion = atof(argv[i+1]);
	}
	
	alig(seq1,seq2,db);

	
	// OUTPUT(stdout)
	for (int i = 0; i < seq1_name.length()-4; i++)
	  cout << seq1_name[i];
	cout << " ";
	
	for (int j = 0; j < seq2_name.length()-4; j++)
	  cout << seq2_name[j];
	cout << " ";

        cout << db[seq2.size][seq1.size] << endl;
	return 0;
}
