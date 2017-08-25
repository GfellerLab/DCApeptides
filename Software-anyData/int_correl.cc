#include <iostream>
#include <stdlib.h>
#include <fstream>
#include <math.h>
#include <algorithm>
#include <vector>
#include <string>
#include <sstream>
#include <cstdio>
#include <cstring>
#include <sys/stat.h>

//Compile with:
// g++ -O2 -o int_correl.x int_correl.cc -llapack -lblas -lm -w

using namespace std;
extern "C" {
    // LU decomoposition of a symmetric matrix
    void dsptrf_(char *UPLO, int *N, double* AP, int* IPIV, int* INFO);

    // generate inverse of a matrix given its LU decomposition
    void dsptri_(char *uplo, int *n, double *a, int *ipivot, double* WORK, int *info);
    
}
int getArguments (int argc, char **argv);

void import_domain_seq();
void invert();
void direct_coupling();
int position(char s);
void import_ligands(int pref, int lref, int tLProteinDomain);
void import_ligands_sep(int pref, int lref, int tLProteinDomain);


void weights();
void init_param();
void compute_freq();

void mutual_information();
void mutual_information_pval();
void mutual_information_rc();


void remove_columns(int **tseq,int tLProteinDomain);
//void remove_columns_ref(int **tseq, int pref, int tM, int tLProteinDomain);
int readFastatAlignmentFile();
int readFastatAlignmentFileLigand();
int printAlignment ();
int printAlignmentLigand ();
void freeSequenceList ();
void freeSequenceListLigand ();
//Global variables

int Q;   //Number of amino acid types (20 or 21)
int tQ;
int tLProteinDomain;  //Number of positions before removing some columns in the alignment
int tLpeptide;  //Number of positions before removing some columns in the alignment

int tM;  //Number of sequences before filtering
int tMp; //Number of sequences in protein before filtering
//int tMl; //Number of sequences in ligand before filtering
int lengthM;  //Number of lines in ligand alignment file
int lengthMpeptide; //Number of lines in ligand alignment file
int tref;//temporary Position of the ref sequence in the initial alignment
int pref;
int lref;

int L; //Number of positions  
int LProteinDomain; //Number of positions in the domain alignment
int Lpeptide; //Number of positions in the ligand alignment
int M;  //Number of final protein (domain) sequences
int Mligand;  //Number of final protein (domain) sequences
int DNA; 


int comp_weight;///if we compute weight

int **seq;  //Sequences (nummerical)
int **tseq;  //Initial sequences
int **tseq2;
int ***lseq;

int *str_pos;
int *seq_pos;  //Keep track of the initial position of the sequences (before removing columns) 

char **tname;
char **name;
//char ***lname;


double **fr;   //amino acid frequencies at each position
double **correl;  //amino acid pair frequencies (N*21 x N*21)
double **inv;  //inverse of the correlation matrix
double *w;   //weights
double **MI;  //mutual information
double **DI;  //Direct information
double RC;
char *letter;


int str_ct;
int lstr_ct;

int *nl;


char    ref_name[100000], family_dir[100000], family_name[100000],family_name_Ligand[100000],weight_name[4096];
struct sequence
{
    char            *label, *sequence;
    struct sequence *next;
};
struct sequence seqList, seqListLigand;
#define MAX_LINE_LEN 10000


//#include "mutual_info.h"

int main(int argc, char **argv){
    int  status;
    if (getArguments (argc, argv) == -1){
        status = 1;
        cout<<"ERROR IN READING ARGUMENTS\n";
    }
    cout<<"DEBUG1: read protein\n";
    readFastatAlignmentFile();
    tLProteinDomain = strlen (seqList.sequence);
    
    cout<<"DEBUG2b: import protein domain\n";
    import_domain_seq();

    cout<<"DEBUG2a: read ligand\n";
    readFastatAlignmentFileLigand();
    tLpeptide=strlen (seqListLigand.sequence);
    lref=0;
    
    cout<<"DEBUG2b: import_ligands\n";
    import_ligands(pref, lref, LProteinDomain);
     
    cout<<"DEBUG3: weights\n";
    weights();
 
    cout<<"DEBUG4: Compute frequency\n";
    compute_freq();
    invert();

    cout<<"DEBUG5: Direct coupling analysis\n";
    direct_coupling();
     
    cout<<"DEBUG5: Free data\n";
    freeSequenceList();
    freeSequenceListLigand();
    
    for(int m=0; m<M; m++){
        delete[] name[m];
    }
    delete[] name;
    for(int l1=0; l1<L; l1++){
        delete[] DI[l1];;
    }
    delete[] DI;
    for(int l=0; l<L; l++){
        for(int q=0; q<tQ; q++){
            delete[] correl[l*tQ+q];
        }
    }
    delete[]correl;
    for(int i=0; i<tM; i++){
        delete[] tseq[i];
        delete[] tname[i];
    }
    for(int m=0; m<tM; m++){
        delete[] tseq2[m];
    }
    delete[]tseq2;
    delete[]tseq;
    delete[]tname;
    delete[]letter;
    delete[]str_pos;
    for(int m=0; m<M; m++){
        delete[] seq[m];
        
    }
    delete[] seq;
    delete[] w;
    for (int l=0; l<L; l++){
        delete[] fr[l];
    }
    delete[] fr;
    for(int i=0; i<L*Q; i++){
        delete[] inv[i];
    }
    delete[] inv;
    delete[] nl;
    
}
void freeSequenceList (){
    struct sequence *seq, *next;
    
    /*
     ** Free up all list elements and their contents.
     */
    seq = seqList.next;
    free (seqList.label);
    free (seqList.sequence);
    while (seq != NULL)
    {
        next = seq->next;
        free (seq->label);
        free (seq->sequence);
        free (seq);
        seq = next;
    }
}

void freeSequenceListLigand (){
    struct sequence *seq, *next;
    
    /*
     ** Free up all list elements and their contents.
     */
    seq = seqListLigand.next;
    free (seqListLigand.label);
    free (seqListLigand.sequence);
    while (seq != NULL)
    {
        next = seq->next;
        free (seq->label);
        free (seq->sequence);
        free (seq);
        seq = next;
    }
}
int readFastatAlignmentFile(){
    // cout<<"In readFastatAlignmentFile\n";
    int status=0;
    char             line[MAX_LINE_LEN];
    int              num_chars, len;
    struct sequence *seq, *last;
    FILE            *fp;
    /*
     ** Read the sequences. FASTA file format is assumed.
     */
    /*
     ** Open the alignment file for reading.
     */
    fp = fopen (family_name, "r");
    if (fp == NULL)
    {
        num_chars = -1;
        fprintf (stderr, "Could not open the alignment file '%s'.\n", family_name);
        status=-1;
        goto End_of_Routine;
    }
    
    seq = NULL;
    while (fgets (line, MAX_LINE_LEN, fp) != NULL)
    {
        /*
         ** Remove the trailing newline character, if present.
         */
        lengthM++;
        len = strlen (line);
        if (line[len-1] == '\n')
        {
            line[len-1] = '\0';
        }
        /*
         ** If a new sequence is encountered, create a new list element and save
         ** the sequence label.
         */
        if (line[0] == '>')
        {
            tM++;
            last = seq;
            if (seq == NULL)
            {
                seq = &seqList;
            }
            else
            {
                seq = (struct sequence*)malloc (sizeof (struct sequence));
                last->next = seq;
            }
            seq->label = strdup (&(line[1]));
            seq->sequence = NULL;
            seq->next = NULL;
        }
        /*
         ** Otherwise, save the actual sequence.
         */
        else
        {
            if (seq->sequence == NULL)
            {
                seq->sequence = strdup (line);
            }
            else
            {
                len = strlen (seq->sequence) + strlen (line) + 1;
                seq->sequence = (char *)realloc (seq->sequence, len*sizeof (char));
                strcat (seq->sequence, line);
            }
        }
    }
    fclose (fp);
    
    tLProteinDomain = strlen (seqList.sequence);
   // cout<<"In readFastatAlignmentFile: tLProteinDomain="<<tLProteinDomain<<endl;
    End_of_Routine:
    /*
     ** Close the file and return the length of the Newick tree representation.
     */
    return (status);
}
int readFastatAlignmentFileLigand(){
    int status=0;
    char             line[MAX_LINE_LEN];
    int              num_chars, len;
    struct sequence *seq, *last;
    FILE            *fp;
    /*
     ** Read the sequences. FASTA file format is assumed.
     */
    /*
     ** Open the alignment file for reading.
     */
    fp = fopen (family_name_Ligand, "r");
    if (fp == NULL)
    {
        num_chars = -1;
        fprintf (stderr, "Could not open the alignment file '%s'.\n", family_name_Ligand);
        status=-1;
        goto End_of_Routine;
    }
    
    seq = NULL;
    while (fgets (line, MAX_LINE_LEN, fp) != NULL)
    {
        /*
         ** Remove the trailing newline character, if present.
         */
        lengthMpeptide++;
        len = strlen (line);
        if (line[len-1] == '\n')
        {
            line[len-1] = '\0';
        }
        /*
         ** If a new sequence is encountered, create a new list element and save
         ** the sequence label.
         */
        if (line[0] == '>')
        {
            Mligand++;
            last = seq;
            if (seq == NULL)
            {
                seq = &seqListLigand;
            }
            else
            {
                seq = (struct sequence*)malloc (sizeof (struct sequence));
                last->next = seq;
            }
            seq->label = strdup (&(line[1]));
            seq->sequence = NULL;
            seq->next = NULL;
        }
        /*
         ** Otherwise, save the actual sequence.
         */
        else
        {
            if (seq->sequence == NULL)
            {
                seq->sequence = strdup (line);
            }
            else
            {
                len = strlen (seq->sequence) + strlen (line) + 1;
                seq->sequence = (char *)realloc (seq->sequence, len*sizeof (char));
                strcat (seq->sequence, line);
            }
        }
    }
    fclose (fp);
    
    tLpeptide = strlen (seqListLigand.sequence);
    End_of_Routine:
    /*
     ** Close the file and return the length of the Newick tree representation.
     */
    return (status);
}

int getArguments (int argc, char **argv){
    int status, i;
    cout<<"DEBUG0: getArguments\n";
    strcpy (family_name,"MHC_seq.txt");  //protein File from PFAM
    strcpy (family_name_Ligand,"MHC_ligand.txt");  //peptides
    strcpy (family_dir,"MHC_bs/");              //Directory output
    strcpy (ref_name,"C8CJC1/1-178");  //Reference sequence used to count the residue and map to the strucutre
    strcpy (weight_name,"weight.txt");//str_ct=1;  //Start position on the structure: la premiere position  (indiec 0) de l alignment de la sequence de reference  correspond a la position 1 dans structure pdb
    comp_weight=1;//yes we compute weight
    

    DNA=0;  //0: protein, 1: DNA sequence
    if(DNA==0){
        Q=21;  //This is the total number of letters
    }
    tLProteinDomain=0;
    tLpeptide=0;
    tM=0;
    lengthM=0;
    lengthMpeptide=0;
    tQ=Q-1;
    status = 0;
    i = 1;
    
    
    

    
    while (i < argc)
    {
        if (strcmp (argv[i], "-ref") == 0)
        {
            strcpy (ref_name, argv[++i]);
            i++;
            // cout<<"-ref\n";
        }
        else if (strcmp (argv[i], "-alignProtein") == 0)
        {
            strcpy (family_name, argv[++i]);
            i++;
            // cout<<"-align\n";
        }
        
        else if (strcmp (argv[i], "-alignPeptide") == 0)
        {
            strcpy (family_name_Ligand, argv[++i]);
            i++;
            // cout<<"-align\n";
        }
        else if (strcmp (argv[i], "-out") == 0)
        {
            strcpy (family_dir, argv[++i]);
            i++;
            // cout<<"-out\n";
        }
        else if (strcmp (argv[i], "-weight") == 0)
        {
            strcpy (weight_name, argv[++i]);
            comp_weight=0;
            i++;
            // cout<<"-out\n";
        }
        else
        {
            status = -1;
            fprintf (stderr, "Unknown argument '%s'.\n", argv[i]);
            goto End_of_Routine;
        }
    }
    
End_of_Routine:
    /*
     ** Return the status.
     */
    return (status);
}

void direct_coupling(){
    
    int D=2*Q;
    double **pr;
    pr=new double*[Q];
    for(int q=0; q<Q; q++){
        pr[q]=new double[Q];
    }
    
    //DI values
    DI=new double*[L];
    for(int l1=0; l1<L; l1++){
        DI[l1]=new double[L];
    }
    
    double *eh;
    eh=new double[D]; //exp(-h)
    double *neh;
    neh=new double[D];
    
    double diff;
    double thresh=0.00005;
    double tot;
    double max;
    
    double **mat; //inverse of the correlation matrix for (i,j)
    mat=new double*[Q];
    for(int q1=0; q1<Q; q1++){
        mat[q1]=new double[Q];
    }
    int count;
    //L=2;
    int ct;
    double *rtot;
    rtot=new double[Q];
    FILE *F;
    char file [4096];
    
    for(int l1=0; l1<L; l1++){
        for(int l2=l1+1; l2<L; l2++){
            diff=1;
            
            for(int q1=0; q1<Q; q1++){
                for(int q2=0; q2<Q; q2++){
                    mat[q1][q2]=inv[l1*Q+q1][l2*Q+q2];
                }
            }
            
            for(int q=0; q<D; q++){
                eh[q]=1.0/Q;
            }
            ct=0;
            while(diff>thresh && ct<5000){
                ct++;
                
                tot=0;
                for(int q1=0; q1<Q; q1++){
                    tot=0;
                    for(int q2=0; q2<Q; q2++){
                        tot=tot+mat[q1][q2]*eh[q2+Q];
                    }
                    neh[q1]=fr[l1][q1]/(tot);
                }
                
                //Normalize
                tot=0;
                for(int q1=0; q1<Q; q1++){
                    tot=tot+neh[q1];
                }
                for(int q1=0; q1<Q; q1++){
                    neh[q1]=neh[q1]/tot;
                }
                
                
                for(int q1=0; q1<Q; q1++){
                    tot=0;
                    for(int q2=0; q2<Q; q2++){
                        tot=tot+mat[q2][q1]*eh[q2];
                    }
                    neh[q1+Q]=fr[l2][q1]/(tot);
                }
                
                //Normalize
                tot=0;
                for(int q1=0; q1<Q; q1++){
                    tot=tot+neh[q1+Q];
                }
                for(int q1=0; q1<Q; q1++){
                    neh[q1+Q]=neh[q1+Q]/tot;
                }
                
                
                diff=0;
                for(int q1=0; q1<D; q1++){
                    if(diff < sqrt((eh[q1]-neh[q1])*(eh[q1]-neh[q1]))){
                        diff = sqrt((eh[q1]-neh[q1])*(eh[q1]-neh[q1]));
                    }
                    eh[q1]=neh[q1];
                }
            }
            if(diff>thresh){
                cout<<"No covergence after 500 steps (pos "<<l1<<" and "<<l2<<")\n";
            }
            
            tot=0;
            for(int q1=0; q1<Q; q1++){
                for(int q2=0; q2<Q; q2++){
                    pr[q1][q2]=mat[q1][q2]*eh[q1]*eh[q2+Q];
                    // pr[q1][q2]=mat[q1][q2];
                    tot=tot+pr[q1][q2];
                }
            }
            for(int q1=0; q1<Q; q1++){
                for(int q2=0; q2<Q; q2++){
                    pr[q1][q2]=pr[q1][q2]/tot;
                }
            }
            
            //Check the marginal
            //cout<<"Validation\n";
            for(int q1=0; q1<Q; q1++){
                tot=0;
                for(int q2=0; q2<Q; q2++){
                    tot=tot+pr[q1][q2];
                }
                if(tot - fr[l1][q1] > 0.001 || tot - fr[l1][q1] < -0.001 ){
                    cout<<"Bad marginal for "<<l1<<" and "<<l2<<": "<<tot<<"\t"<<fr[l1][q1]<<endl;
                }
            }
            for(int q1=0; q1<Q; q1++){
                tot=0;
                for(int q2=0; q2<Q; q2++){
                    tot=tot+pr[q2][q1];
                }
                if(tot - fr[l2][q1] > 0.001 || tot - fr[l2][q1] < -0.001 ){
                    cout<<"Bad marginal for "<<l1<<" and "<<l2<<": "<<tot<<"\t"<<fr[l2][q1]<<endl;
                }
            }
            
            DI[l1][l2]=0;
            for(int q1=0; q1<Q; q1++){
                for(int q2=0; q2<Q; q2++){
                    if(pr[q1][q2]>0){
                        DI[l1][l2]=DI[l1][l2]+pr[q1][q2]*log(pr[q1][q2]/(fr[l1][q1]*fr[l2][q2]));
                    }
                }
            }
            if(l1 < LProteinDomain && l2 >= LProteinDomain){
                sprintf(file, "%s/EC_matrices/%d_%d.txt", family_dir, seq_pos[str_pos[l1]], seq_pos[str_pos[l2]]);
                
                F=fopen(file, "w");
                fprintf(F, "%.5f\t", DI[l1][l2]);
                for(int q1=0; q1<Q; q1++){
                    fprintf(F, "%c\t", letter[q1]);
                }
                fprintf(F, "\n");
                for(int q1=0; q1<Q; q1++){
                    fprintf(F, "%c\t", letter[q1]);
                    for(int q2=0; q2<Q; q2++){
                        if(pr[q1][q2]>0){
                            fprintf(F, "%.5f\t", pr[q1][q2]*log(pr[q1][q2]/(fr[l1][q1]*fr[l2][q2])));
                        } else{
                            fprintf(F, "0");
                        }
                    }
                    fprintf(F, "\n");
                }
                fclose(F);
            }
            
            DI[l2][l1]=DI[l1][l2];
        }
    }
    
    int pos1, pos2;
    char let1, let2;
    char chain1, chain2;
    
    sprintf(file, "%s/EC_list.txt", family_dir);
    F=fopen(file, "w");
    for(int l1=0; l1<L; l1++){
        for(int l2=l1+1; l2<L; l2++){
                if(l1<LProteinDomain){
                    let1=letter[seq[pref][l1]];
                    pos1=seq_pos[str_pos[l1]];
                    chain1='A';
                }
                if(l1>=LProteinDomain){
                    let1=letter[lseq[pref][lref][l1-LProteinDomain]];
                    pos1=seq_pos[str_pos[l1]];
                    chain1='B';
                }
                if(l2<LProteinDomain){
                    let2=letter[seq[pref][l2]];
                    pos2=seq_pos[str_pos[l2]];
                    chain2='A';
                }
                if(l2>=LProteinDomain){
                    let2=letter[lseq[pref][lref][l2-LProteinDomain]];
                    pos2=seq_pos[str_pos[l2]];
                    chain2='B';
                }
                if(chain1==chain2){
                if (pos1 >pos2){
                cout<<"ERROR!!! in chain "<<chain1<<" where "<<pos1 <<"> "<<pos2<<endl;
                }
                }
                fprintf(F, "%c\t%c\t%c%d\t%c%d\t%.5f\n",  let1, let2, chain1, pos1, chain2, pos2,  DI[l1][l2]);
            } 
    }
    fclose(F);
    delete[] rtot;
    delete[] neh;
    delete[] eh;
    for(int q1=0; q1<Q; q1++){
        delete[]  mat[q1];
    }
    delete[] mat;
    for(int l1=0; l1<Q; l1++){
        delete[]  pr[l1];
    }
    delete[] pr;
    
}


//Inverse the correlation matrix
void invert(){
    
    int D=L*tQ;
    double* a;
    a=new double[D*(D+1)/2];
    
    int ct=0;
    for(int j=0; j<D; j++){
        for(int i=0; i<=j; i++){
            a[ct]=correl[i][j];
            ct++;
        }
    }
    
    char uplo;
    uplo='U';
    int *ipiv = new int[D];
    int lwork = D*D;
    double *work = new double[lwork];
    int info;
    
    dsptrf_(&uplo, &D, a, ipiv, &info);
    if(info != 0)
        cout<<"LU factorization: "<<info<<endl;
    dsptri_(&uplo, &D, a, ipiv, work, &info);
    if(info != 0)
        cout<<"Inversion: "<<info<<endl;

    //Check the inverse
    double **ti;
    ti=new double*[D];
    ct=0;
    for(int i=0; i<D; i++){
        ti[i]=new double[D];
    }
    for(int i=0; i<D; i++){
        for(int j=0; j<=i; j++){
            ti[i][j]=a[ct];
            ti[j][i]=a[ct];
            ct++;
        }
    }
    
    double tp;
    for(int i=0; i<D; i++){
        //cout<<i<<endl;
        for(int j=0; j<D; j++){
            tp=0;
            for(int t=0; t<D; t++){
                tp=tp+correl[i][t]*ti[t][j];
            }
            if(i==j){
                if(tp>1.0001 || tp < 0.9999){
                }
            }
            else{
                if(tp>0.0001 || tp < -0.0001){
                }
            }
        }
    }
    
    inv =new double*[L*Q];
    for(int i=0; i<L*Q; i++){
        inv[i]=new double[L*Q];
        for(int j=0; j<L*Q; j++){
            inv[i][j]=0;
        }
    }
    
    //Take the exponent of the matrix entries and include the lines corresponding to gaps
    ct=0;
    for(int l1=0; l1<L; l1++){
        for(int l2=0; l2<L; l2++){  //No need to compute the diagonal elements since they are not used later on
            if(l1 != l2){
                for(int q1=0; q1<tQ; q1++){
                    for(int q2=0; q2<tQ; q2++){
                        inv[l1*Q+q1][l2*Q+q2]=exp(-ti[l1*tQ+q1][l2*tQ+q2]);
                    }
                }
            }
            // Fill with 1 the values corresponding to the gaps
            for(int q=0; q<Q; q++){
                inv[l1*Q+q][l2*Q+Q-1]=1;
                inv[l1*Q+Q-1][l2*Q+q]=1;
            }
        }
    }
    delete[] ipiv;
    for(int i=0; i<D; i++){
        delete[]  ti[i];
    }
    delete[] ti;
    delete[] a; 
    delete[] work;  
}

//Compute the frequencies and the paired frequencies. Include a random count
void compute_freq(){
    
    RC=0;
    for(int m=0; m<M; m++){
        RC=RC+w[m];
    }
    cout<<"Effective number of sequences: "<<RC<<endl;
    FILE *F;
    char file [4096];
    
    //L=LProteinDomain+Lpeptide;
    cout<<"Compute the aa frequencies at each positions"<<endl;
    fr=new double*[L];
    for(int l=0; l<L; l++){
        fr[l]=new double[Q];
        for(int q=0; q<Q; q++){
            fr[l][q]=0;
        }
    }
    
   // cout<<"Here?"<<endl;
    for(int m=0; m<M; m++){
        //cout<<"m="<<m<<endl; 
        for(int l=0; l<LProteinDomain; l++){
        // cout<<"w[m]="<<w[m]<<endl;  
        // cout<<"l= "<<l<<"; fr[l][seq[m][l]]="<<fr[l][seq[m][l]]<<endl;  
            fr[l][seq[m][l]]=fr[l][seq[m][l]]+w[m];
        }
        
        //cout<<"L="<<L<<endl;
        //cout<<"m="<<m<<endl;
        for(int l=LProteinDomain; l<L; l++){
            for(int n=0; n<nl[m]; n++){
               // cout<<"n="<<n<<endl;
               // cout<<"l="<<l<<endl;
               // cout<<"w[m]="<<w[m]<<endl;
               //cout<<"Number ligand per sequence(m)= nl[m]="<<nl[m]<<endl;
               //cout<<"l-LProteinDomain="<<l-LProteinDomain<<endl;
               //cout<<"lseq[m][n][l-LProteinDomain]="<<lseq[m][n][l-LProteinDomain]<<endl;
               if(nl[m] ==0){ cout<<"Sequence has no ligand, m="<<m<<endl;}
               if(w[m] ==0){ cout<<"Weight is equal to 0 for m="<<m<<endl;}
               fr[l][lseq[m][n][l-LProteinDomain]]=fr[l][lseq[m][n][l-LProteinDomain]]+w[m]*1.0/nl[m];
               
            }
        }
    }
    
    //cout<<"Add pseudo count"<<endl;
    for(int l=0; l<L; l++){
        for (int q=0; q<Q; q++){
            fr[l][q]=fr[l][q]+RC/(1.0*Q);
        }
    }
   // cout<<"Normalize the frequencies"<<endl;
   
    double tot=0;
    for(int l=0; l<L; l++){
        tot=0;
        for (int q=0; q<Q; q++){
            tot=tot+fr[l][q];
        }
        for (int q=0; q<Q; q++){
            fr[l][q]=fr[l][q]/tot;
        }
    }

    //cout<<"Compute the pairwise frequencies"<<endl;
    //******************************//
    //Compute the pairwise frequencies
    //******************************//
    
    double **f12;
    f12=new double*[Q];
    for(int q=0; q<Q; q++){
        f12[q]=new double[Q];
    }
    
    double **fMI;
    fMI=new double*[L];
    for(int l=0; l<L; l++){
        fMI[l]=new double[L];
    }

    correl=new double*[L*tQ];
    for(int l=0; l<L; l++){
        for(int q=0; q<tQ; q++){
            correl[l*tQ+q]=new double[L*tQ];
        }
    }
    
    int mkd;
    sprintf(file, "%s/EC_matrices", family_dir);
    mkd = mkdir(file, S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
    sprintf(file, "%s/MI_matrices", family_dir);
    mkd = mkdir(file, S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
    cout<<"Number of MHC sequence: M= "<<M<<endl;

    for(int l1=0; l1<L; l1++){
        for(int l2=l1; l2<L; l2++){
          
            for(int q1=0; q1<Q; q1++){
                for(int q2=0; q2<Q; q2++){
                    f12[q1][q2]=0;
                }
            }
            if(l1<LProteinDomain && l2<LProteinDomain){
                for(int m=0; m<M; m++){
                    f12[seq[m][l1]][seq[m][l2]]=f12[seq[m][l1]][seq[m][l2]]+w[m];
                }
            }
            if(l1<LProteinDomain && l2>=LProteinDomain){
                for(int m=0; m<M; m++){
                    for(int i=0; i<nl[m]; i++){
                         f12[seq[m][l1]][lseq[m][i][l2-LProteinDomain]]=f12[seq[m][l1]][lseq[m][i][l2-LProteinDomain]]+w[m]*1.0/nl[m];
                      }
                }
            }
            if(l1>=LProteinDomain && l2>=LProteinDomain){
                for(int m=0; m<M; m++){
                    for(int i=0; i<nl[m]; i++){
                        f12[lseq[m][i][l1-LProteinDomain]][lseq[m][i][l2-LProteinDomain]]=f12[lseq[m][i][l1-LProteinDomain]][lseq[m][i][l2-LProteinDomain]]+w[m]*1.0/nl[m];
                    }  
                }
            }
            //Add the random count
            if(l1 != l2){
                for(int q1=0; q1<Q; q1++){
                    for(int q2=0; q2<Q; q2++){
                        f12[q1][q2]=f12[q1][q2]+RC/(1.0*Q*Q);
                    }
                }
            }
            if(l1==l2){
                //Only add the random count on q1=q2, since f12[q1][q2]==0 for q1 != q2
                for(int q1=0; q1<Q; q1++){
                    for(int q2=0; q2<Q; q2++){
                        if(q1 != q2 && f12[q1][q2] > 0){
                            cout<<"Problem\n";
                        }
                    }
                    f12[q1][q1]=f12[q1][q1]+RC/(1.0*Q);
                }
            }
            //  cout<<"Normalize the pairwise frequencies"<<endl;
            //Normalize the pairwise frequencies
            tot=0;
            for(int q1=0; q1<Q; q1++){
                for(int q2=0; q2<Q; q2++){
                    tot=tot+f12[q1][q2];
                }
            }
            
            for(int q1=0; q1<Q; q1++){
                for(int q2=0; q2<Q; q2++){
                    f12[q1][q2]=f12[q1][q2]/tot;
                } 
            }  
            
            if(l1 != l2){
                fMI[l1][l2]=0;
                for(int q1=0; q1<Q; q1++){
                    for(int q2=0; q2<Q; q2++){
                        if(f12[q1][q2]>0){
                            fMI[l1][l2]=fMI[l1][l2]+f12[q1][q2]*log(f12[q1][q2]/(fr[l1][q1]*fr[l2][q2]));
                        }
                    }
                }
            }
            
            if(l1 < LProteinDomain && l2 >= LProteinDomain){
                sprintf(file, "%s/MI_matrices/%d_%d.txt", family_dir, seq_pos[str_pos[l1]], seq_pos[str_pos[l2]]);
                F=fopen(file, "w");
                fprintf(F, "%.5f\t", fMI[l1][l2]);
                for(int q1=0; q1<Q; q1++){
                    fprintf(F, "%c\t", letter[q1]);
                }
                fprintf(F, "\n");
                for(int q1=0; q1<Q; q1++){
                    fprintf(F, "%c\t", letter[q1]);
                    for(int q2=0; q2<Q; q2++){
                        if(f12[q1][q2]>0){
                            fprintf(F, "%.5f\t", f12[q1][q2]*log(f12[q1][q2]/(fr[l1][q1]*fr[l2][q2])));
                        } else{
                            fprintf(F, "0");
                        }
                    }
                    fprintf(F, "\n");
                }
                fclose(F);
            }
            for(int q1=0; q1<tQ; q1++){
                for(int q2=0; q2<tQ; q2++){
                    correl[l1*tQ+q1][l2*tQ+q2]=f12[q1][q2]-fr[l1][q1]*fr[l2][q2];
                    correl[l2*tQ+q2][l1*tQ+q1]=correl[l1*tQ+q1][l2*tQ+q2];
                }
            } 
        }
    }
    
    int pos1, pos2;
    char let1, let2;
    char chain1, chain2;
    
    sprintf(file, "%s/MI_list.txt", family_dir);
    F=fopen(file, "w");
    
    for(int l1=0; l1<L; l1++){
        for(int l2=l1+1; l2<L; l2++){
            cout <<"l1="<<l1<<" l2="<<l2<<endl;
            if(l1<LProteinDomain){
                let1=letter[seq[pref][l1]];
                pos1=seq_pos[str_pos[l1]];
                chain1='A';
            }
            if(l1>=LProteinDomain){
                let1=letter[lseq[pref][lref][l1-LProteinDomain]];
                pos1=seq_pos[str_pos[l1]];
                chain1='B';
            }
            if(l2<LProteinDomain){
                let2=letter[seq[pref][l2]];
                pos2=seq_pos[str_pos[l2]];
                chain2='A';
            }
            if(l2>=LProteinDomain){
                let2=letter[lseq[pref][lref][l2-LProteinDomain]];
                pos2=seq_pos[str_pos[l2]];
                chain2='B';
            }  
            fprintf(F, "%c\t%c\t%c%d\t%c%d\t%.5f\n",  let1, let2, chain1, pos1, chain2, pos2, fMI[l1][l2]);
        }	
    }
    fclose(F);
    for(unsigned i=0; i<Q; i++){
        delete[]  f12[i];
    }
    for(unsigned i=0; i<L; i++){
      delete[] fMI[i];
    }
    delete[] fMI;
    delete[] f12;
}

void import_domain_seq(){
    
    char cline[4096];
    char *str, *str1;
    int mkd;
    str=new char[4096];
    str1=new char[4096];
    mkd = mkdir(family_dir, S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
    cout<<"Store data in implemented data structures for protein tM="<<tM<<endl;
    tseq=new int*[tM];
    tname=new char*[tM];
    
    for(int i=0; i<tM; i++){
        tseq[i]=new int[tLProteinDomain];
        tname[i]=new char[4096];
    }
    for(int m=0; m<tM; m++){
    for(int l=0; l<4096; l++){
        tname[m][l]='\0';
        }
    }
    int tref;
    int tm=0;
    struct sequence *seqReader;
    seqReader = &(seqList);
    while (seqReader !=NULL)
    {
        if(strcmp(seqReader->label,ref_name)==0)
        {
            tref=tm;
            pref=tref;
            tLProteinDomain=strlen(seqReader->sequence);
        }
        for(int i=0; i<strlen(seqReader->label); i++){
            tname[tm][i]=seqReader->label[i];
        }
        for(int l=0; l<tLProteinDomain; l++){
            tseq[tm][l]=position(seqReader->sequence[l]);
       }
        tm++;
        seqReader = seqReader->next;
    }
    /*for(int i=0; i<tM; i++){
        cout<<"DEBUG:i :"<<i<<"="<<tname[i]<<endl;
    }*/
    letter=new char[Q+1];
    if(DNA==1){
        strcpy(letter, "ACGT-");
    }
    if(DNA==0){
        strcpy(letter, "ACDEFGHIKLMNPQRSTVWY-");
    }
    for(int i=0; i<tLProteinDomain; i++){
        cout<<letter[tseq[tref][i]];
    }
    cout<<endl;
    str_pos=new int[100000];
    for(int l=0; l<tLProteinDomain; l++){
        str_pos[l]=l;
    }
    remove_columns(tseq,tLProteinDomain);
    for(int m=0; m<M; m++){
        for(int l=0; l<LProteinDomain; l++){
            if(seq[m][l] < 0 ){
                cout<<"Problem "<<l<<" "<<name[m]<<" "<<seq[m][l]<<"\n";
            }
        }
    }
    cout<<"Number of final positions in protein: "<<LProteinDomain<<endl;
    cout<<"Number of final protein sequences: "<<M<<endl;
    delete[]str;
    delete[]str1;
}

void remove_columns(int **tseq, int tLProteinDomain){
    seq=new int*[tM];
    tseq2=new int*[tM];
    name=new char*[tM];
    int stop=0;
    for(int m=0; m<tM; m++){
        tseq2[m]=new int[tLProteinDomain];
    }
    double gap_fr=0.7; //maximal frequency of gaps
    double tot;
    int ct=0;
    int elim=0;
    M=0;
    for(int m=0; m<tM; m++){
        stop=0;
        for(int l=0; l<tLProteinDomain; l++){
            if(tseq[m][l] == -1){
                stop=1;
            }
        }
        if(stop==0 ){
            name[M]=new char[4096];
            for(int l=0; l<4096; l++){
            name[M][l]=tname[m][l];
            }
            ct=0;
    
            for(int l=0; l<tLProteinDomain; l++){
                tseq2[M][ct]=tseq[m][l];
                ct++;
            }
            M++;
        }
    }
    
    int cl=0;
    for(int l=0; l<tLProteinDomain; l++){
        //Compute the gap frequency
        tot=0;
        for(int m=0; m<M; m++){
            if(tseq2[m][l]==tQ){
                tot++;
            }
        }
        if(tot/(1.0*M) <= gap_fr){
            cl++;
        }else{
        }
    }
    
    LProteinDomain=cl;
    cout<<"DEBUG: remove_columns-Part2 where number of positions is LProteinDomain="<<LProteinDomain<<" out of tLProteinDomain="<<tLProteinDomain<<endl;
    for(int m=0; m<M; m++){
        seq[m]=new int[LProteinDomain];
    }
    seq_pos=new int[100000];
    cout <<"remove_columns: number of protein sequences tM= "<<tM<<endl;
    cout  << "tLProteinDomain (?)= "<<tLProteinDomain<<endl;
    cout  << "tM (Number of sequences before filtering)= "<<tM<<endl;
    cout  << "M (Number of sequences) ="<<M<<endl;

    //RMK: tL; //Number of positions before removing some columns in the alignment
    //RMK: tM; //Number of sequences before filtering
    //RMK: lengthM; //Number of lines in alignment file
    
    ct=0;
    for(int l=0; l<tLProteinDomain; l++){
        tot=0;
        elim=0;
        for(int m=0; m<M; m++){
            //if(m==tref){pref=m;}
            if(tseq2[m][l]==tQ){
                tot++;
            }
        }
        if(tot/(1.0*M) <= gap_fr){
            for(int m=0; m<M; m++){
                seq[m][ct]=tseq2[m][l];
                seq_pos[ct]=l+1;
            }
            ct++;
        }
    }
     
    cout<<"DEBUG:LProteinDomain="<<LProteinDomain<<endl;
    for(int l1=0; l1<LProteinDomain; l1++){
      cout<<"seq_pos["<<l1<<"]="<< seq_pos[l1]<<endl;
    }
    
     cout<<"DEBUG:LProteinDomain="<<LProteinDomain<<endl;
    for(int l1=0; l1<LProteinDomain; l1++){
      cout<<"str_pos["<<l1<<"]="<< str_pos[l1]<<endl;
    }
    
    
}

void import_ligands(int pref, int lref, int LProteinDomain){
    int *dom;
    char *cs;
    string line;
    char *file;
    string nada="";
    int *tk;
    
    nl    =new int[M];
    lseq=new int**[M];  //ligand sequence
    tk=new int[M];      //tokenizer
    file=new char[4096]; 
    for(int m=0; m<M; m++){
        nl[m]=0;
    }
    dom=new  int[Mligand];
    cs =new  char[Mligand];
    for(int n=0; n<Mligand; n++){
        dom[n]=0;
    }    
    struct sequence *seqReaderLigand;
    seqReaderLigand = &(seqListLigand);
    Mligand=0;
    while (seqReaderLigand !=NULL)
    {
        int tst=0;
        for(int m=0; m<M ; m++){
            if(strcmp(name[m], seqReaderLigand->label) == 0 ){
                
                dom[Mligand]=m;
                nl[m]= nl[m]+1;
                tst=1;  
              }
        }
        if(tst==0){ 
        cout <<"Error: not added! "<<seqReaderLigand->label<<" "<<seqReaderLigand->sequence <<endl;
        }
        else{
                Mligand++;
        }
        seqReaderLigand = seqReaderLigand->next;
        
    }
    for(int m=0; m<M; m++){
        lseq[m]=new int*[nl[m]];
	tk[m]=0;
    }
    
    seqReaderLigand = &(seqListLigand);
    cout<<"Mligand:"<<Mligand<< endl;
    for(int i=0; i<Mligand; i++){
        int HLASeqID=dom[i];
       
        int temp=tk[HLASeqID];
        lseq[HLASeqID][temp]=new int[tLpeptide];
        for(int l=0; l<strlen(seqReaderLigand->sequence); l++){
            lseq[HLASeqID][temp][l]=position(seqReaderLigand->sequence[l]);
        }
        seqReaderLigand = seqReaderLigand->next;
        tk[HLASeqID]=temp+1;
    }
  
    Lpeptide=tLpeptide;
    for(int l=0; l<Lpeptide; l++){
    	seq_pos[LProteinDomain+l]=l+1;
    }
    for(int l=0; l<Lpeptide; l++){
        str_pos[l+LProteinDomain]=LProteinDomain+l;

    }
    L=LProteinDomain+Lpeptide;
    cout<<"DEBUG:L="<<L<<endl;
    for(int l1=0; l1<L; l1++){
      cout<<"seq_pos["<<l1<<"]="<< seq_pos[l1]<<endl;
    }
    
     cout<<"DEBUG:L="<<L<<endl;
    for(int l1=0; l1<L; l1++){
      cout<<"str_pos["<<l1<<"]="<< str_pos[l1]<<endl;
    }
   /* 
   cout<<"DEBUG..."<<endl;
    for(int i=0; i<M; i++){
        for(int iii=0; iii<9; iii++){
               cout<<"lseq["<<i<<"][0]["<<iii<<"]="<<lseq[i][0][iii]<<endl;
             
        }
    }*/
    
    cout<<"Number of final positions in peptide: "<<Lpeptide<<endl;
    cout<<"Number of final peptide sequences: "<<Mligand<<endl;
    delete[]dom;
    delete[]cs;
    delete[]file;
  
}

void weights(){
    w=new double[M];
    double sim_thresh=0.7;
    double tot;
    double sim;
    int comp_weight=1;
    if(comp_weight==1){
        for(int m=0; m<M; m++){
            w[m]=1;
        }
        for(int m=0; m<M; m++){
            for(int m1=m+1; m1<M; m1++){
                tot=0;
                sim=0;
                for(int l=0; l<LProteinDomain; l++){
                    tot++;
                    if(seq[m][l] == seq[m1][l] ){
                        sim++;
                    }
                }
                if(1.0*sim/tot >= sim_thresh){
                    w[m]++;
                    w[m1]++;
                }
            }
        }
        for(int m=0; m<M; m++){
            w[m]=1.0/w[m];
        }
        FILE *F;
        char file [4096];
        sprintf(file, "%s/weight.txt", family_dir);
        F=fopen(file, "w");
        for(int m=0; m<M; m++){
            fprintf(F, "%.5f\n", w[m]);
        }
        fclose(F); 
    }
    else{
        fstream afile;
        char file [4096];
        sprintf(file, "%s/weight.txt", family_dir);
        afile.open(file, ios::in);
        for(int m=0; m<M; m++){
            afile>>w[m];
        }
        afile.close();
    }
}

int position(char s){
    int pp=0;
    //We are working with peptides
    if(DNA==0){
        if(s == 'A' || s == 'a')
            pp=0;
        if(s == 'C' || s == 'c')
            pp=1;
        if(s == 'D' || s == 'd')
            pp=2;
        if(s == 'E' || s == 'e')
            pp=3;
        if(s == 'F' || s == 'f')
            pp=4;
        if(s == 'G' || s == 'g')
            pp=5;
        if(s == 'H' || s == 'h')
            pp=6;
        if(s == 'I' || s == 'i')
            pp=7;
        if(s == 'K' || s == 'k')
            pp=8;
        if(s == 'L' || s == 'l')
            pp=9;
        if(s == 'M' || s == 'm')
            pp=10;
        if(s == 'N' || s == 'n')
            pp=11;
        if(s == 'P' || s == 'p')
            pp=12;
        if(s == 'Q' || s == 'q')
            pp=13;
        if(s == 'R' || s == 'r')
            pp=14;
        if(s == 'S' || s == 's')
            pp=15;
        if(s == 'T' || s == 't')
            pp=16;
        if(s == 'V' || s == 'v')
            pp=17;
        if(s == 'W' || s == 'w')
            pp=18;
        if(s == 'Y' || s == 'y')
            pp=19;
        if(s == 'X' || s == '-' || s == 'x' || s == '*' || s == '.')
            pp=tQ;

    }
    if(DNA==1){
        
        if(s == 'A' || s == 'a')
            pp=0;
        if(s == 'C' || s == 'c')
            pp=1;
        if(s == 'G' || s == 'g')
            pp=2;
        if(s == 'T' || s == 't')
            pp=3;
        if(s == 'X' || s == '-' || s == 'x' || s == '*' || s == '.')
            pp=tQ;
    }
    
    if(DNA==2){
        
        if(s == 'A' || s == 'a')
            pp=0;
        if(s == 'C' || s == 'c')
            pp=1;
        if(s == 'G' || s == 'g')
            pp=2;
        if(s == 'U' || s == 'u')
            pp=3;
        if(s == 'X' || s == '-' || s == 'x' || s == '*' || s == '.')
            pp=tQ;
    }
    return(pp);
    
}