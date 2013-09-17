#include <cstdio>
#include <cstdlib>

char FOLDER[50];
char JOB[50];
char RFFN[100];

void CreateJobFile(double U, int StartFromSC)
{
  char FN[50];
  sprintf(FN, "%s/%s",FOLDER,JOB);
  FILE* ParamsFile = fopen(FN,"w");
   
  fprintf(ParamsFile, "#!/bin/bash\n");
  fprintf(ParamsFile, "#PBS -q hpsee\n");
  fprintf(ParamsFile, "#PBS -l nodes=1:ppn=1\n");
  fprintf(ParamsFile, "#PBS -e U%.3f.err\n",U);
  fprintf(ParamsFile, "#PBS -o U%.3f.out\n",U);
  fprintf(ParamsFile, "\n");
  fprintf(ParamsFile, "cd %s/run\n",FOLDER);
  fprintf(ParamsFile, "chmod u+x main\n");
  fprintf(ParamsFile, "./main %.3f %d\n",U,StartFromSC);

  fclose(ParamsFile);
}

int main()
{
  sprintf(FOLDER,"/nfs/jaksa/WL_IPT");
  
  for (double U=0.05; U<0.5; U+=0.05)  
  {
    sprintf(JOB,"ipt.U%.3f.pbs",U);
    char command[50];
    CreateJobFile(U,0); 
    sprintf(command,"qsub %s/%s",FOLDER, JOB);
    system(command);
  }
/*
  for (double U=2.5; U<3.5; U+=0.05)  
  {
    sprintf(JOB,"ipt_FSC.U%.3f.pbs",U);
    char command[50];
    CreateJobFile(U,1); 
    sprintf(command,"qsub %s/%s",FOLDER, JOB);
    system(command);
  }
*/
  return 0;

}
