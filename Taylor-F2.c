/*
Taylor-F2 Waveform Generator
Furqan Dar
20th March, 2016
*/




#include "Taylor-F2-Functions.h"






int main()
{
    printf("\n\n************************************************\n\nStarting Program.\n");


/*
Dynamically allocate memory to each of the four global arrays using malloc. 
Using sizeof(*Array) makes it so that malloc will automatically adjust for the type of array initialized with the pointer in Taylor-F2-Functions.h
Each array has it's own error check to hard confirm if malloc was successful or not.
*/

PSDArr = malloc(sizeof(*PSDArr)*DAT_ARR);
if (PSDArr != NULL)
{
    printf("Memory allocated for the Power Spectral Density (PSD) array.\n");
}
else
{
    printf("Memory allocation failed for PSD Array! Exiting program.\n");
    exit(1);
}


DatArr = malloc(sizeof(fComp)*DAT_ARR);
if (DatArr != NULL)
{
printf("Memory allocated for the imported data array.\n");
}
else
{
    printf("Memory allocation failed for imported data Array! Exiting program.\n");
    exit(1);
}

MCMCArr1 = malloc(sizeof(fComp)*DAT_ARR);
if (MCMCArr1 != NULL)
{
printf("Memory allocated for the first MCMC array.\n");
}
else
{
    printf("Memory allocation failed for first MCMC Array! Exiting program.\n");
    exit(1);
}

MCMCArr2 = malloc(sizeof(fComp)*DAT_ARR);
if (MCMCArr2 != NULL)
{
printf("Memory allocated for the second MCMC array.\n\n\n");
}
else
{
    printf("Memory allocation failed for second MCMC Array! Exiting program.\n");
    exit(1);
}



    fNum FarR = 350.;


    //Generate the fake data with desired parameters and export to a file: FakeData.txt
    TF2_Printer(1.4, 1.4, 40., FarR, 0., 0., 0., 7);
    
    //Import the generated fake data from FakeData.txt
    DataImporter();

    PSDImporter();

    NoiseImporter();

    NoiseMixer();
    

    MassMCMC(50000, 1.4, 1.4, FarR, 7);
    
    free(PSDArr);
    free(DatArr);
    free(MCMCArr1);
    free(MCMCArr2);
    printf("\n\n************************************************\n\nProgram ran successfully.\n");

    return 0;
}




/*

    fNum DatArr[DAT_ARR] = {0.};
    FILE *fPoint;
    fPoint = fopen("InitDat.txt","r");
    fNum amp = 0.;
    fNum phase = 0.;
    for (int i = 0; i < DAT_ARR; i++)
    {
        fscanf(fPoint, "%le\t%le", &amp, &phase);
        DatArr[i] = amp*cexp(-I*phase);
        //printf("%le %le %le %le\n", carg(DatArr[i]), phase,2.*PI, cabs(DatArr[i]), amp );
    }
    
    fclose(fPoint);

    printf("Imported 'Fake' Data. Initialized and filled'fake' data array.\n");

    fNum PSDINComp[8192] = {0.};
    fNum PSDINRed[DAT_ARR] = {0.};
    FILE *psdPoint;
    psdPoint = fopen("AvePSD.txt","r");

    for (int i = 0; i < 8192; i++)
    {
        fscanf(psdPoint, "%le", &PSDINComp[i]);
        //printf("%le\n", PSDINComp[i]);
    }

    printf("Imported full ASD.\n");

    for (int i = 0; i < DAT_ARR; i++)
    {
        PSDINRed[i] = PSDINComp[319+i];
    }

printf("Initialized and filled reduced ASD array.\n");


    
    MassMCMC(10, 1.5,1.3, DatArr, PSDINRed, 1);

*/