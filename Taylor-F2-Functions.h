/*
Taylor-F2 Waveform Generator - Functions File
Furqan Dar
1st March, 2016
*/

/*
Import the required libraries. If list gets too big, will simply make personal .h file that has all the libraries.
*/

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <time.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <complex.h>
#include <omp.h>


#define DAT_ARR 7872


typedef complex double fComp; //Defines the global complex float type being used
typedef double fNum;//Defines the global float type being used
typedef int fIdx; //Defines the global int type being used

/*
Define global constants for rescalings. If list gets too big, might make .h files containing constants. (Especially since multiple files might need to share these)
Using doubles and not long doubles because gsl-functions is not fully compatible with long doubles
*/



const fNum PI = 3.14159265359;//Might change to use M_PI
const fNum LIGHT_C = 1.;
const fNum LIGHT_C_SI = 299792458.;
const fNum GRAV_G = 1.;
const fNum GRAV_G_SI = 6.67384e-11;
const fNum MPC = 1.0292712502794183e14;
const fNum MPC_SI = 3.086e+22;
const fNum SolMass_KG = 1.9891e+30;
const fNum SolMass = 4.925492321898863e-6;
const fNum EulerGamma = 0.577215664902;
const fNum sample_rate = 2048.; //Might just use the Euler-Gamma from gsl-functions when that phase begins
const fNum dF = 0.125;//0.031250 is the actual dF in the total data. Undersampled for now.;

/*
Declare four GLOBAL arrays. One for the PSD, one for the Data and two for the generated model during the MCMC routine.
Using GLOBAL arrays so that functions don't need to be redundantly passed these arrays.
Furthermore, changing values of these arrays in specific functions becomes much easier.
Remember that all four arrays are fComp type.
(Apparently, not good computer science practice.)
*/

fNum (*PSDArr);
fComp (*DatArr);
fComp (*MCMCArr1);
fComp (*MCMCArr2);




/*
Function definitions
Based HEAVILY on CBCGen-RK.ipynb
*/

fNum TotMass(fNum mass1, fNum mass2)
{//Returns the value for the total mass. Might be redundant, but makes reading easier.
	return mass1 + mass2;

}

fNum RedMass(fNum mass1, fNum mass2)
{//Returns reduced mass ratio.
	return mass1*mass2/TotMass(mass1,mass2);
}

fNum SymMass(fNum mass1, fNum mass2)
{//Returns the value for the symmetric mass ratio.
	return RedMass(mass1,mass2)/TotMass(mass1,mass2);
}

fNum VelFromFreq( fNum freq, fNum mass1, fNum mass2)
{//Returns v as defined by the frequency
	return pow(TotMass(mass1,mass2)*PI*freq,1./3.);
}

fNum fISCO_Calc(fNum mass1, fNum mass2)
{//Returns f_ISCO. Completely determined by masses.
	return 1./pow(6.,3./2.)/PI/TotMass(mass1,mass2);
}

fNum TF2_Amp(fNum m1, fNum m2, fNum FarR, fNum freq, fNum inc)
{//Evaluates and returns the amplitude of TF2 Stationary Phase waveform in frequency space
	fNum mu = SymMass(m1, m2);
	fNum TM = TotMass(m1, m2);
	fNum nCONST = sqrt(5.*PI/24.);
	fNum x = pow(PI*TM*freq, 2./3.);
	fNum incNorm = -(1. + cos(inc)*cos(inc))/2.;
	fNum TotalConst = TM*TM*sqrt(mu)*nCONST*incNorm/FarR;
	fNum Amp = TotalConst/pow(x, 7./4.);
	return Amp;
}

fNum TF2_PsiF(fNum m1, fNum m2, fNum freq, fNum freq0, fIdx ORDER)
{
	fNum mu = SymMass(m1, m2);
	fNum TM = TotMass(m1, m2);
	fNum nCONST = - PI/4.;
	fNum x = pow(PI*TM*freq, 2./3.);
	fNum x0 = pow(PI*TM*fISCO_Calc(m1,m2), 2./3.);
	fNum VarX = 3./128./mu/pow(x, 2.5);
	switch(ORDER)
	{
		case 1:
			return VarX*1. + nCONST;
		case 2:
			return VarX*(1. + (3715./756. + 55.*mu/9.)*x ) + nCONST;
		case 3:
			return VarX*(1. + (3715./756. + 55.*mu/9.)*x - 16.*PI*pow(x,1.5)) + nCONST;
		case 4:
			return VarX*(1. + (3715./756. + 55.*mu/9.)*x - 16.*PI*pow(x,1.5) + (15293365./508032. + 27145.*mu/504. + 3085.*mu*mu/72.)*x*x ) + nCONST;
		case 5:
			return VarX*(1. + (3715./756. + 55.*mu/9.)*x - 16.*PI*pow(x,1.5) + (15293365./508032. + 27145.*mu/504. + 3085.*mu*mu/72.)*x*x + (38645./756. - 65.*mu/9.)*(1. + 3.*log(x/x0)/2.)*PI*pow(x,2.5)) + nCONST;
		case 6:
			return VarX*(1. + (3715./756. + 55.*mu/9.)*x - 16.*PI*pow(x,1.5) + (15293365./508032. + 27145.*mu/504. + 3085.*mu*mu/72.)*x*x + (38645./756. - 65.*mu/9.)*(1. + 3.*log(x/x0)/2.)*PI*pow(x,2.5) + 
								(11583231236531./4694215680. - 640.*PI*PI/3. - 6848.*EulerGamma/21. - 3424.*log(16.*x)/21. + (-15737765635./3048192. + 2255.*PI*PI/12.)*mu + 76055.*mu*mu/1728. - 127825.*mu*mu*mu/1296.)*x*x*x) + nCONST;
		case 7:
			return VarX*(1. + (3715./756. + 55.*mu/9.)*x - 16.*PI*pow(x,1.5) + (15293365./508032. + 27145.*mu/504. + 3085.*mu*mu/72.)*x*x + (38645./756. - 65.*mu/9.)*(1. + 3.*log(x/x0)/2.)*PI*pow(x,2.5) + 
								(11583231236531./4694215680. - 640.*PI*PI/3. - 6848.*EulerGamma/21. - 3424.*log(16.*x)/21. + (-15737765635./3048192. + 2255.*PI*PI/12.)*mu + 76055.*mu*mu/1728. - 127825.*mu*mu*mu/1296.)*x*x*x
									+ (77096675./254016. + 378515.*mu/1512. - 74045.*mu*mu/756.)*PI*pow(x,3.5) ) + nCONST;
		default:
			return 0.;
	}
}

void TF2_Printer(fNum m1, fNum m2, fNum freq0, fNum FarR, fNum inc, fNum phi_c, fNum t_c, fIdx ORDER)
{
	m1 *= SolMass;
	m2 *= SolMass;
	FarR *= MPC;
	fNum nCONST =  2.*PI*(phi_c - t_c);
	fNum freq = freq0;
	fNum phase;
	fNum amp;
	fNum f_ISCO = fISCO_Calc(m1, m2);
	FILE *fPoint;
	fPoint = fopen("FakeData.txt", "w+");
	fIdx i = 0;
	//printf("ISCO Frequency is %f\n", f_ISCO);
	for (i = 0; i< DAT_ARR; i++)
	{
		fNum j = (fNum) i;
		freq = freq0 + j*dF;
		nCONST =  -2.*PI*(phi_c - freq*t_c);
		phase = TF2_PsiF(m1, m2, freq, freq0, ORDER)+nCONST;
		amp = TF2_Amp(m1, m2, FarR, freq, inc);
		
		fprintf(fPoint, "%le\t%le\n", amp, phase);
		DatArr[i]=amp*cexp(-I*phase);
		//printf("%le %le\n",-amp, phase );
	}
	//printf("%d\n", i);
	fclose(fPoint);
	printf("Generated 'Fake' Data.\n");
}

void DataImporter()
{
/*
Function that imports data to DatArr[]. NOTE THAT FILE NAME SHOULD BE FakeData.txt for now. (Might add user flexibility in the future.)
The format of the file being imported from is expected as such
amp(z) \t arg(z)\n
Thus, each line has the amp and phase corresponding to a frequency, delimited by a tab. See TF2_Printer for reference.
For now, the frequency range expedted is from 40 Hz to 1024 Hz.
*/

//Declaring two fNum types to store amplitude and phase.
	fNum amp;
	fNum phase; 

	FILE *fp;
	fp = fopen("FakeData.txt", "r+");
	fIdx i;
	for (i = 0; i < DAT_ARR; i++)
	{
		fscanf(fp,"%le\t%le\n", &amp, &phase);
		DatArr[i]=amp*cexp(-I*phase);
		//printf("%le %le\n", amp,phase);
	}
	fclose(fp);
	printf("Data successfully imported.\n");

}

void NoiseImporter()
{
/*
Function that imports noise in fourier space to PSDArr. Note that PSDArr is being used as a dummy array here and will actually be used to store the PSD. 
NOTE THAT FILE NAMSE SHOULD BE InjNoise.txt. (Might add user flexibility in the future.)
If this function is not called in the main file, the imported data will be 'clean'.
The format of the file being imported from is expected as such
amp(z) \t arg(z)\n
Thus, each line has the amp and phase corresponding to a frequency, delimited by a tab. See TF2_Printer for reference.
For now, the frequency range expedted is from 40 Hz to 1024 Hz.
NOTE THIS NOISE IS TAKEN FROM PSDGen-RandGen juPyter notebook.
*/

//Declaring two fNum types to store amplitude and phase.
fNum amp;
fNum phase; 

FILE *fp;
fp = fopen("InjNoise.txt", "r");
fIdx i;
for (i = 0; i < DAT_ARR; i++)
{
	fscanf(fp,"%le %le\n", &amp, &phase);
	MCMCArr1[i]=amp*cexp(I*phase);
	//printf("%le %le - %le %le\n", amp,cabs(MCMCArr1[i]), phase,carg(MCMCArr1[i]));
}
fclose(fp);
printf("Noise successfully imported.\n");

}

void NoiseMixer()
{
/*
Function that injects the data into the already imported noise.
Order of functions matters since the arrays are global.
PSDImporter should be AFTER this function.
*/
fIdx i;
for (i = 0; i < DAT_ARR; i++)
{	printf("%le %le --->>>", cabs(DatArr[i]), carg(DatArr[i]) );
	DatArr[i]=DatArr[i]+PSDArr[i];
	printf("%le %le\n", cabs(DatArr[i]), carg(DatArr[i]) );
}
printf("Noise successfully mixed with the imported data.\n");

}

void PSDImporter()
{

/*
Function that imports PSD to PSDArr[]. NOTE THAT FILE NAME SHOULD BE PSDRed.txt for now. (Might add user flexibility in the future.)
The format of the file being imported from is expected as such
PSD
Thus, each line has the PSD.
For now, the frequency range expedted is from 40 Hz to 1024 Hz.
*/


FILE *fp;
fp = fopen("PSDRed.txt", "r+");
fIdx i; 
for (i = 0; i < DAT_ARR; i++)
{
	fscanf(fp,"%le\n", &PSDArr[i]);
}
fclose(fp);
printf("PSD successfully imported.\n");

}

fNum SnR_Model_Injected()
{
	fNum L = 0.;
	int i ;
	#pragma omp parallel for private(i) reduction  (+:L)
	for (i = 0; i < DAT_ARR; i++)
	{
		L += cabs((DatArr[i])*conj(MCMCArr1[i])/(PSDArr[i]*PSDArr[i]));
	}
	return L;
}

fNum FreqSnR (fComp xArr[])
{
	fNum L = 0.;
	int i ;
	#pragma omp parallel for private(i) reduction  (+:L)
	for (i = 0; i < DAT_ARR; i++)
	{
		L += cabs((xArr[i])*conj(xArr[i])/(PSDArr[i]*PSDArr[i]));
	}
	return L;
}






fNum WFLike(fComp ModArr[])
{
	fNum L = 0.;
	int i ;
	#pragma omp parallel for private(i) reduction  (+:L)
	for (i = 0; i < DAT_ARR; i++)
	{
		//L += ((-ModArr[i]+DatArr[i])*conj((-ModArr[i]+DatArr[i])))/(PSDArr[i]*PSDArr[i]);
		L += cabs((DatArr[i]-ModArr[i])*conj(DatArr[i]-ModArr[i])/(PSDArr[i]*PSDArr[i]));
		//printf("%le\n", L);
	}
	//printf("%le\n", exp(-2.*L));
	return exp(-2.*L);
}

fNum MassPrior(fNum mass)
{
	if((mass >= 1.3) && (mass <= 1.5))
		{
			return 1.;
		}
	else
	{
		return 0.;
	}
}

void MassMCMC(fIdx N, fNum m1, fNum m2, fNum FarR, fIdx ORDER)
{

	printf("Starting MCMC Routine.\n\n\n");
	//Making two arrays that will hold the unscaled masses. The first component is the initial mass, and the second component is the proposed mass.
	fNum mass1[2] = {m1, 0.};
	fNum mass2[2] = {m2, 0.};
	
	//Making two arrays that contain scaled masses. These are the masses that are input into the WF Generators.
	fNum m1_scaled[2] = {m1*SolMass, 0.};
	fNum m2_scaled[2] = {m2*SolMass, 0.};
	
	//Defining constants that are put into the WF Generators.
	fNum freq = 0.;
	fNum freq0 = 40.;
	FarR *= MPC;
	fNum inc = 0.0;
	fNum phi_c = 0.; 
	fNum t_c = 0.;
	fNum nCONST = 0.;

	//Defining two numbers to store the Posterior. LikeNew is the posterior of the proposed parameters.
	//Also definning the Metropolis ratio number and initializing as 0
	fNum LikeOld = 0.;
	fNum LikeNew = 0.;
	fNum MetR = 0.;

	//Defining and setting up the Random number generating routine.
	//Sigma is the std. dev of the gaussian used to sample random paramaters around the old parameter.
	const gsl_rng_type * T;
	gsl_rng * r;
	r = gsl_rng_alloc (gsl_rng_mt19937);
	fNum RandNum = 0.;
	fNum sigma = 1.0/1.0;

	
	//Initializing iterator (perhaps necessary for OpenMp)
	fIdx k ;
	fNum j ;

	//Opening up a file to which write selected masses to.
	FILE *fPoint;
	fPoint = fopen("MassDist.txt" ,"w");

	//Generating a WF using provided initial conditions.
	#pragma omp parallel for private(k, j, freq, nCONST)
		for (k = 0; k < DAT_ARR; k++)
		{
			j = (fNum) k;
			//printf("%lf %lf\n", j*dF, k*dF );
			freq = freq0 + j*dF;
			nCONST = -2.*PI*(phi_c - freq*t_c);
			MCMCArr1[k] = TF2_Amp(m1_scaled[0], m2_scaled[0], FarR, freq, inc)*cexp(-I*(TF2_PsiF(m1_scaled[0], m2_scaled[0], freq, freq0, ORDER)+nCONST));
			//printf("%lf %lf\n", creal(DatArr[k]-ModArr1[k]) ,cimag(DatArr[k]-ModArr1[k]) );
		}


	//Performing SnR Calculations to see if the proposed initial conditions and imported data are 'good'
	printf("Injected Data with Model Data SnR is %lf\n", SnR_Model_Injected());
	printf("Imported Noisy data SnR is %lf\n", FreqSnR(DatArr));
	printf("Initial conditions Model SnR with itself is %lf\n", FreqSnR(MCMCArr1));
	
	//Calculating the posterior before beginning the MCMC loop to reduce redundant calculations.
	LikeOld = WFLike(MCMCArr1)*MassPrior(mass1[0])*MassPrior(mass2[0]);
	//printf("%lf %lf %lf\n", LikeOld, mass1[0],mass2[0] );

	//Put in flag for Initial Values
	fIdx i;
	for (i = 0; i < N; i++)
	{
		//Generating a normally distrubed random number.
		//Using that random number to propose a new mass, and then scaling it to be used in WF Generation.
		RandNum = gsl_rng_uniform(r);
		if (RandNum >= 0.5)
		{
		RandNum = gsl_ran_gaussian(r, sigma);
		mass1[1] = mass1[0] + RandNum;
		m1_scaled[1] = mass1[1]*SolMass;
		}

		//Read above comment. Doing the same for the second mass.
		else
		{
		RandNum = gsl_ran_gaussian(r, sigma);
		mass2[1] = mass2[0] + RandNum;
		m2_scaled[1] = mass2[1]*SolMass;
		}

	//Generating WF with proposed masses.
	#pragma omp parallel for private(k, j, freq, nCONST)
		for (k = 0; k < DAT_ARR; k++)
		{
			j = (fNum) k;
			freq = freq0 + j*dF;
			nCONST = -2.*PI*(phi_c - freq*t_c);
			MCMCArr2[k] = TF2_Amp(m1_scaled[1], m2_scaled[1], FarR, freq, inc)*cexp(-I*(TF2_PsiF(m1_scaled[1], m2_scaled[1], freq, freq0, ORDER)+nCONST));
			//printf("%lf %lf\n", creal(ModArr1[k]-ModArr2[k]) ,cimag(ModArr1[k]-ModArr2[k]) );
		}


	//Calculating the posterior of the proposed masses.
	LikeNew = WFLike(MCMCArr2)*MassPrior(mass1[1])*MassPrior(mass2[1]);
	//printf("%le %le\n", LikeOld, LikeNew );

	//Calculating the Metropolis ratio.
	MetR = LikeNew/LikeOld;
	//printf("%lf\n",MetR);
	
	//Generating a random number from U(0,1)
	RandNum = gsl_rng_uniform(r);

	//Performing the check.
	if (MetR > RandNum)
		{
			//printf("%lf %lf\n", mass1[0], mass2[0]);
			mass1[0] = mass1[1];
			mass2[0] = mass2[1];
			LikeOld = LikeNew;
			fprintf(fPoint, "%lf %lf\n", mass1[0],mass2[0]);
			//printf("%lf %lf\n", mass1[0],mass2[0]);
			sigma *= 8.;
		}
	else
		{
			sigma /= 2.;
		}
		
		//printf("%lf %lf %lf %lf %lf\n",mass1[0], mass1[1], mass2[0], mass2[1], MetR );
	}
		fclose(fPoint);
		printf("\n\n\nMCMC Routine Finished.\n");

}


/*


fNum WFLike(fNum ModArr[], fNum DArr[], fNum PSD[])
{
	fNum L = 0.;
	for (int i = 0; i < DAT_ARR; i++)
	{
		L += ((-ModArr[i]+DArr[i])*conj((-ModArr[i]+DArr[i])))/(PSD[i]);
		//printf("%le %le %le %le %le %le\n", creal(ModArr[i]), cimag(ModArr[i]),creal(DArr[i]), cimag(DArr[i]), PSD[i], L );
	}
	printf("%le\n", -log(2.*L) );
	return -(2.*L);
}

fNum MassPrior(fNum mass)
{
	if((mass <= 2.0) && (mass >= 1.0))
		{
			return 1.;
		}
	else
	{
		return 0.;
	}
}



void MassMCMC(fIdx N, fNum m1, fNum m2, fNum DatArr[], fNum PSD[], fIdx ORDER)
{

	printf("Starting MCMC Routine.\n");
	fNum mass1[2] = {m1, 0.};
	fNum mass2[2] = {m2, 0.};
	fNum m1_scaled[2] = {m1*SolMass,0.};
	fNum m2_scaled[2] = {m2*SolMass,0.};
	fNum ModArr1[DAT_ARR] = {0.};
	fNum ModArr2[DAT_ARR] = {0.};
	fNum freq = 0.;
	
	fNum freq0 = 40.;
	fNum FarR = 1.0;
	FarR *= MPC;
	fNum inc = 0.0;
	fNum phi_c = 0.; 
	fNum t_c = 0.;
	fNum nCONST = 0;
	fNum LikeOld = 0.;
	fNum LikeNew = 0.;

	const gsl_rng_type * T;
	gsl_rng * r;
	r = gsl_rng_alloc (gsl_rng_mt19937);
	fNum RandNum = 0.;

	fNum sigma = 1.0;
	fNum MetR = 0.;

	FILE *fPoint;
	fPoint = fopen("MassDist.txt" ,"w");

	for (int k = 0; k < DAT_ARR; k++)
		{
			fNum j = (fNum)k;
			freq = freq0 + j*dF;
			nCONST = -2.*PI*(phi_c - freq*t_c);
			ModArr1[k] = TF2_Amp(m1_scaled[0], m2_scaled[0], FarR, freq, inc)*cexp(-I*(TF2_PsiF(m1_scaled[0], m2_scaled[0], freq, freq0, ORDER)+nCONST));
			//printf("%lf %lf\n", creal(DatArr[k]-ModArr1[k]) ,cimag(DatArr[k]-ModArr1[k]) );
		}
		LikeOld = WFLike(ModArr1, DatArr, PSD)*MassPrior(mass1[0])*MassPrior(mass2[0]);
		//printf("%lf %lf %lf\n", LikeOld, mass1[0],mass2[0] );

	for (int i = 0; i < N; i++)
	{
		RandNum = gsl_ran_gaussian(r, sigma);
		mass1[1] = mass1[0] + RandNum;
		m1_scaled[1] = mass1[1]*SolMass;
		RandNum = gsl_ran_gaussian(r, sigma);
		mass2[1] = mass2[0] + RandNum;
		m2_scaled[1] = mass2[1]*SolMass;

		for (int k = 0; k < DAT_ARR; k++)
		{
			fNum j = (fNum)k;
			freq = freq0 + j*dF;
			nCONST = -2.*PI*(phi_c - freq*t_c);
			ModArr2[k] = TF2_Amp(m1_scaled[1], m2_scaled[1], FarR, freq, inc)*cexp(-I*(TF2_PsiF(m1_scaled[1], m2_scaled[1], freq, freq0, ORDER)+nCONST));
			//printf("%lf %lf\n", creal(ModArr1[k]-ModArr2[k]) ,cimag(ModArr1[k]-ModArr2[k]) );
		}
		RandNum = gsl_rng_uniform(r);
		LikeNew = WFLike(ModArr2, DatArr, PSD)*MassPrior(mass1[1])*MassPrior(mass2[1]);
		//printf("%le %le\n", LikeOld, LikeNew );
		MetR = LikeNew/LikeOld;
		//printf("%lf\n",MetR);

		if (MetR > RandNum)
		{
			//printf("%lf %lf\n", mass1[0], mass2[0]);
			mass1[0] = mass1[1];
			mass2[0] = mass2[1];
			LikeOld = LikeNew;
			fprintf(fPoint, "%lf %lf\n", mass1[0],mass2[0]);
			sigma *= 8.;
		}
		sigma /= 2.;
		//printf("%lf %lf %lf %lf %lf\n",mass1[0], mass1[1], mass2[0], mass2[1], MetR );
	}
		fclose(fPoint);
		printf("MCMC Routine Finished.\n");

}

/*

fNum factorial(int x)
{
	fNum DumV = 1.;
	for (int i = x; i > 0; i--)
	{
		DumV *= i;
	}
	return DumV;

}



fNum PoissonLike(int x, int y, int l)
{	if (y >= 0)
	{
		return pow(l, y-x)*factorial(x)/factorial(y);
	}
	else
	{
		return 0.;
	}
}


void PoissonEx(int N, int lambda, int x0)
{
	const gsl_rng_type * T;
	gsl_rng * r;
	r = gsl_rng_alloc (gsl_rng_mt19937);
	fNum RandNum = 0.;

	fIdx x = x0;
	fIdx y = 0;
	fNum MetR = 0.;
	for (int i = 0; i < N; i++)
	{
		RandNum = gsl_rng_uniform(r);
		if (RandNum > 0.5)
		{
			y = x + 1;
		}
		else 
		{
			y = x - 1;
		}

		MetR = PoissonLike(x , y, lambda);
		RandNum = gsl_rng_uniform(r);

		if ( RandNum <= MetR)
		{
			x = y;
			printf("%d\n", x);
		}
	}

}


fNum BiGaussLike(fNum xnew, fNum xold, fNum ynew, fNum yold)
{
	return exp((-((-4. + xnew)*xnew) + (-4. + xold)*xold - (ynew - yold)*(2. + ynew + yold))/2.);
}

fNum BiGaussLike2(fNum xnew, fNum xold, fNum ynew, fNum yold)
{
	return exp((-xnew*(4.+xnew)+xold*(4.+xold)-(ynew-yold)*(2+ynew+yold))/2)*(exp(2.+4.*xnew)+exp(0.5+3.*ynew))/(exp(2.+4.*xold)+exp(0.5+3.*yold));
}


void BiGaussEx(int N, fNum x0, fNum y0)
{
	const gsl_rng_type * T;
	gsl_rng * r;
	r = gsl_rng_alloc (gsl_rng_mt19937);
	fNum RandNum = 0.;

	fNum x = x0;
	fNum y = y0;
	fNum xn = 0.;
	fNum yn = 0.;

	fNum MetR = 0.;
	for (int i = 0; i < N; i++)
	{
	RandNum = gsl_ran_gaussian(r, 1.0);
	xn = x + RandNum;
	RandNum = gsl_ran_gaussian(r, 1.0);
	yn = y + RandNum;

	MetR = BiGaussLike2(xn, x, yn, y);
	RandNum = gsl_rng_uniform(r);

		if ( RandNum < MetR)
		{
			x = xn;
			y = yn;
			printf("%lf %lf\n", x, y);
		}
	}

}

*/

