#include <iostream>
#include <vector>
#include <math.h>
#include <memory.h>
#include <pngwriter.h>

using namespace std;

#define access_mtx_value(mtx, x, y, xsize) mtx[y*xsize + x]

#ifndef BOOL
#define BOOL int
#endif

#ifndef TRUE
#define TRUE -1
#endif

#ifndef FALSE
#define FALSE 0
#endif

// helper functions

double square(double fValue);

inline double square(double fValue)
{
	return fValue*fValue;
}



//! \brief This class represents a gaussian chirplet
class CGaussianChirplet
{
public:
	//! \brief Contructor
	CGaussianChirplet();

	//! Calculate the value of the chirplet at time fTime
	//!
	//! \param[in] fTime the time to calculate the chirplet value
	//! \return the chriplet value at time fTime
	double value(double fTime) const;

	//! this will be used to render the 2d image freq x time
	//! still not implemented
	double value2d(double fTime, double fFrequency) const;

public:
	double m_fAmplitude;

	// gaussian parameters {{

	//! center of the gaussian
	double m_fCentralTime;

	//! sigma, the width of the gaussian
	double m_fSigma;

	// }} gaussian parameters

	// the chirplet/oscilatory parameters {{

	//! the radial frequency
	double m_fw;

	//! the phase of the oscilatory part
	double m_fPhase;

	//! the chrirp factor
	double m_fc;

	// }} the chirplet parameters
};

typedef vector<CGaussianChirplet> CGaussianChirpletVector;

CGaussianChirplet::CGaussianChirplet()
{

}

double CGaussianChirplet::value(double fTime) const
{
	return m_fAmplitude*exp(-square((fTime - m_fCentralTime)/m_fSigma))*cos((m_fw + m_fc*fTime)*fTime + m_fPhase);
}

double CGaussianChirplet::value2d(double fTime, double fFrequency) const
{
	return m_fAmplitude*exp(-(square((fTime - m_fCentralTime)/m_fSigma) + square((fFrequency - (m_fw + m_fc*((fTime - m_fCentralTime))))*m_fSigma)));
}

void render_image(unsigned int uiTimeSamples, unsigned int uiFrequencySamples, const double* pfImage, const char *pcszOutFile, BOOL bInverse = TRUE);

void render_image(unsigned int uiTimeSamples, unsigned int uiFrequencySamples, const double* pfImage, const char *pcszOutFile, BOOL bInverse)
{
	pngwriter pngImage((int)uiTimeSamples, (int)uiFrequencySamples, 1.0, pcszOutFile);

	unsigned int i, j;

	for(i = 0; i < uiTimeSamples; i++)
	{
		for(j = 0; j < uiFrequencySamples; j++)
		{
			double fValue = access_mtx_value(pfImage, i, j, uiTimeSamples);

			// secure a max of 1.0 for value
			fValue = fValue <= 1.0 ? fValue : 1.0;

			if(bInverse)
				fValue = 1.0 - fValue;

			pngImage.plot((int)i, (int)j, fValue, fValue, fValue);

		}
	}

	pngImage.close();
}

//! Scales a double vector so that the maximum value is set to 1.0
void scale_to_max(double* pfValues, unsigned int nSize);

void scale_to_max(double* pfValues, unsigned int nSize)
{
	unsigned int i;
	double fBuffer = 0.0;

	// find the max value
	for(i = 0; i < nSize; i++)
	{
		if(pfValues[i] > fBuffer)
		{
			fBuffer = pfValues[i];
		}
	}

	if(fBuffer != 0.0)
	{
		fBuffer = 1.0/fBuffer;

		for(i = 0; i < nSize; i++)
		{
			pfValues[i] *= fBuffer;
		}
	}
}

//! this scales bla bla bla, implement later
double gaussian_adjust_intensity(double fValue)
{
	return exp(-square((fValue-0.5)*10.0));
}

int main(int argc, const char *argv[])
{
    if(argc < 1)
		return 0;

    unsigned int i, j, k;
    unsigned int uiTimeSamples = 800;
    unsigned int uiFrequencySamples = 480;
    unsigned int uiMaxChirplets = 1024;
    BOOL bUseAmplitude = TRUE;
    BOOL bGaussianAdjustIntensity = FALSE;
    BOOL bInverse = TRUE;
    BOOL bScale = TRUE;
    double fSamplingFrequency = 128.0;
    double fTimeMin = 0.0;
    double fTimeMax = 16.0;
	double fFreqMin = 0.0;
    double fFreqMax = 32.0;
    FILE* pInFile;

	CGaussianChirpletVector gcvector;

	pInFile = fopen(argv[1], "rt");

	if(!pInFile)
	{
		fprintf(stdout, "Impossible to open input file.\n");
	}

	// read the file {{

	CGaussianChirpletVector::iterator it;
	CGaussianChirplet gchirp;

	float frbuffer[4];
	int nrbuffer[2];

	// while(fscanf(pInFile, "%e %e %e %e %e %e", &gchirp.m_fAmplitude, &gchirp.m_fSigma, &gchirp.m_fCentralTime, &gchirp.m_fw, &dummy[0], &gchirp.m_fc )!=EOF)
	for(i = 0; fscanf(pInFile, "%e %d %d %e %e %e", &frbuffer[0], &nrbuffer[0], &nrbuffer[1], &frbuffer[1], &frbuffer[2], &frbuffer[3] )!=EOF; i++)
	{
		gchirp.m_fSigma = ((double)nrbuffer[0])/fSamplingFrequency;
		gchirp.m_fCentralTime = ((double)nrbuffer[1])/fSamplingFrequency; // is sampling 128? probably

		if(bUseAmplitude)
		{
			// gchirp.m_fAmplitude = sqrt(frbuffer[0]); // lets pretend everybody is equal for tests
			gchirp.m_fAmplitude = sqrt(frbuffer[0]); // lets pretend everybody is equal for tests
		}
		else
		{
			gchirp.m_fAmplitude = 1.0; // lets pretend everybody is equal for tests
		}

		gchirp.m_fw = fSamplingFrequency*frbuffer[1]/(2*M_PI);
		gchirp.m_fPhase = frbuffer[2];
		gchirp.m_fc = fSamplingFrequency*fSamplingFrequency*frbuffer[3]/(2*M_PI);

		/*
		CGaussianChirplet test_chirp;

		test_chirp.m_fAmplitude = 1.0;
		test_chirp.m_fw = 32.0;
		test_chirp.m_fc = 10;
		test_chirp.m_fCentralTime = 4.0;
		test_chirp.m_fSigma = 1.0;
		*/

		gcvector.push_back(gchirp);

		i++;

		if(i >= uiMaxChirplets)
			break;
	}

	fclose(pInFile);

	// }} read the file

	// allocate the double bitmap
	double *pfImage = new double[uiTimeSamples*uiFrequencySamples];

    // reset the bitmap values
    memset(pfImage, 0, uiTimeSamples*uiFrequencySamples*sizeof(double));

    // steps in time
    double fTimeStep = (fTimeMax - fTimeMin)/((double)uiTimeSamples);

    // steps in frequency
    double fFreqStep = (fFreqMax - fFreqMin)/((double)uiFrequencySamples);

	// gaussian chirplet loop
	// render to the float bitmap
    for(k = 0; k < gcvector.size(); k++)
	{
		// time loop
		for(i = 0; i < uiTimeSamples; i++)
		{
			// frequency loop
			for(j = 0; j < uiFrequencySamples; j++)
			{
				double fTime = fTimeMin + ((double)i)*fTimeStep;
				double fFreq = fFreqMin + ((double)j)*fFreqStep;

				// fTime =
				// fFreq =
				double fValue = gcvector[k].value2d(fTime, fFreq);
				double fScaledValue = fValue / gcvector[k].m_fAmplitude;

				// interesting cosmetic, just draw values around 0.5 of the max
				// if(fNormalValue >= 0.4 && fNormalValue <= 0.6)

				// access_mtx_value(pfImage, i, j, uiTimeSamples) += fValue;
				access_mtx_value(pfImage, i, j, uiTimeSamples) += bGaussianAdjustIntensity ? (bUseAmplitude ? gcvector[k].m_fAmplitude : 1.0) * gaussian_adjust_intensity(fScaledValue) : fValue;
			}
		}
	}

	// scale
	if(bScale)
	{
		scale_to_max(pfImage, uiTimeSamples*uiFrequencySamples);
	}

	render_image(uiTimeSamples, uiFrequencySamples, pfImage, argv[2], bInverse);

	// TODO: here you need to write the bitmap

    return 0;
}
