#include <iostream>
#include <vector>
#include <math.h>
#include <memory.h>
#include <pngwriter.h>

using namespace std;

#define access_mtx_value(mtx, x, y, xsize) mtx[y*xsize + x]

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

void render_image(unsigned int nTimeSamples, unsigned int nFrequencySamples, const double* pfImage, const char *pcszOutFile)
{
	pngwriter pngImage((int)nTimeSamples, (int)nFrequencySamples, 1.0, pcszOutFile);

	unsigned int i, j;

	for(i = 0; i < nTimeSamples; i++)
	{
		for(j = 0; j < nFrequencySamples; j++)
		{
			double fValue = access_mtx_value(pfImage, i, j, nTimeSamples);
			pngImage.plot((int)i, (int)j, fValue, fValue, fValue);
		}
	}

	pngImage.close();
}

void normalize(double* pfValues, unsigned int nSize)
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

int main(int argc, const char *argv[])
{
    // cout << "Hello world!" << endl;

    if(argc < 1)
		return 0;

    unsigned int i, j, k;
    unsigned int nTimeSamples = 800;
    unsigned int nFrequencySamples = 480;
    // unsigned int nSamplesShift = 1;
    double fTimeMin = 0.0;
    double fTimeMax = 8.0;
	double fFreqMin = 0.0;
    double fFreqMax = 32.0;

	CGaussianChirpletVector gcvector;

	FILE* pInFile = fopen(argv[1], "rt");

	if(!pInFile)
	{
		fprintf(stdout, "Impossible to open input file.\n");
	}

	CGaussianChirpletVector::iterator it;

	// it = gcvector.begin();

	CGaussianChirplet gchirp;

	float frbuffer[4];
	int nrbuffer[2];

	// while(fscanf(pInFile, "%e %e %e %e %e %e", &gchirp.m_fAmplitude, &gchirp.m_fSigma, &gchirp.m_fCentralTime, &gchirp.m_fw, &dummy[0], &gchirp.m_fc )!=EOF)
	while(fscanf(pInFile, "%e %d %d %e %e %e", &frbuffer[0], &nrbuffer[0], &nrbuffer[1], &frbuffer[1], &frbuffer[2], &frbuffer[3] )!=EOF)
	{
		gchirp.m_fSigma = ((double)nrbuffer[0])/128.0;;
		gchirp.m_fCentralTime = ((double)nrbuffer[1])/128.0; // is sampling 128? probably
		gchirp.m_fAmplitude = 1.0; // lets pretend everybody is equal for tests
		gchirp.m_fw = 128*frbuffer[1]/6.8;
		gchirp.m_fc = 128*128*frbuffer[3]/6.8;

		/*
		CGaussianChirplet test_chirp;

		test_chirp.m_fAmplitude = 1.0;
		test_chirp.m_fw = 32.0;
		test_chirp.m_fc = 10;
		test_chirp.m_fCentralTime = 4.0;
		test_chirp.m_fSigma = 1.0;
		*/

		gcvector.push_back(gchirp);
	}

	fclose(pInFile);

	// TODO: here you need to read the atoms and fill gcvector

    double *pfImage = new double[nTimeSamples*nFrequencySamples];

    // reset the bitmap values
    memset(pfImage, 0, nTimeSamples*nFrequencySamples*sizeof(double));

    // steps in time
    double fTimeStep = (fTimeMax - fTimeMin)/((double)nTimeSamples);

    // steps in frequency
    double fFreqStep = (fFreqMax - fFreqMin)/((double)nFrequencySamples);

	// gaussian chirplet loop
    for(k = 0; k < gcvector.size(); k++)
	{
		// time loop
		for(i = 0; i < nTimeSamples; i++)
		{
			// frequency loop
			for(j = 0; j < nFrequencySamples; j++)
			{
				double fTime = ((double)i)*fTimeStep;
				double fFreq = ((double)j)*fFreqStep;

				// fTime =
				// fFreq =

				access_mtx_value(pfImage, i, j, nTimeSamples) += gcvector[k].value2d(fTime, fFreq);
			}
		}
	}

	normalize(pfImage, nTimeSamples*nFrequencySamples);

	render_image(nTimeSamples, nFrequencySamples, pfImage, argv[2]);

	// TODO: here you need to write the bitmap

    return 0;
}
