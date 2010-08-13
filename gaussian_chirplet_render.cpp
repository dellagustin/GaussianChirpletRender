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

double CGaussianChirplet::value(double fTime) const
{
	return m_fAmplitude*exp(-square((fTime - m_fCentralTime)/m_fSigma))*cos((m_fw + m_fc*fTime)*fTime + m_fPhase);
}

double CGaussianChirplet::value2d(double fTime, double fFrequency) const
{
	return 0;
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
	int i;
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

int main()
{
    // cout << "Hello world!" << endl;
    unsigned int i, j, k;
    unsigned int nTimeSamples = 800;
    unsigned int nFrequencySamples = 480;
    // unsigned int nSamplesShift = 1;
    double fTimeMin = 0.0;
    double fTimeMax = 8.0;
	double fFreqMin = 2.0;
    double fFreqMax = 64.0;

	CGaussianChirpletVector gcvector;

	CGaussianChirpletVector::iterator it;

	// it = gcvector.begin();

	// gcvector.insert( ,it);


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

				access_mtx_value(pfImage, i, j, nTimeSamples) = gcvector[i].value2d(fTime, fFreq);
			}
		}
	}

	normalize(pfImage, nTimeSamples*nFrequencySamples);

	render_image(nTimeSamples, nFrequencySamples, pfImage, "out.png");

	// TODO: here you need to write the bitmap

    return 0;
}
