#include <iostream>
#include <vector>
#include <math.h>
#include <memory.h>

using namespace std;

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

	// TODO: here you need to read the atoms and fill gcvector

    double *pfBitmap = new double[nTimeSamples*nFrequencySamples];

    // reset the bitmap values
    memset(pfBitmap, 0, nTimeSamples*nFrequencySamples*sizeof(double));

    // steps in time
    double fTimeStep = (fTimeMax - fTimeMin)/((double)nTimeSamples);

    // steps in frequency
    double fFreqStep = (fTimeMax - fTimeMin)/((double)nTimeFrequency);

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

				pfBitmap[j*nFrequencySamples + i] = gcvector[i].value2d(fTime, fFreq);
			}
		}
	}

	// TODO: here you need to write the bitmap

    return 0;
}
