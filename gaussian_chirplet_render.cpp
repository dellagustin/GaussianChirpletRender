#include <iostream>
#include <vector>
#include <math.h>

using namespace std;

// helper functions

double square(double fValue);

inline double square(double fValue)
{
	return fValue*fValue;
}

class CGaussianChirplet
{
public:
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
    return 0;
}
