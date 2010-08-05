#include <iostream>
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

public
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
}

double CGaussianChirplet::value(double fTime) const
{
	return m_fAmplitude*exp(-square((fTime - m_fCentralTime)/m_fSigma))*cos((m_fw + m_fc*fTime)*fTime + m_fPhase);
}

int main()
{
    cout << "Hello world!" << endl;
    return 0;
}
