/*-----------------------------------------------------------------------*\
|    ___                   ____  __  __  ___  _  _______                  |
|   / _ \ _ __   ___ _ __ / ___||  \/  |/ _ \| |/ / ____| _     _         |
|  | | | | '_ \ / _ \ '_ \\___ \| |\/| | | | | ' /|  _| _| |_ _| |_       |
|  | |_| | |_) |  __/ | | |___) | |  | | |_| | . \| |__|_   _|_   _|      |
|   \___/| .__/ \___|_| |_|____/|_|  |_|\___/|_|\_\_____||_|   |_|        |
|        |_|                                                              |
|                                                                         |
|   Author: Alberto Cuoci <alberto.cuoci@polimi.it>                       |
|   CRECK Modeling Group <http://creckmodeling.chem.polimi.it>            |
|   Department of Chemistry, Materials and Chemical Engineering           |
|   Politecnico di Milano                                                 |
|   P.zza Leonardo da Vinci 32, 20133 Milano                              |
|                                                                         |
\*-----------------------------------------------------------------------*/

// OpenSMOKE++
#include "OpenSMOKEpp"
#include "maps/Maps_CHEMKIN"

// OpenFOAM Ode Solvers
#include "fvCFD.H"
#include "ODESolver.H"

class testOdeSystem : public ODESystem
{

public:

	testOdeSystem(OpenSMOKE::ThermodynamicsMap_CHEMKIN& thermoMap_, OpenSMOKE::KineticsMap_CHEMKIN& kineticsMap_, const double T);

	label nEqns() const;

	void derivatives( const scalar t, const scalarField& c, scalarField& dcdt ) const;

	void jacobian( const scalar t, const scalarField& c, scalarField& dfdt, scalarSquareMatrix& dfdc ) const;

private:

	OpenSMOKE::ThermodynamicsMap_CHEMKIN& 	thermoMap_;
	OpenSMOKE::KineticsMap_CHEMKIN&		kineticsMap_;
	
	double T_;
};

testOdeSystem::testOdeSystem(OpenSMOKE::ThermodynamicsMap_CHEMKIN& thermoMap, OpenSMOKE::KineticsMap_CHEMKIN& kineticsMap, const double T) : 
	ODESystem(),
	thermoMap_(thermoMap),
	kineticsMap_(kineticsMap),
	T_(T)
{}

label testOdeSystem::nEqns() const
{
	return thermoMap_.NumberOfSpecies();
}

void testOdeSystem::derivatives( const scalar t, const scalarField& c, scalarField& dcdt ) const
{
	// Calculate pressure
	const double cTot = std::accumulate(c.begin(), c.end(), 0.);
	const double P = cTot*(PhysicalConstants::R_J_kmol*T_);	

	// Set maps
	thermoMap_.SetTemperature(T_);
	thermoMap_.SetPressure(P);
	kineticsMap_.SetTemperature(T_);
	kineticsMap_.SetPressure(P);
	
	// Calculates kinetics
	Eigen::VectorXd R(thermoMap_.NumberOfSpecies());
	kineticsMap_.ReactionRates(c.cdata());
	kineticsMap_.FormationRates(R.data());
	
	// Species equations
	for (unsigned int i=0;i<thermoMap_.NumberOfSpecies();i++)	
		dcdt[i] = R(i);
}

void testOdeSystem::jacobian( const scalar t, const scalarField& c, scalarField& dfdt, scalarSquareMatrix& dfdc ) const
{
	// Calculate pressure
	const double cTot = std::accumulate(c.begin(), c.end(), 0.);
	const double P = cTot*(PhysicalConstants::R_J_kmol*T_);	

	// Set maps
	thermoMap_.SetTemperature(T_);
	thermoMap_.SetPressure(P);
	kineticsMap_.SetTemperature(T_);
	kineticsMap_.SetPressure(P);

	// Derivative of formation rates with respect to concentrations (in 1/s)
	Eigen::MatrixXd dR_over_dC(thermoMap_.NumberOfSpecies(), thermoMap_.NumberOfSpecies());
	kineticsMap_.DerivativesOfFormationRates(c.cdata(), &dR_over_dC);

	for (unsigned int i=0;i<thermoMap_.NumberOfSpecies();i++)
		dfdt[i] = 0.;	

	for (unsigned int i=0;i<thermoMap_.NumberOfSpecies();i++)
		for (unsigned int j=0;j<thermoMap_.NumberOfSpecies();j++)
			dfdc[i][j] = dR_over_dC(i,j);
}

