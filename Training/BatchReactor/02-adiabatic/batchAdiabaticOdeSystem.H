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

class batchAdiabaticOdeSystem : public ODESystem
{

public:

	batchAdiabaticOdeSystem(OpenSMOKE::ThermodynamicsMap_CHEMKIN& thermoMap, OpenSMOKE::KineticsMap_CHEMKIN& kineticsMap);

	void setInitialTemperature(const double TInitial);
	void setInitialPressure(const double PInitial);
	void setEnthalpy(const double H);

	label nEqns() const;

	void derivatives( const scalar t, const scalarField& c, scalarField& dcdt ) const;

	void jacobian( const scalar t, const scalarField& c, scalarField& dfdt, scalarSquareMatrix& dfdc ) const;

private:

	OpenSMOKE::ThermodynamicsMap_CHEMKIN& 	thermoMap_;
	OpenSMOKE::KineticsMap_CHEMKIN&		kineticsMap_;
	
	double TInitial_;
	double PInitial_;
	double H_;
};

batchAdiabaticOdeSystem::batchAdiabaticOdeSystem(OpenSMOKE::ThermodynamicsMap_CHEMKIN& thermoMap, OpenSMOKE::KineticsMap_CHEMKIN& kineticsMap) : 
	ODESystem(),
	thermoMap_(thermoMap),
	kineticsMap_(kineticsMap)
{}

void batchAdiabaticOdeSystem::setInitialTemperature(const double TInitial)
{
	TInitial_ = TInitial;
}

void batchAdiabaticOdeSystem::setInitialPressure(const double PInitial)
{
	PInitial_ = PInitial;
}

void batchAdiabaticOdeSystem::setEnthalpy(const double H)
{
	H_ = H;
}

label batchAdiabaticOdeSystem::nEqns() const
{
	return thermoMap_.NumberOfSpecies();
}

void batchAdiabaticOdeSystem::derivatives( const scalar t, const scalarField& c, scalarField& dcdt ) const
{
	// Calculate temperature from enthalpy
	Eigen::VectorXd x(thermoMap_.NumberOfSpecies());
	const double cTot = std::accumulate(c.begin(), c.end(), 0.);
	for (unsigned int i=0;i<thermoMap_.NumberOfSpecies();i++)
		x(i) = c[i]/cTot;
	const double T = thermoMap_.GetTemperatureFromEnthalpyAndMoleFractions(H_, PInitial_, x.data(), TInitial_);
	const double P = cTot*(PhysicalConstants::R_J_kmol*T);	

	// Calculates thermodynamic properties
	thermoMap_.SetTemperature(T);
	thermoMap_.SetPressure(P);
	
	// Calculates kinetics
	Eigen::VectorXd R(thermoMap_.NumberOfSpecies());
	kineticsMap_.SetTemperature(T);
	kineticsMap_.SetPressure(P);
	kineticsMap_.ReactionRates(c.cdata());
	kineticsMap_.FormationRates(R.data());
	
	// Species equations
	for (unsigned int i=0;i<thermoMap_.NumberOfSpecies();i++)	
		dcdt[i] = R(i);
}

void batchAdiabaticOdeSystem::jacobian( const scalar t, const scalarField& c, scalarField& dfdt, scalarSquareMatrix& dfdc ) const
{
	// Calculate temperature from enthalpy
	Eigen::VectorXd x(thermoMap_.NumberOfSpecies());
	const double cTot = std::accumulate(c.begin(), c.end(), 0.);
	for (unsigned int i=0;i<thermoMap_.NumberOfSpecies();i++)
		x(i) = c[i]/cTot;
	const double T = thermoMap_.GetTemperatureFromEnthalpyAndMoleFractions(H_, PInitial_, x.data(), TInitial_);
	const double P = cTot*(PhysicalConstants::R_J_kmol*T);	

	// Set maps
	thermoMap_.SetTemperature(T);
	thermoMap_.SetPressure(P);
	kineticsMap_.SetTemperature(T);
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

