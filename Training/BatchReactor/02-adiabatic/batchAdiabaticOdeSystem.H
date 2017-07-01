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
	void setInternalEnergy(const double U);

	label nEqns() const;

	void derivatives( const scalar t, const scalarField& c, scalarField& dcdt ) const;

	void jacobian( const scalar t, const scalarField& c, scalarField& dfdt, scalarSquareMatrix& dfdc ) const;

private:

	OpenSMOKE::ThermodynamicsMap_CHEMKIN& 	thermoMap_;
	OpenSMOKE::KineticsMap_CHEMKIN&		kineticsMap_;
	
	double TInitial_;	// initial temperature (in K)
	double PInitial_;	// initial pressure (in Pa)
	double U_;		// mass internal energy, constant (in J/kg)

	int maxIterations_;	// max number iterations for estimating the pressure
};

batchAdiabaticOdeSystem::batchAdiabaticOdeSystem(OpenSMOKE::ThermodynamicsMap_CHEMKIN& thermoMap, OpenSMOKE::KineticsMap_CHEMKIN& kineticsMap) : 
	ODESystem(),
	thermoMap_(thermoMap),
	kineticsMap_(kineticsMap)
{
	maxIterations_ = 10;
}

void batchAdiabaticOdeSystem::setInitialTemperature(const double TInitial)
{
	TInitial_ = TInitial;
}

void batchAdiabaticOdeSystem::setInitialPressure(const double PInitial)
{
	PInitial_ = PInitial;
}

void batchAdiabaticOdeSystem::setInternalEnergy(const double U)
{
	U_ = U;
}

label batchAdiabaticOdeSystem::nEqns() const
{
	return thermoMap_.NumberOfSpecies();
}

void batchAdiabaticOdeSystem::derivatives( const scalar t, const scalarField& cc, scalarField& dcdt ) const
{
	// Reconstruct concentrations and mole fractions
	Eigen::VectorXd c(thermoMap_.NumberOfSpecies());
	for (unsigned int i=0;i<thermoMap_.NumberOfSpecies();i++)
		c(i) = std::max(cc[i],0.);
	Eigen::VectorXd x(thermoMap_.NumberOfSpecies());
	const double cTot = c.sum();
	for (unsigned int i=0;i<thermoMap_.NumberOfSpecies();i++)
		x(i) = c[i]/cTot;
	const double mw = thermoMap_.MolecularWeight_From_MoleFractions(x.data());
	
	// Calculate temperature and pressure (by successive substitutions)
	double P = PInitial_;
	double T = TInitial_;
	for(int i=0;i<maxIterations_;i++)
	{	
		const double Pold = P;	
		const double H = U_+P/(cTot*mw);
		T = thermoMap_.GetTemperatureFromEnthalpyAndMoleFractions(H*mw, P, x.data(), T);
		P = cTot*(PhysicalConstants::R_J_kmol*T);
		if (std::fabs(P-Pold)/P<1e-4)	break;	
	}

	// Calculates thermodynamic properties
	thermoMap_.SetTemperature(T);
	thermoMap_.SetPressure(P);
	
	// Calculates kinetics
	Eigen::VectorXd R(thermoMap_.NumberOfSpecies());
	kineticsMap_.SetTemperature(T);
	kineticsMap_.SetPressure(P);
	kineticsMap_.ReactionRates(c.data());
	kineticsMap_.FormationRates(R.data());
	
	// Species equations
	for (unsigned int i=0;i<thermoMap_.NumberOfSpecies();i++)	
		dcdt[i] = R(i);
}

void batchAdiabaticOdeSystem::jacobian( const scalar t, const scalarField& cc, scalarField& dfdt, scalarSquareMatrix& dfdc ) const
{
	// Reconstruct concentrations and mole fractions
	Eigen::VectorXd c(thermoMap_.NumberOfSpecies());	// concentrations of species (in kmol/m3)
	for (unsigned int i=0;i<thermoMap_.NumberOfSpecies();i++)
		c(i) = std::max(cc[i],0.);
	Eigen::VectorXd x(thermoMap_.NumberOfSpecies());
	const double cTot = c.sum();
	for (unsigned int i=0;i<thermoMap_.NumberOfSpecies();i++)
		x(i) = c[i]/cTot;
	const double mw = thermoMap_.MolecularWeight_From_MoleFractions(x.data());

	// Calculate temperature and pressure (by successive substitutions)
	double P = PInitial_;
	double T = TInitial_;
	for(int i=0;i<maxIterations_;i++)
	{	
		const double Pold = P;	
		const double H = U_+P/(cTot*mw);
		T = thermoMap_.GetTemperatureFromEnthalpyAndMoleFractions(H*mw, P, x.data(), T);
		P = cTot*(PhysicalConstants::R_J_kmol*T);
		if (std::fabs(P-Pold)/P<1e-4)	break;	
	}

	// Set maps
	thermoMap_.SetTemperature(T);
	thermoMap_.SetPressure(P);
	kineticsMap_.SetTemperature(T);
	kineticsMap_.SetPressure(P);

	// Derivative of formation rates with respect to concentrations (in 1/s)
	Eigen::MatrixXd dR_over_dC(thermoMap_.NumberOfSpecies(), thermoMap_.NumberOfSpecies());
	kineticsMap_.DerivativesOfFormationRates(c.data(), &dR_over_dC);

	for (unsigned int i=0;i<thermoMap_.NumberOfSpecies();i++)
		dfdt[i] = 0.;	

	for (unsigned int i=0;i<thermoMap_.NumberOfSpecies();i++)
		for (unsigned int j=0;j<thermoMap_.NumberOfSpecies();j++)
			dfdc[i][j] = dR_over_dC(i,j);
}

