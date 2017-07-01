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

	void setTemperature(const double TStart);
	void setPressure(const double P);
	void setEnthalpy(const double Hfixed);

	label nEqns() const;

	void derivatives( const scalar t, const scalarField& cc, scalarField& dcdt ) const;

	void jacobian( const scalar t, const scalarField& cc, scalarField& dfdt, scalarSquareMatrix& dfdc ) const;

private:

	OpenSMOKE::ThermodynamicsMap_CHEMKIN& 	thermoMap_;
	OpenSMOKE::KineticsMap_CHEMKIN&		kineticsMap_;
	
	double TStart_;		// initial temperature (in K)
	double P0_;		// pressure, constant (in Pa)
	double Hfixed_;		// mass enthalpy, constant (in J/kg)
};

batchAdiabaticOdeSystem::batchAdiabaticOdeSystem(OpenSMOKE::ThermodynamicsMap_CHEMKIN& thermoMap, OpenSMOKE::KineticsMap_CHEMKIN& kineticsMap) : 
	ODESystem(),
	thermoMap_(thermoMap),
	kineticsMap_(kineticsMap)
{
}

void batchAdiabaticOdeSystem::setTemperature(const double TStart)
{
	TStart_ = TStart;
}

void batchAdiabaticOdeSystem::setPressure(const double P0)
{
	P0_ = P0;
}

void batchAdiabaticOdeSystem::setEnthalpy(const double Hfixed)
{
	Hfixed_ = Hfixed;
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
	
	// Calculate temperature (pressure is assumed constant)
	const double T_ = thermoMap_.GetTemperatureFromEnthalpyAndMoleFractions(Hfixed_*mw, P0_, x.data(), TStart_);
	
	// Calculates thermodynamic properties
	thermoMap_.SetTemperature(T_);
	thermoMap_.SetPressure(P0_);
	
	// Calculates kinetics
	Eigen::VectorXd R(thermoMap_.NumberOfSpecies());
	kineticsMap_.SetTemperature(T_);
	kineticsMap_.SetPressure(P0_);
	kineticsMap_.ReactionRates(c.data());
	kineticsMap_.FormationRates(R.data());
	
	// Species equations
	for (unsigned int i=0;i<thermoMap_.NumberOfSpecies();i++)	
		dcdt[i] = R(i);
}

void batchAdiabaticOdeSystem::jacobian( const scalar t, const scalarField& cc, scalarField& dfdt, scalarSquareMatrix& dfdc ) const
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
	
	// Calculate temperature (pressure is assumed constant)
	const double T_ = thermoMap_.GetTemperatureFromEnthalpyAndMoleFractions(Hfixed_*mw, P0_, x.data(), TStart_);

	// Set maps
	thermoMap_.SetTemperature(T_);
	thermoMap_.SetPressure(P0_);
	kineticsMap_.SetTemperature(T_);
	kineticsMap_.SetPressure(P0_);

	// Derivative of formation rates with respect to concentrations (in 1/s)
	Eigen::MatrixXd dR_over_dC(thermoMap_.NumberOfSpecies(), thermoMap_.NumberOfSpecies());
	kineticsMap_.DerivativesOfFormationRates(c.data(), &dR_over_dC);

	for (unsigned int i=0;i<thermoMap_.NumberOfSpecies();i++)
		dfdt[i] = 0.;	

	for (unsigned int i=0;i<thermoMap_.NumberOfSpecies();i++)
		for (unsigned int j=0;j<thermoMap_.NumberOfSpecies();j++)
			dfdc[i][j] = dR_over_dC(i,j);
}

