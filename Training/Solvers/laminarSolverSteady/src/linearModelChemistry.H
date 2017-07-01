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

class linearModelChemistry
{

public:

	linearModelChemistry(OpenSMOKE::ThermodynamicsMap_CHEMKIN& thermoMap, OpenSMOKE::KineticsMap_CHEMKIN& kineticsMap);

	void reactionSourceTerms( OpenSMOKE::ThermodynamicsMap_CHEMKIN& thermoMap, OpenSMOKE::KineticsMap_CHEMKIN& kineticsMap,
				  const Eigen::VectorXd& y, const double P0, Eigen::VectorXd& S);

	void reactionJacobian( 	OpenSMOKE::ThermodynamicsMap_CHEMKIN& thermoMap, OpenSMOKE::KineticsMap_CHEMKIN& kineticsMap,
		       		const Eigen::VectorXd& y, const double P0,
				Eigen::VectorXd &J );

private:

	unsigned int NC_;
	unsigned int NE_; 

	Eigen::VectorXd omega_;
	Eigen::VectorXd x_;
	Eigen::VectorXd R_;
	Eigen::VectorXd c_;

	Eigen::VectorXd y_plus_;
     	Eigen::VectorXd dy_plus_;
     	Eigen::VectorXd dy_original_;
};

linearModelChemistry::linearModelChemistry(OpenSMOKE::ThermodynamicsMap_CHEMKIN& thermoMap, OpenSMOKE::KineticsMap_CHEMKIN& kineticsMap)
{
	NC_ = thermoMap.NumberOfSpecies();
	NE_ = NC_+1;
	omega_.resize(NC_);
	x_.resize(NC_);
	R_.resize(NC_);
	c_.resize(NC_);
	y_plus_.resize(NE_);
	dy_plus_.resize(NE_);
	dy_original_.resize(NE_);
}

void linearModelChemistry::reactionSourceTerms(	OpenSMOKE::ThermodynamicsMap_CHEMKIN& thermoMap, OpenSMOKE::KineticsMap_CHEMKIN& kineticsMap,
						const Eigen::VectorXd& y, const double P0, Eigen::VectorXd& S)
{
	for(unsigned int i=0;i<NC_;++i)
		omega_(i) = std::max(y(i), 0.);
	const double T = y(NC_);

	// Calculates the pressure and the concentrations of species
	double MW = 0.;
	thermoMap.MoleFractions_From_MassFractions(x_.data(), MW, omega_.data());
	const double cTot = P0/PhysicalConstants::R_J_kmol/T;
	c_ = cTot*x_;

	// Calculates thermodynamic properties
	thermoMap.SetTemperature(T);
	thermoMap.SetPressure(P0);

	// Calculates kinetics
	kineticsMap.SetTemperature(T);
	kineticsMap.SetPressure(P0);
	kineticsMap.KineticConstants();
	kineticsMap.ReactionRates(c_.data());
	kineticsMap.FormationRates(R_.data());

	// Species
	for (unsigned int i=0;i<NC_;++i)	
		S(i) = R_(i)*thermoMap.MW(i);

	// Energy
	const double QR_ = kineticsMap.HeatRelease(R_.data());
	S(NC_) = QR_;
}

void linearModelChemistry::reactionJacobian( 	OpenSMOKE::ThermodynamicsMap_CHEMKIN& thermoMap, OpenSMOKE::KineticsMap_CHEMKIN& kineticsMap,
		       				const Eigen::VectorXd& y, const double P0,
						Eigen::VectorXd &J ) 
{
     // Calculated as suggested by Buzzi (private communication)
     const double ZERO_DER = std::sqrt(OPENSMOKE_TINY_FLOAT);
     const double ETA2 = std::sqrt(OpenSMOKE::OPENSMOKE_MACH_EPS_DOUBLE);
     const double TOLR = 1.e-7;
     const double TOLA = 1.e-12;

     y_plus_ = y;

     // Call equations
     reactionSourceTerms(thermoMap, kineticsMap, y, P0, dy_original_);     

     // Derivatives with respect to y(kd)
     for(int kd=0;kd<NE_;kd++)
     {
         double hf = 1.e0;
         double error_weight = 1./(TOLA+TOLR*fabs(y(kd)));
         double hJ = ETA2 * fabs(std::max(y(kd), 1./error_weight));
         double hJf = hf/error_weight;
         hJ = std::max(hJ, hJf);
         hJ = std::max(hJ, ZERO_DER);

         double dy = std::min(hJ, 1.e-3 + 1e-3*fabs(y(kd)));
         double udy = 1. / dy;
         y_plus_(kd) += dy;

	 reactionSourceTerms(thermoMap, kineticsMap, y_plus_, P0, dy_plus_);

         J(kd) = (dy_plus_(kd)-dy_original_(kd)) * udy;

         y_plus_(kd) = y(kd);
     }
}

