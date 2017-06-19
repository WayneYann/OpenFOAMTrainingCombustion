class linearModelChemistry
{

public:

	linearModelChemistry(OpenSMOKE::ThermodynamicsMap_CHEMKIN& thermoMap_, OpenSMOKE::KineticsMap_CHEMKIN& kineticsMap_)
	{
		NC_ = thermoMap_.NumberOfSpecies();
		NE_ = NC_+1;
		ChangeDimensions(NC_, &omega_, true);
		ChangeDimensions(NC_, &x_, true);
		ChangeDimensions(NC_, &R_, true);
		ChangeDimensions(NC_, &c_, true);
		ChangeDimensions(NE_, &y_plus_, true);
		ChangeDimensions(NE_, &dy_plus_, true);
		ChangeDimensions(NE_, &dy_original_, true);
		Jdiagonal_.resize(NE_);
	}

	void reactionSourceTerms(	OpenSMOKE::ThermodynamicsMap_CHEMKIN& thermoMap_, OpenSMOKE::KineticsMap_CHEMKIN& kineticsMap_,
					const OpenSMOKE::OpenSMOKEVectorDouble& y, const double P0, OpenSMOKE::OpenSMOKEVectorDouble& S);

	void reactionJacobian( 	OpenSMOKE::ThermodynamicsMap_CHEMKIN& thermoMap_, OpenSMOKE::KineticsMap_CHEMKIN& kineticsMap_,
		       		const OpenSMOKE::OpenSMOKEVectorDouble& y, const double P0,
		       		Eigen::VectorXd &J);

private:

	unsigned int NC_;
	unsigned int NE_; 
	OpenSMOKE::OpenSMOKEVectorDouble omega_;
	OpenSMOKE::OpenSMOKEVectorDouble x_;
	OpenSMOKE::OpenSMOKEVectorDouble R_;
	OpenSMOKE::OpenSMOKEVectorDouble c_;

	OpenSMOKE::OpenSMOKEVectorDouble y_plus_;
     	OpenSMOKE::OpenSMOKEVectorDouble dy_plus_;
     	OpenSMOKE::OpenSMOKEVectorDouble dy_original_;

	Eigen::VectorXd Jdiagonal_;
};

void linearModelChemistry::reactionSourceTerms(	OpenSMOKE::ThermodynamicsMap_CHEMKIN& thermoMap_, OpenSMOKE::KineticsMap_CHEMKIN& kineticsMap_,
						const OpenSMOKE::OpenSMOKEVectorDouble& y, const double P0, OpenSMOKE::OpenSMOKEVectorDouble& S)
{
	for(unsigned int i=1;i<=NC_;++i)
		omega_[i] = std::max(y[i], 0.);
	const double T = y[NC_+1];

	// Calculates the pressure and the concentrations of species
	double MW_ = 0.;
	thermoMap_.MoleFractions_From_MassFractions(x_.GetHandle(), MW_, omega_.GetHandle());
	const double cTot_ = P0/PhysicalConstants::R_J_kmol/T;
    	const double rho_ = cTot_*MW_;
	Product(cTot_, x_, &c_);

	// Calculates thermodynamic properties
	thermoMap_.SetTemperature(T);
	thermoMap_.SetPressure(P0);

	// Calculates kinetics
	kineticsMap_.SetTemperature(T);
	kineticsMap_.SetPressure(P0);
	kineticsMap_.KineticConstants();
	kineticsMap_.ReactionRates(c_.GetHandle());
	kineticsMap_.FormationRates(R_.GetHandle());

	// Species
	for (unsigned int i=1;i<=NC_;++i)	
		S[i] = R_[i]*thermoMap_.MW(i-1);

	// Energy
    	{
		const double CpMixMolar = thermoMap_.cpMolar_Mixture_From_MoleFractions(x_.GetHandle());
		const double CpMixMass_ = CpMixMolar / MW_;
		const double QR_ = kineticsMap_.HeatRelease(R_.GetHandle());
		
		S[NC_+1]  = QR_;
	}
}

void linearModelChemistry::reactionJacobian( 	OpenSMOKE::ThermodynamicsMap_CHEMKIN& thermoMap_, OpenSMOKE::KineticsMap_CHEMKIN& kineticsMap_,
		       				const OpenSMOKE::OpenSMOKEVectorDouble& y, const double P0,
		       				Eigen::VectorXd &J ) 
{
     // Calculated as suggested by Buzzi (private communication)
     const double ZERO_DER = std::sqrt(OPENSMOKE_TINY_FLOAT);
     const double ETA2 = std::sqrt(OpenSMOKE::OPENSMOKE_MACH_EPS_DOUBLE);
     const double BETA = 1.e3 * OpenSMOKE::OPENSMOKE_MACH_EPS_DOUBLE;
     const double TOLR = 1.e-7;
     const double TOLA = 1.e-12;

     for(unsigned int i=1;i<=NE_;i++)
		y_plus_[i] = y[i];

     // Call equations
     reactionSourceTerms(thermoMap_, kineticsMap_, y, P0, dy_original_);     

     // Derivatives with respect to y[kd]
     for(int kd=1;kd<=NE_;kd++)
     {
         double hf = 1.e0;
         double error_weight = 1./(TOLA+TOLR*fabs(y[kd]));
         double hJ = ETA2 * fabs(std::max(y[kd], 1./error_weight));
         double hJf = hf/error_weight;
         hJ = std::max(hJ, hJf);
         hJ = std::max(hJ, ZERO_DER);

         double dy = std::min(hJ, 1.e-3 + 1e-3*fabs(y[kd]));
         double udy = 1. / dy;
         y_plus_[kd] += dy;

	 reactionSourceTerms(thermoMap_, kineticsMap_, y_plus_, P0, dy_plus_);

         J(kd-1) = (dy_plus_[kd]-dy_original_[kd]) * udy;

         y_plus_[kd] = y[kd];
     }
}

