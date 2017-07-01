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

class sparkModel
{
public: 

	sparkModel()
	{
		active_ = false;
		ignitionTime_ = 0.;
		duration_ = 0.;
		ignitionTemperature_ = 0.;
		ignitionDiameter_ = 0.;
		nIgnCells_ = 0;
	}
	
	bool active() const { return active_; }
	Foam::vector position() const { return position_; }
	scalar ignitionTime() const { return ignitionTime_; }
	scalar duration() const { return duration_; }
	scalar ignitionTemperature() const { return ignitionTemperature_; }
	scalar ignitionDiameter() const { return ignitionDiameter_; }
	
	void Read(const dictionary& sparkDictionary)
	{
		active_ = Switch(sparkDictionary.lookup(word("spark")));
		
		if(active_ == true)
		{
			position_ = sparkDictionary.lookup("position");
			ignitionTime_ = readScalar(sparkDictionary.lookup("time"));
			duration_ = readScalar(sparkDictionary.lookup("duration"));
			ignitionTemperature_ = readScalar(sparkDictionary.lookup("temperature"));
			ignitionDiameter_ = readScalar(sparkDictionary.lookup("diameter"));   
		}
	}
	
	void Setup(const fvMesh& mesh)
	{
		if (active_ == true)
		{
			const volVectorField& centres = mesh.C();

			label ignitionCell;
			ignitionCell = mesh.findCell(position_);
			
			if (ignitionCell == -1)
			{
				Info << "Warning: theorical ignition position does not correspond to any cell." << endl;
			}
			else 
			{
				Info << "Ignition cell center: " << centres[ignitionCell] << ". Original position: " << position_ <<endl;

				scalar radius = ignitionDiameter_/2.0;
				ignCellsList_.setSize(1);
				ignCellsList_[0] = ignitionCell;

				forAll(centres, celli)
				{
					scalar dist = mag(centres[celli] - position_);

					if (dist < radius && celli != ignitionCell)
					{
						ignCellsList_.setSize(nIgnCells_+1);
						ignCellsList_[nIgnCells_] = celli;
						nIgnCells_++;
					}
				}

				if (ignCellsList_.size() > 0) 
					Info << "Added " << nIgnCells_ << " cells to ignition point" << endl;
				else
					Info << "The solver wasn't able to find more than one cell in the ignition point" <<endl;
			}
		}
	}
	
	void SetSpark(const double t, volScalarField& T)
	{
		if( (active_ == true) && (t >= ignitionTime_) && (t <= (ignitionTime_+duration_)) ) 
			for(int i = 0; i<nIgnCells_; i++)
					T[ignCellsList_[i]] = ignitionTemperature_;
	}
	
private:

	Switch active_;
	Foam::vector position_;
	scalar ignitionTime_;
	scalar duration_;
	scalar ignitionTemperature_;
	scalar ignitionDiameter_;
	labelList ignCellsList_;
    	label nIgnCells_;
};
