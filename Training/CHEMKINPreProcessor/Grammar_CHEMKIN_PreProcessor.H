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
|-------------------------------------------------------------------------|
|                                                                         |
|   This file is part of OpenSMOKE++ framework.                           |
|                                                                         |
|	License                                                           |
|                                                                         |
|   Copyright(C) 2017-2012  Alberto Cuoci                                 |
|   OpenSMOKE++ is free software: you can redistribute it and/or modify   |
|   it under the terms of the GNU General Public License as published by  |
|   the Free Software Foundation, either version 3 of the License, or     |
|   (at your option) any later version.                                   |
|                                                                         |
|   OpenSMOKE++ is distributed in the hope that it will be useful,        |
|   but WITHOUT ANY WARRANTY; without even the implied warranty of        |
|   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the         |
|   GNU General Public License for more details.                          |
|                                                                         |
|   You should have received a copy of the GNU General Public License     |
|   along with OpenSMOKE++. If not, see <http://www.gnu.org/licenses/>.   |
|                                                                         |
\*-----------------------------------------------------------------------*/

#include "dictionary/OpenSMOKE_DictionaryManager.h"
#include "dictionary/OpenSMOKE_DictionaryGrammar.h"
#include "dictionary/OpenSMOKE_DictionaryKeyWord.h"

class GrammarKineticInterpreter : public OpenSMOKE::OpenSMOKE_DictionaryGrammar
{
protected:

	virtual void DefineRules()
	{
		AddKeyWord( OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@Thermodynamics", 
															OpenSMOKE::SINGLE_PATH, 
															"Name of the file (ASCII) containing the thermodynamic data (CHEMKIN format)", 
															true) );

		AddKeyWord( OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@Output", 
															OpenSMOKE::SINGLE_PATH, 
															"Name of the folder where the pre-processed data are written", 
															true) );
		
		AddKeyWord( OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@Transport", 
															OpenSMOKE::SINGLE_PATH, 
															"Name of the file (ASCII) containing the transport data (CHEMKIN format)", 
															false) );

		AddKeyWord( OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@Kinetics", 
															OpenSMOKE::SINGLE_PATH, 
															"Name of the file (ASCII) containing the kinetic mechanism (CHEMKIN format)", 
															false) );

		AddKeyWord( OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@CheckThermodynamics", 
															OpenSMOKE::SINGLE_BOOL, 
															"The thermodynamic data are checked and additional files are written in output (together with a consistent reformulation of thermodynamic data)", 
															false) );

		AddKeyWord(OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@SpeciesBundling",
															OpenSMOKE::SINGLE_BOOL,
															"Species bundling is applied according to the following sequence of relative errors: { 0.01, 0.025, 0.050, 0.075, 0.1, 0.25, 0.50 }",
															false,
															"none",
															"@Transport",
															"none"));
		
		AddKeyWord( OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@OutputOldStyle", 
															OpenSMOKE::SINGLE_BOOL, 
															"The output file are written also using the old format used by the OpenSMOKE framework (before 2013)", 
															false,
															"none",
															"@Transport @Kinetics",
															"none") );

		AddKeyWord( OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@TransportFittingCoefficients", 
															OpenSMOKE::SINGLE_BOOL, 
															"The fitting coefficients for the transport properties are written on a file. Please consider that this operation is very slow for large kinetic mechanisms (more than 1000 species) and produces huge files!", 
															false,
															"none",
															"@Transport",
															"none") );

		AddKeyWord( OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@ReactionTables", 
															OpenSMOKE::SINGLE_BOOL, 
															"For each reaction detailed information is reported on a file (kinetic constants, change of moles, etc.)", 
															false,
															"none",
															"@Kinetics",
															"none") );

		AddKeyWord(OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@ReactionTablesListOfTemperatures",
															OpenSMOKE::VECTOR_STRING,
															"List of temperatures at which the detailed information is reported on a file for each reaction (kinetic constants, change of moles, etc.)",
															false,
															"none",
															"@Kinetics @ReactionTables",
															"none"));

		AddKeyWord( OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@ReverseFitting", 
															OpenSMOKE::SINGLE_BOOL, 
															"For each reversible reaction the reverse kinetic constants are estimated assuming the Arrhenius' law", 
															false,
															"none",
															"@Kinetics",
															"none") );

		AddKeyWord(OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@SparsityPatternAnalysis",
															OpenSMOKE::SINGLE_BOOL,
															"The analysis of sparsity pattern of kinetic mechanism is automatically performed after the pre-processing",
															false,
															"none",
															"@Kinetics",
															"none"));

		AddKeyWord( OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@Comments", 
															OpenSMOKE::SINGLE_DICTIONARY, 
															"Additional data (author name, comments, etc.) can be added to the pre-processed kinetic mechanism", 
															false) );

		AddKeyWord(OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@RewriteCHEMKIN",
															OpenSMOKE::SINGLE_BOOL,
															"Rewrites the pre-processed kinetic mechanism in CHEMKIN format after correction on atomic balances (default: true)",
															false) );
		
	}
};

