/*

Copyright (c) 2005-2014, University of Oxford.
All rights reserved.

University of Oxford means the Chancellor, Masters and Scholars of the
University of Oxford, having an administrative office at Wellington
Square, Oxford OX1 2JD, UK.

This file is part of Chaste.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:
 * Redistributions of source code must retain the above copyright notice,
   this list of conditions and the following disclaimer.
 * Redistributions in binary form must reproduce the above copyright notice,
   this list of conditions and the following disclaimer in the documentation
   and/or other materials provided with the distribution.
 * Neither the name of the University of Oxford nor the names of its
   contributors may be used to endorse or promote products derived from this
   software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE
GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT
OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

*/

#ifndef ARBITRARYRACDEPENDENTADHESIONUPDATERULE_HPP_
#define ARBITRARYRACDEPENDENTADHESIONUPDATERULE_HPP_

#include "ChasteSerialization.hpp"
#include <boost/serialization/base_object.hpp>

#include "AbstractPottsUpdateRule.hpp"
#include "PottsBasedCellPopulation.hpp"

// Needed here to avoid serialization errors (on Boost<1.37)
#include "CellLabel.hpp"

/**
 * An adhesion update rule for use in cell-based simulations
 * using the cellular Potts model.
 */
template<unsigned DIM>
class ArbitraryRacDependentAdhesionPottsUpdateRule : public AbstractPottsUpdateRule<DIM>
{
friend class TestPottsUpdateRules;

//EXPLICAR APROPIADAMENTE
private:
/**
     * Node-node adhesion energy parameter.
     * Set to the default value 0.1 in the constructor.
     * \todo provide units
     */
    double mNodeNode_RacHighHigh_AdhesionEnergyParameter;

    /**
     * Labelled node-node adhesion energy parameter.
     * Set to the default value 0.2 in the constructor.
     * \todo provide units
     */
    double mLabelledNodeNode_OR_RacHighLow_AdhesionEnergyParameter;

    /**
     * Labelled node-labelled node adhesion energy parameter.
     * Set to the default value 0.4 in the constructor.
     * \todo provide units
     */
    double mLabelledNodeLabelledNode_OR_RacLowLow_AdhesionEnergyParameter;

    /**
     * Labelled node-labelled node adhesion energy parameter.
     * Set to the default value 0.4 in the constructor.
     * \todo provide units
     */
    double mOne_LabelledNode_OR_RacLow_AdhesionEnergyParameter;

    /**
     * Labelled node-labelled node adhesion energy parameter.
     * Set to the default value 0.4 in the constructor.
     * \todo provide units
     */
    double mLabelledNodeLabelledNode_AND_RacLowLow_AdhesionEnergyParameter;


    friend class boost::serialization::access;
    /**
     * Boost Serialization method for archiving/checkpointing.
     * Archives the object and its member variables.
     *
     * @param archive  The boost archive.
     * @param version  The current version of this class.
     */
    template<class Archive>
    void serialize(Archive & archive, const unsigned int version)
    {
        archive & boost::serialization::base_object<AbstractPottsUpdateRule<DIM> >(*this);
		archive & mNodeNode_RacHighHigh_AdhesionEnergyParameter;
		archive & mLabelledNodeNode_OR_RacHighLow_AdhesionEnergyParameter;
		archive & mLabelledNodeLabelledNode_OR_RacLowLow_AdhesionEnergyParameter;
		archive & mOne_LabelledNode_OR_RacLow_AdhesionEnergyParameter;
		archive & mLabelledNodeLabelledNode_AND_RacLowLow_AdhesionEnergyParameter;
    }

public:

    /**
     * Constructor.
     */
    ArbitraryRacDependentAdhesionPottsUpdateRule();

    /**
     * Destructor.
     */
    virtual ~ArbitraryRacDependentAdhesionPottsUpdateRule();


    double GetRacValueForNode(unsigned nodeIndex, PottsBasedCellPopulation<DIM>& rCellPopulation);


    /**
     * Overridden EvaluateHamiltonianContribution() method
     *
     * Uses  sum_adjacentsites (1-delta(spin(i),spin(j))) gamma(spin(i),spin(j))
     *
     * @param currentNodeIndex The index of the current node/lattice site
     * @param targetNodeIndex The index of the target node/lattice site
     * @param rCellPopulation The cell population
     *
     * @return The difference in the Hamiltonian with the configuration of the target node
     * having the same spin as the current node with the current configuration. i.e H_1-H_0
     */
    double EvaluateHamiltonianContribution(unsigned currentNodeIndex,
                                           unsigned targetNodeIndex,
                                           PottsBasedCellPopulation<DIM>& rCellPopulation);

    /**
           * Method to calculate the adhesion type of a node.
           *
           * @param pNodeA pointer to the 1st node
           * @param pNodeB pointer to the 2nd node
           *
           * @return The adhesion type of the node
           */

    virtual double GetNodeNodeAdhesionEnergy(double nodeAdhesionTypeA, unsigned nodeIndex, PottsBasedCellPopulation<DIM>& rCellPopulation);

     /**
         * Method to calculate the adhesion type of a node.
         *
         * @param pNodeA pointer to the 1st node
         * @param pNodeB pointer to the 2nd node
         *
         * @return The adhesion type of the node
         */

      virtual double GetNodeAdhesionType(unsigned nodeIndex, PottsBasedCellPopulation<DIM>& rCellPopulation );

       /**
        * Method to calculate the specific interaction between cell and medium can be overridden in
        * child classes to  implement differential adhesion .etc.
        *
        * @param pCell pointer to the cell
        *
        * @return Cell boundary interaction adhesion energy for the cell
        */
       virtual double GetCellBoundaryAdhesionEnergy(CellPtr pCell);




       /** GET
         * @return mCellBoundaryAdhesionEnergyParameter
         */
       double GetNodeNode_RacHighHigh_AdhesionEnergyParameter();

       /**
        * @return mLabelledNodeNode_OR_RacHighLow_AdhesionEnergyParameter
        */
       double GetLabelledNodeNode_OR_RacHighLow_AdhesionEnergyParameter();

       /**
             * @return mLabelledNodeLabelledNode_OR_RacLowLow_AdhesionEnergyParameter
             */
        double GetLabelledNodeLabelledNode_OR_RacLowLow_AdhesionEnergyParameter();

       /**
             * @return mOne_LabelledNode_OR_RacLow_AdhesionEnergyParameter
             */
        double GetOne_LabelledNode_OR_RacLow_AdhesionEnergyParameter();

        /**
                   * @return mLabelledNodeLabelledNode_AND_RacLowLow_AdhesionEnergyParameter
                   */
       double GetLabelledNodeLabelledNode_AND_RacLowLow_AdhesionEnergyParameter();




       /** SET
               * Set mNodeNode_RacHighHigh_AdhesionEnergyParameter
               *
               */

       double SetNodeNode_RacHighHigh_AdhesionEnergyParameter();

       /**
               * Set mLabelledNodeNode_OR_RacHighLow_AdhesionEnergyParameter
        */
       double SetLabelledNodeNode_OR_RacHighLow_AdhesionEnergyParameter();

       /**
               * Set mLabelledNodeLabelledNode_OR_RacLowLow_AdhesionEnergyParameter
             */
        double SetLabelledNodeLabelledNode_OR_RacLowLow_AdhesionEnergyParameter();

       /**
               * Set mOne_LabelledNode_OR_RacLow_AdhesionEnergyParameter
             */
        double SetOne_LabelledNode_OR_RacLow_AdhesionEnergyParameter();

        /**
               * Set mLabelledNodeLabelledNode_AND_RacLowLow_AdhesionEnergyParameter
                   */
       double SetLabelledNodeLabelledNode_AND_RacLowLow_AdhesionEnergyParameter();




    /**
     * Overridden OutputUpdateRuleParameters() method.
     *
     * @param rParamsFile the file stream to which the parameters are output
     */
    void OutputUpdateRuleParameters(out_stream& rParamsFile);
};

#include "SerializationExportWrapper.hpp"
EXPORT_TEMPLATE_CLASS_SAME_DIMS(ArbitraryRacDependentAdhesionPottsUpdateRule)

#endif /*ARBITRARYRACDEPENDENTADHESIONUPDATERULE_HPP_*/
