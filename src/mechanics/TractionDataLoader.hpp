

#ifndef TRACTIONDATALOADER_HPP_
#define TRACTIONDATALOADER_HPP_

#include <vector>
#include <map>
#include "UblasCustomFunctions.hpp"
#include "PottsArbitrarySurfaceIn3DMesh.hpp"
#include "CellLabel.hpp"






/**
 * This class reads a HemeLB extracted properties file containing information
 * about traction and tangetial component of traction on the surface of a vessel
 * and stores it as data on the lattice sites of a PottsArbitrarySurfaceIn3DMesh.
 */

class TractionDataLoader
{

public:
	/**
	 * Constructor. Takes the name of the property extraction file containing information about tractions
	 *
	 * @param tractionFilename file name
	 */
	TractionDataLoader(const std::string& tractionFilename);

	/**
	 * Maps traction information onto the lattice sites of a Potts mesh
	 *
	 * @param pPottsMesh potts mesh
	 */
	void UpdateLatticeSiteData(PottsArbitrarySurfaceIn3DMesh<3>* pPottsMesh);

private:
	/**
	 * 	Store the loaded traction data; this is the applied traction at a given point.
	 */
	std::vector < c_vector<double,3> > mAppliedTractions;

	/**
	 * 	Store the loaded tangential component of the traction data; this is the applied tangential traction at a given point.
	 */
	std::vector < c_vector<double,3> > mAppliedTangentTractions;

	/**
	 * 	Store the loaded traction data; this is the location that the traction is defined at.
	 */
	std::vector < c_vector<double,3> > mAppliedPosition;

};

#endif /* TRACTIONDATALOADER_HPP_ */
