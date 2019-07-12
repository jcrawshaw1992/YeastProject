#ifndef FLOWBASEDCDC42STIMULUSPROTOCOL_HPP
#define FLOWBASEDCDC42STIMULUSPROTOCOL_HPP

#include "AbstractCDC42StimulusProtocol.hpp"
#include <map>

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
class FlowBasedCDC42StimulusProtocol : public AbstractCDC42StimulusProtocol<ELEMENT_DIM, SPACE_DIM>
{
public:
    FlowBasedCDC42StimulusProtocol(const TetrahedralMesh<ELEMENT_DIM, SPACE_DIM>& rMesh, double stimulusDuration, const c_vector<double,ELEMENT_DIM>& flowDirection);

    virtual ~FlowBasedCDC42StimulusProtocol();

    double GetCDC42ActivationRate(const Node<SPACE_DIM>* pNode, double t) const;

private:
    typedef std::map<const Node<SPACE_DIM>*, c_vector<double,ELEMENT_DIM> > EdgeMapType;

    void PopulateEdgeNormalForNodesOnBoundary(const TetrahedralMesh<ELEMENT_DIM,SPACE_DIM>& rMesh);

    c_vector<double,ELEMENT_DIM> mFlowDirection;

    EdgeMapType mEdgeNormalForNodesOnBoundary;
};

#endif
