/*
Jess's abstract wrapped update rule

*/

#include "AbstractWrappedPottsUpdateRule.hpp"

template<unsigned DIM>
AbstractWrappedPottsUpdateRule<DIM>::AbstractWrappedPottsUpdateRule()
    : AbstractUpdateRule<DIM>()
{
}

template<unsigned DIM>
AbstractWrappedPottsUpdateRule<DIM>::~AbstractWrappedPottsUpdateRule()
{
}

template<unsigned DIM>
void AbstractWrappedPottsUpdateRule<DIM>::OutputUpdateRuleParameters(out_stream& rParamsFile)
{
    // Call method on direct parent class
    AbstractUpdateRule<DIM>::OutputUpdateRuleParameters(rParamsFile);
}

// Explicit instantiation
template class AbstractWrappedPottsUpdateRule<1>;
template class AbstractWrappedPottsUpdateRule<2>;
template class AbstractWrappedPottsUpdateRule<3>;
