// PluginPulsePTR.hpp
// Rylee Alanza Lyman (ryleealanza@gmail.com)

#pragma once

#include "SC_PlugIn.hpp"

namespace PulsePTR {

class PulsePTR : public SCUnit {
public:
    PulsePTR();

    // Destructor
    // ~PulsePTR();

private:
    // Calc function
    void next(int nSamples);

    // Member variables
};

} // namespace PulsePTR
