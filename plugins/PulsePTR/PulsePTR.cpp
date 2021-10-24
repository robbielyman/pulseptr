// PluginPulsePTR.cpp
// Rylee Alanza Lyman (ryleealanza@gmail.com)

#include "SC_PlugIn.hpp"
#include "PulsePTR.hpp"

static InterfaceTable* ft;

namespace PulsePTR {

PulsePTR::PulsePTR() {
    mCalcFunc = make_calc_function<PulsePTR, &PulsePTR::next>();
    next(1);
}

void PulsePTR::next(int nSamples) {
    const float* input = in(0);
    const float* gain = in(1);
    float* outbuf = out(0);

    // simple gain function
    for (int i = 0; i < nSamples; ++i) {
        outbuf[i] = input[i] * gain[i];
    }
}

} // namespace PulsePTR

PluginLoad(PulsePTRUGens) {
    // Plugin magic
    ft = inTable;
    registerUnit<PulsePTR::PulsePTR>(ft, "PulsePTR", false);
}
