// PluginPulsePTR.cpp
// Rylee Alanza Lyman (ryleealanza@gmail.com)

#include "SC_PlugIn.hpp"
#include "PulsePTR.hpp"

static InterfaceTable* ft;

namespace PulsePTR {

PulsePTR::PulsePTR() {
    if (isAudioRateIn(0)) {
        if (isAudioRateIn(1)) {
            if (isAudioRateIn(2)) {
                if (isAudioRateIn(3))   mCalcFunc = make_calc_function<PulsePTR, &PulsePTR::next_aaaa>();
                else                    mCalcFunc = make_calc_function<PulsePTR, &PulsePTR::next_aaak>();
            }
            else {
                if (isAudioRateIn(3))   mCalcFunc = make_calc_function<PulsePTR, &PulsePTR::next_aaka>();
                else                    mCalcFunc = make_calc_function<PulsePTR, &PulsePTR::next_aakk>();
            }
        } 
        else {
            if (isAudioRateIn(2)) {
                if (isAudioRateIn(3))   mCalcFunc = make_calc_function<PulsePTR, &PulsePTR::next_akaa>();
                else                    mCalcFunc = make_calc_function<PulsePTR, &PulsePTR::next_akak>();
            }
            else {
                if (isAudioRateIn(3))   mCalcFunc = make_calc_function<PulsePTR, &PulsePTR::next_akka>();
                else                    mCalcFunc = make_calc_function<PulsePTR, &PulsePTR::next_akkk>();
            }
        }
    }
    else {
        if (isAudioRateIn(1)) {
            if (isAudioRateIn(2)) {
                if (isAudioRateIn(3))   mCalcFunc = make_calc_function<PulsePTR, &PulsePTR::next_kaaa>();
                else                    mCalcFunc = make_calc_function<PulsePTR, &PulsePTR::next_kaak>();
            }
            else {
                if (isAudioRateIn(3))   mCalcFunc = make_calc_function<PulsePTR, &PulsePTR::next_kaka>();
                else                    mCalcFunc = make_calc_function<PulsePTR, &PulsePTR::next_kakk>();
            }
        } 
        else {
            if (isAudioRateIn(2)) {
                if (isAudioRateIn(3))   mCalcFunc = make_calc_function<PulsePTR, &PulsePTR::next_kkaa>();
                else                    mCalcFunc = make_calc_function<PulsePTR, &PulsePTR::next_kkak>();
            }
            else {
                if (isAudioRateIn(3))   mCalcFunc = make_calc_function<PulsePTR, &PulsePTR::next_kkka>();
                else                    mCalcFunc = make_calc_function<PulsePTR, &PulsePTR::next_kkkk>();
            }
        }
    }
    next_kkkk(1);
}

// pulse wave algorithm: Peter McCulloch
inline double PulsePTR::algorithm(double p, double t0, double t2, double t3,
        double w, double p1, double p2, double p3) {
    if (p < w) {
        if (p < t0) {
            // y = -P3*x^3 + w
            return -p3 * p * p * p + w;
        }
        else if (p < t2) {
            // y = 2*P3*x^3 - 3*P2*x^2 + 1.5*P1*x + w - 0.5
            return 2.f * p3 * p * p * p - 3.f * p2 * p * p + 1.5f * p1 * p + w - 0.5f;
        }
        else if (p < t3) {
            // y = -P3*x^3 + 3*P2*x^2 - 4.5*P1*x + w + 3.5
            return -p3 * p * p * p + 3.f * p2 * p * p - 4.5f * p1 * p + w + 3.5f;
        }
        else {
            return w - 1.f;
        }
    }
    else {
        double pw = p - w;
        if (pw < t0) {
            // y = P3*x^3 + w - 1
            return p3 * pw * pw * pw + w - 1;
        }
        else if (pw < t2) {
            // y = -2*P3*x^3 + 3*P2*x^2 - 1.5*P1*x + w - 0.5
            return -2.f * p3 * pw * pw * pw + 3.f * p2 * pw * pw - 1.5f * p1 * pw + w - 0.5f;
        }
        else if (pw < t3) {
            // y = P3*x^3 - 3*P2*x^2 + 4.5*P1*x + w - 4.5
            return p3 * pw * pw * pw - 3.f * p2 * pw * pw + 4.5f * p1 * pw + w - 4.5f;
        }
        else {
            return w;
        }
    }
}

void PulsePTR::next_aaaa(int nSamples) {
    const float* freq   = in(0);
    const float* phase  = in(1);
    const float* sync   = in(2);
    const float* width  = in(3);
    float* outbuf   = out(0);
    float* syncOut  = out(1);

    double pos      = mPhase;
    double rateinv  = 1 / sampleRate();
    double phasein  = mPhaseIn;
    double lastsync = mSync;
    for (int i = 0; i < nSamples; ++i) {
        double step     = freq[i] * rateinv;
        double step2    = 2.f*step;
        double step3    = 3.f*step;
        double samples  = 1.f / step;
        double w        = width[i] < 0.01f ? 0.01f : width[i] > 0.99f ? 0.99f : width[i];
        double p2       = samples * samples / 2.f;
        double p3       = samples * samples * samples / 6.f;
        double out1     = 0.f;
        double out2     = 0.f;
        if (sync[i] > 0.f && lastsync == 0.f) {
            pos = 0.f;
        }
        else {
            pos = pos + (phase[i] - phasein) + step;
            while (pos >= 1.f) {
                pos -= 1.f;
            }
            while (pos < 0.f) {
                pos += 1.f;
            }
        }
        lastsync    = sync[i];
        phasein     = phase[i];
        if (pos < step) out2 = 1.f;
        out1 = 2.f * algorithm(pos, step, step2, step3, w, samples, p2, p3);

        syncOut[i] = out2;
        outbuf[i]  = out1;
    }
    mPhase      = pos;
    mSync       = lastsync;
    mPhaseIn    = phasein;
}

void PulsePTR::next_aaak(int nSamples) {
    const float* freq   = in(0);
    const float* phase  = in(1);
    const float* sync   = in(2);
    const float width   = in0(3);
    float* outbuf   = out(0);
    float* syncOut  = out(1);

    double pos      = mPhase;
    double rateinv  = 1 / sampleRate();
    double phasein  = mPhaseIn;
    double lastsync = mSync;
    double w        = width < 0.01f ? 0.01f : width > 0.99f ? 0.99f : width;
    for (int i = 0; i < nSamples; ++i) {
        double step     = freq[i] * rateinv;
        double step2    = 2.f*step;
        double step3    = 3.f*step;
        double samples  = 1.f / step;
        double p2       = samples * samples / 2.f;
        double p3       = samples * samples * samples / 6.f;
        double out1     = 0.f;
        double out2     = 0.f;
        if (sync[i] > 0.f && lastsync == 0.f) {
            pos = 0.f;
        }
        else {
            pos = pos + (phase[i] - phasein) + step;
            while (pos >= 1.f) {
                pos -= 1.f;
            }
            while (pos < 0.f) {
                pos += 1.f;
            }
        }
        lastsync    = sync[i];
        phasein     = phase[i];
        if (pos < step) out2 = 1.f;
        out1 = 2.f * algorithm(pos, step, step2, step3, w, samples, p2, p3);

        syncOut[i] = out2;
        outbuf[i]  = out1;
    }
    mPhase      = pos;
    mSync       = lastsync;
    mPhaseIn    = phasein;
}

void PulsePTR::next_aaka(int nSamples) {
    const float* freq   = in(0);
    const float* phase  = in(1);
    const float sync    = in0(2);
    const float* width  = in(3);
    float* outbuf   = out(0);
    float* syncOut  = out(1);

    double pos      = mPhase;
    double rateinv  = 1 / sampleRate();
    double phasein  = mPhaseIn;
    double lastsync = mSync;
    for (int i = 0; i < nSamples; ++i) {
        double step     = freq[i] * rateinv;
        double step2    = 2.f*step;
        double step3    = 3.f*step;
        double samples  = 1.f / step;
        double w        = width[i] < 0.01f ? 0.01f : width[i] > 0.99f ? 0.99f : width[i];
        double p2       = samples * samples / 2.f;
        double p3       = samples * samples * samples / 6.f;
        double out1     = 0.f;
        double out2     = 0.f;
        if (sync > 0.f && lastsync == 0.f) {
            pos = 0.f;
        }
        else {
            pos = pos + (phase[i] - phasein) + step;
            while (pos >= 1.f) {
                pos -= 1.f;
            }
            while (pos < 0.f) {
                pos += 1.f;
            }
        }
        lastsync    = sync;
        phasein     = phase[i];
        if (pos < step) out2 = 1.f;
        out1 = 2.f * algorithm(pos, step, step2, step3, w, samples, p2, p3);

        syncOut[i] = out2;
        outbuf[i]  = out1;
    }
    mPhase      = pos;
    mSync       = lastsync;
    mPhaseIn    = phasein;
}

void PulsePTR::next_aakk(int nSamples) {
    const float* freq   = in(0);
    const float* phase  = in(1);
    const float sync    = in0(2);
    const float width   = in0(3);
    float* outbuf   = out(0);
    float* syncOut  = out(1);

    double pos      = mPhase;
    double rateinv  = 1 / sampleRate();
    double phasein  = mPhaseIn;
    double lastsync = mSync;
    double w        = width < 0.01f ? 0.01f : width > 0.99f ? 0.99f : width;
    for (int i = 0; i < nSamples; ++i) {
        double step     = freq[i] * rateinv;
        double step2    = 2.f*step;
        double step3    = 3.f*step;
        double samples  = 1.f / step;
        double p2       = samples * samples / 2.f;
        double p3       = samples * samples * samples / 6.f;
        double out1     = 0.f;
        double out2     = 0.f;
        if (sync > 0.f && lastsync == 0.f) {
            pos = 0.f;
        }
        else {
            pos = pos + (phase[i] - phasein) + step;
            while (pos >= 1.f) {
                pos -= 1.f;
            }
            while (pos < 0.f) {
                pos += 1.f;
            }
        }
        lastsync    = sync;
        phasein     = phase[i];
        if (pos < step) out2 = 1.f;
        out1 = 2.f * algorithm(pos, step, step2, step3, w, samples, p2, p3);

        syncOut[i] = out2;
        outbuf[i]  = out1;
    }
    mPhase      = pos;
    mSync       = lastsync;
    mPhaseIn    = phasein;
}

void PulsePTR::next_akaa(int nSamples) {
    const float* freq   = in(0);
    const float phase   = in0(1);
    const float* sync   = in(2);
    const float* width  = in(3);
    float* outbuf   = out(0);
    float* syncOut  = out(1);

    double pos      = mPhase;
    double rateinv  = 1 / sampleRate();
    double phasein  = mPhaseIn;
    double lastsync = mSync;
    for (int i = 0; i < nSamples; ++i) {
        double step     = freq[i] * rateinv;
        double step2    = 2.f*step;
        double step3    = 3.f*step;
        double samples  = 1.f / step;
        double w        = width[i] < 0.01f ? 0.01f : width[i] > 0.99f ? 0.99f : width[i];
        double p2       = samples * samples / 2.f;
        double p3       = samples * samples * samples / 6.f;
        double out1     = 0.f;
        double out2     = 0.f;
        if (sync[i] > 0.f && lastsync == 0.f) {
            pos = 0.f;
        }
        else {
            pos = pos + (phase - phasein) + step;
            while (pos >= 1.f) {
                pos -= 1.f;
            }
            while (pos < 0.f) {
                pos += 1.f;
            }
        }
        lastsync    = sync[i];
        phasein     = phase;
        if (pos < step) out2 = 1.f;
        out1 = 2.f * algorithm(pos, step, step2, step3, w, samples, p2, p3);

        syncOut[i] = out2;
        outbuf[i]  = out1;
    }
    mPhase      = pos;
    mSync       = lastsync;
    mPhaseIn    = phasein;
}

void PulsePTR::next_akak(int nSamples) {
    const float* freq   = in(0);
    const float phase   = in0(1);
    const float* sync   = in(2);
    const float width   = in0(3);
    float* outbuf   = out(0);
    float* syncOut  = out(1);

    double pos      = mPhase;
    double rateinv  = 1 / sampleRate();
    double phasein  = mPhaseIn;
    double lastsync = mSync;
    double w        = width < 0.01f ? 0.01f : width > 0.99f ? 0.99f : width;
    for (int i = 0; i < nSamples; ++i) {
        double step     = freq[i] * rateinv;
        double step2    = 2.f*step;
        double step3    = 3.f*step;
        double samples  = 1.f / step;
        double p2       = samples * samples / 2.f;
        double p3       = samples * samples * samples / 6.f;
        double out1     = 0.f;
        double out2     = 0.f;
        if (sync[i] > 0.f && lastsync == 0.f) {
            pos = 0.f;
        }
        else {
            pos = pos + (phase - phasein) + step;
            while (pos >= 1.f) {
                pos -= 1.f;
            }
            while (pos < 0.f) {
                pos += 1.f;
            }
        }
        lastsync    = sync[i];
        phasein     = phase;
        if (pos < step) out2 = 1.f;
        out1 = 2.f * algorithm(pos, step, step2, step3, w, samples, p2, p3);

        syncOut[i] = out2;
        outbuf[i]  = out1;
    }
    mPhase      = pos;
    mSync       = lastsync;
    mPhaseIn    = phasein;
}

void PulsePTR::next_akka(int nSamples) {
    const float* freq   = in(0);
    const float phase   = in0(1);
    const float sync    = in0(2);
    const float* width  = in(3);
    float* outbuf   = out(0);
    float* syncOut  = out(1);

    double pos      = mPhase;
    double rateinv  = 1 / sampleRate();
    double phasein  = mPhaseIn;
    double lastsync = mSync;
    for (int i = 0; i < nSamples; ++i) {
        double step     = freq[i] * rateinv;
        double step2    = 2.f*step;
        double step3    = 3.f*step;
        double samples  = 1.f / step;
        double w        = width[i] < 0.01f ? 0.01f : width[i] > 0.99f ? 0.99f : width[i];
        double p2       = samples * samples / 2.f;
        double p3       = samples * samples * samples / 6.f;
        double out1     = 0.f;
        double out2     = 0.f;
        if (sync > 0.f && lastsync == 0.f) {
            pos = 0.f;
        }
        else {
            pos = pos + (phase - phasein) + step;
            while (pos >= 1.f) {
                pos -= 1.f;
            }
            while (pos < 0.f) {
                pos += 1.f;
            }
        }
        lastsync    = sync;
        phasein     = phase;
        if (pos < step) out2 = 1.f;
        out1 = 2.f * algorithm(pos, step, step2, step3, w, samples, p2, p3);

        syncOut[i] = out2;
        outbuf[i]  = out1;
    }
    mPhase      = pos;
    mSync       = lastsync;
    mPhaseIn    = phasein;
}

void PulsePTR::next_akkk(int nSamples) {
    const float* freq   = in(0);
    const float phase   = in0(1);
    const float sync    = in0(2);
    const float width   = in0(3);
    float* outbuf   = out(0);
    float* syncOut  = out(1);

    double pos      = mPhase;
    double rateinv  = 1 / sampleRate();
    double phasein  = mPhaseIn;
    double lastsync = mSync;
    double w        = width < 0.01f ? 0.01f : width > 0.99f ? 0.99f : width;
    for (int i = 0; i < nSamples; ++i) {
        double step     = freq[i] * rateinv;
        double step2    = 2.f*step;
        double step3    = 3.f*step;
        double samples  = 1.f / step;
        double p2       = samples * samples / 2.f;
        double p3       = samples * samples * samples / 6.f;
        double out1     = 0.f;
        double out2     = 0.f;
        if (sync > 0.f && lastsync == 0.f) {
            pos = 0.f;
        }
        else {
            pos = pos + (phase - phasein) + step;
            while (pos >= 1.f) {
                pos -= 1.f;
            }
            while (pos < 0.f) {
                pos += 1.f;
            }
        }
        lastsync    = sync;
        phasein     = phase;
        if (pos < step) out2 = 1.f;
        out1 = 2.f * algorithm(pos, step, step2, step3, w, samples, p2, p3);

        syncOut[i] = out2;
        outbuf[i]  = out1;
    }
    mPhase      = pos;
    mSync       = lastsync;
    mPhaseIn    = phasein;
}

void PulsePTR::next_kaaa(int nSamples) {
    const float freq    = in0(0);
    const float* phase  = in(1);
    const float* sync   = in(2);
    const float* width  = in(3);
    float* outbuf   = out(0);
    float* syncOut  = out(1);

    double pos      = mPhase;
    double rateinv  = 1 / sampleRate();
    double phasein  = mPhaseIn;
    double lastsync = mSync;
    double step     = freq * rateinv;
    double step2    = 2.f*step;
    double step3    = 3.f*step;
    double samples  = 1.f / step;
    double p2       = samples * samples / 2.f;
    double p3       = samples * samples * samples / 6.f;
    for (int i = 0; i < nSamples; ++i) {
        double w        = width[i] < 0.01f ? 0.01f : width[i] > 0.99f ? 0.99f : width[i];
        double out1     = 0.f;
        double out2     = 0.f;
        if (sync[i] > 0.f && lastsync == 0.f) {
            pos = 0.f;
        }
        else {
            pos = pos + (phase[i] - phasein) + step;
            while (pos >= 1.f) {
                pos -= 1.f;
            }
            while (pos < 0.f) {
                pos += 1.f;
            }
        }
        lastsync    = sync[i];
        phasein     = phase[i];
        if (pos < step) out2 = 1.f;
        out1 = 2.f * algorithm(pos, step, step2, step3, w, samples, p2, p3);

        syncOut[i] = out2;
        outbuf[i]  = out1;
    }
    mPhase      = pos;
    mSync       = lastsync;
    mPhaseIn    = phasein;
}

void PulsePTR::next_kaak(int nSamples) {
    const float freq    = in0(0);
    const float* phase  = in(1);
    const float* sync   = in(2);
    const float width   = in0(3);
    float* outbuf   = out(0);
    float* syncOut  = out(1);

    double pos      = mPhase;
    double rateinv  = 1 / sampleRate();
    double phasein  = mPhaseIn;
    double lastsync = mSync;
    double w        = width < 0.01f ? 0.01f : width > 0.99f ? 0.99f : width;
    double step     = freq * rateinv;
    double step2    = 2.f*step;
    double step3    = 3.f*step;
    double samples  = 1.f / step;
    double p2       = samples * samples / 2.f;
    double p3       = samples * samples * samples / 6.f;
    for (int i = 0; i < nSamples; ++i) {
        double out1     = 0.f;
        double out2     = 0.f;
        if (sync[i] > 0.f && lastsync == 0.f) {
            pos = 0.f;
        }
        else {
            pos = pos + (phase[i] - phasein) + step;
            while (pos >= 1.f) {
                pos -= 1.f;
            }
            while (pos < 0.f) {
                pos += 1.f;
            }
        }
        lastsync    = sync[i];
        phasein     = phase[i];
        if (pos < step) out2 = 1.f;
        out1 = 2.f * algorithm(pos, step, step2, step3, w, samples, p2, p3);

        syncOut[i] = out2;
        outbuf[i]  = out1;
    }
    mPhase      = pos;
    mSync       = lastsync;
    mPhaseIn    = phasein;
}

void PulsePTR::next_kaka(int nSamples) {
    const float freq    = in0(0);
    const float* phase  = in(1);
    const float sync    = in0(2);
    const float* width  = in(3);
    float* outbuf   = out(0);
    float* syncOut  = out(1);

    double pos      = mPhase;
    double rateinv  = 1 / sampleRate();
    double phasein  = mPhaseIn;
    double lastsync = mSync;
    double step     = freq * rateinv;
    double step2    = 2.f*step;
    double step3    = 3.f*step;
    double samples  = 1.f / step;
    double p2       = samples * samples / 2.f;
    double p3       = samples * samples * samples / 6.f;
    for (int i = 0; i < nSamples; ++i) {
        double w        = width[i] < 0.01f ? 0.01f : width[i] > 0.99f ? 0.99f : width[i];
        double out1     = 0.f;
        double out2     = 0.f;
        if (sync > 0.f && lastsync == 0.f) {
            pos = 0.f;
        }
        else {
            pos = pos + (phase[i] - phasein) + step;
            while (pos >= 1.f) {
                pos -= 1.f;
            }
            while (pos < 0.f) {
                pos += 1.f;
            }
        }
        lastsync    = sync;
        phasein     = phase[i];
        if (pos < step) out2 = 1.f;
        out1 = 2.f * algorithm(pos, step, step2, step3, w, samples, p2, p3);

        syncOut[i] = out2;
        outbuf[i]  = out1;
    }
    mPhase      = pos;
    mSync       = lastsync;
    mPhaseIn    = phasein;
}

void PulsePTR::next_kakk(int nSamples) {
    const float freq    = in0(0);
    const float* phase  = in(1);
    const float sync    = in0(2);
    const float width   = in0(3);
    float* outbuf   = out(0);
    float* syncOut  = out(1);

    double pos      = mPhase;
    double rateinv  = 1 / sampleRate();
    double phasein  = mPhaseIn;
    double lastsync = mSync;
    double w        = width < 0.01f ? 0.01f : width > 0.99f ? 0.99f : width;
    double step     = freq * rateinv;
    double step2    = 2.f*step;
    double step3    = 3.f*step;
    double samples  = 1.f / step;
    double p2       = samples * samples / 2.f;
    double p3       = samples * samples * samples / 6.f;
    for (int i = 0; i < nSamples; ++i) {
        double out1     = 0.f;
        double out2     = 0.f;
        if (sync > 0.f && lastsync == 0.f) {
            pos = 0.f;
        }
        else {
            pos = pos + (phase[i] - phasein) + step;
            while (pos >= 1.f) {
                pos -= 1.f;
            }
            while (pos < 0.f) {
                pos += 1.f;
            }
        }
        lastsync    = sync;
        phasein     = phase[i];
        if (pos < step) out2 = 1.f;
        out1 = 2.f * algorithm(pos, step, step2, step3, w, samples, p2, p3);

        syncOut[i] = out2;
        outbuf[i]  = out1;
    }
    mPhase      = pos;
    mSync       = lastsync;
    mPhaseIn    = phasein;
}

void PulsePTR::next_kkaa(int nSamples) {
    const float freq    = in0(0);
    const float phase   = in0(1);
    const float* sync   = in(2);
    const float* width  = in(3);
    float* outbuf   = out(0);
    float* syncOut  = out(1);

    double pos      = mPhase;
    double rateinv  = 1 / sampleRate();
    double phasein  = mPhaseIn;
    double lastsync = mSync;
    double step     = freq * rateinv;
    double step2    = 2.f*step;
    double step3    = 3.f*step;
    double samples  = 1.f / step;
    double p2       = samples * samples / 2.f;
    double p3       = samples * samples * samples / 6.f;
    for (int i = 0; i < nSamples; ++i) {
        double w        = width[i] < 0.01f ? 0.01f : width[i] > 0.99f ? 0.99f : width[i];
        double out1     = 0.f;
        double out2     = 0.f;
        if (sync[i] > 0.f && lastsync == 0.f) {
            pos = 0.f;
        }
        else {
            pos = pos + (phase - phasein) + step;
            while (pos >= 1.f) {
                pos -= 1.f;
            }
            while (pos < 0.f) {
                pos += 1.f;
            }
        }
        lastsync    = sync[i];
        phasein     = phase;
        if (pos < step) out2 = 1.f;
        out1 = 2.f * algorithm(pos, step, step2, step3, w, samples, p2, p3);

        syncOut[i] = out2;
        outbuf[i]  = out1;
    }
    mPhase      = pos;
    mSync       = lastsync;
    mPhaseIn    = phasein;
}

void PulsePTR::next_kkak(int nSamples) {
    const float freq    = in0(0);
    const float phase   = in0(1);
    const float* sync   = in(2);
    const float width   = in0(3);
    float* outbuf   = out(0);
    float* syncOut  = out(1);

    double pos      = mPhase;
    double rateinv  = 1 / sampleRate();
    double phasein  = mPhaseIn;
    double lastsync = mSync;
    double w        = width < 0.01f ? 0.01f : width > 0.99f ? 0.99f : width;
    double step     = freq * rateinv;
    double step2    = 2.f*step;
    double step3    = 3.f*step;
    double samples  = 1.f / step;
    double p2       = samples * samples / 2.f;
    double p3       = samples * samples * samples / 6.f;
    for (int i = 0; i < nSamples; ++i) {
        double out1     = 0.f;
        double out2     = 0.f;
        if (sync[i] > 0.f && lastsync == 0.f) {
            pos = 0.f;
        }
        else {
            pos = pos + (phase - phasein) + step;
            while (pos >= 1.f) {
                pos -= 1.f;
            }
            while (pos < 0.f) {
                pos += 1.f;
            }
        }
        lastsync    = sync[i];
        phasein     = phase;
        if (pos < step) out2 = 1.f;
        out1 = 2.f * algorithm(pos, step, step2, step3, w, samples, p2, p3);

        syncOut[i] = out2;
        outbuf[i]  = out1;
    }
    mPhase      = pos;
    mSync       = lastsync;
    mPhaseIn    = phasein;
}

void PulsePTR::next_kkka(int nSamples) {
    const float freq    = in0(0);
    const float phase   = in0(1);
    const float sync    = in0(2);
    const float* width  = in(3);
    float* outbuf   = out(0);
    float* syncOut  = out(1);

    double pos      = mPhase;
    double rateinv  = 1 / sampleRate();
    double phasein  = mPhaseIn;
    double lastsync = mSync;
    double step     = freq * rateinv;
    double step2    = 2.f*step;
    double step3    = 3.f*step;
    double samples  = 1.f / step;
    double p2       = samples * samples / 2.f;
    double p3       = samples * samples * samples / 6.f;
    for (int i = 0; i < nSamples; ++i) {
        double w        = width[i] < 0.01f ? 0.01f : width[i] > 0.99f ? 0.99f : width[i];
        double out1     = 0.f;
        double out2     = 0.f;
        if (sync > 0.f && lastsync == 0.f) {
            pos = 0.f;
        }
        else {
            pos = pos + (phase - phasein) + step;
            while (pos >= 1.f) {
                pos -= 1.f;
            }
            while (pos < 0.f) {
                pos += 1.f;
            }
        }
        lastsync    = sync;
        phasein     = phase;
        if (pos < step) out2 = 1.f;
        out1 = 2.f * algorithm(pos, step, step2, step3, w, samples, p2, p3);

        syncOut[i] = out2;
        outbuf[i]  = out1;
    }
    mPhase      = pos;
    mSync       = lastsync;
    mPhaseIn    = phasein;
}

void PulsePTR::next_kkkk(int nSamples) {
    const float freq    = in0(0);
    const float phase   = in0(1);
    const float sync    = in0(2);
    const float width   = in0(3);
    float* outbuf   = out(0);
    float* syncOut  = out(1);

    double pos      = mPhase;
    double rateinv  = 1 / sampleRate();
    double phasein  = mPhaseIn;
    double lastsync = mSync;
    double w        = width < 0.01f ? 0.01f : width > 0.99f ? 0.99f : width;
    double step     = freq * rateinv;
    double step2    = 2.f*step;
    double step3    = 3.f*step;
    double samples  = 1.f / step;
    double p2       = samples * samples / 2.f;
    double p3       = samples * samples * samples / 6.f;
    for (int i = 0; i < nSamples; ++i) {
        double out1     = 0.f;
        double out2     = 0.f;
        if (sync > 0.f && lastsync == 0.f) {
            pos = 0.f;
        }
        else {
            pos = pos + (phase - phasein) + step;
            while (pos >= 1.f) {
                pos -= 1.f;
            }
            while (pos < 0.f) {
                pos += 1.f;
            }
        }
        lastsync    = sync;
        phasein     = phase;
        if (pos < step) out2 = 1.f;
        out1 = 2.f * algorithm(pos, step, step2, step3, w, samples, p2, p3);

        syncOut[i] = out2;
        outbuf[i]  = out1;
    }
    mPhase      = pos;
    mSync       = lastsync;
    mPhaseIn    = phasein;
}
} // namespace PulsePTR

PluginLoad(PulsePTRUGens) {
    // Plugin magic
    ft = inTable;
    registerUnit<PulsePTR::PulsePTR>(ft, "PulsePTR", false);
}
