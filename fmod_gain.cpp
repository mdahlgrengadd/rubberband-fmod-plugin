/*==============================================================================
Gain DSP Plugin Example
Copyright (c), Firelight Technologies Pty, Ltd 2004-2021.

This example shows how to create a simple gain DSP effect.
==============================================================================*/

#ifdef WIN32
    #define _CRT_SECURE_NO_WARNINGS
#endif


#include <math.h>
#include <stdio.h>
#include <string.h>
#include <memory>
#include <algorithm>
#include "fmod.hpp"
#include "rubberband\RubberBandStretcher.h"

#include "src\base\RingBuffer.h"

#define FMOD_GAIN_USEPROCESSCALLBACK            /* FMOD plugins have 2 methods of processing data.  
                                                    1. via a 'read' callback which is compatible with FMOD Ex but limited in functionality, or 
                                                    2. via a 'process' callback which exposes more functionality, like masks and query before process early out logic. */

extern "C" {
    F_EXPORT FMOD_DSP_DESCRIPTION* F_CALL FMODGetDSPDescription();
}

const double FMOD_GAIN_PARAM_GAIN_MIN     = -24.0;
const double FMOD_GAIN_PARAM_GAIN_MAX     = 24.0;
const double FMOD_GAIN_PARAM_GAIN_DEFAULT = 0.0;
#define FMOD_GAIN_RAMPCOUNT 1

enum
{
    FMOD_GAIN_PARAM_GAIN = 0,
    FMOD_GAIN_PARAM_INVERT,
    FMOD_GAIN_NUM_PARAMETERS
};

//#define DECIBELS_TO_LINEAR(__dbval__)  ((__dbval__ <= FMOD_GAIN_PARAM_GAIN_MIN) ? 0.0f : powf(10.0f, __dbval__ / 20.0f))
//#define LINEAR_TO_DECIBELS(__linval__) ((__linval__ <= 0.0f) ? FMOD_GAIN_PARAM_GAIN_MIN : 20.0f * log10f((float)__linval__))
#define OCTAVE_TO_PITCH_RATIO(__octval__) 1.0 / pow(2.0, __octval__ / 12.0)

FMOD_RESULT F_CALLBACK FMOD_Gain_dspcreate       (FMOD_DSP_STATE *dsp_state);
FMOD_RESULT F_CALLBACK FMOD_Gain_dsprelease      (FMOD_DSP_STATE *dsp_state);
FMOD_RESULT F_CALLBACK FMOD_Gain_dspreset        (FMOD_DSP_STATE *dsp_state);
#ifdef FMOD_GAIN_USEPROCESSCALLBACK
FMOD_RESULT F_CALLBACK FMOD_Gain_dspprocess      (FMOD_DSP_STATE *dsp_state, unsigned int length, const FMOD_DSP_BUFFER_ARRAY *inbufferarray, FMOD_DSP_BUFFER_ARRAY *outbufferarray, FMOD_BOOL inputsidle, FMOD_DSP_PROCESS_OPERATION op);
#else
FMOD_RESULT F_CALLBACK FMOD_Gain_dspread         (FMOD_DSP_STATE *dsp_state, float *inbuffer, float *outbuffer, unsigned int length, int inchannels, int *outchannels);
#endif
FMOD_RESULT F_CALLBACK FMOD_Gain_dspsetparamfloat(FMOD_DSP_STATE *dsp_state, int index, float value);
FMOD_RESULT F_CALLBACK FMOD_Gain_dspsetparamint  (FMOD_DSP_STATE *dsp_state, int index, int value);
FMOD_RESULT F_CALLBACK FMOD_Gain_dspsetparambool (FMOD_DSP_STATE *dsp_state, int index, FMOD_BOOL value);
FMOD_RESULT F_CALLBACK FMOD_Gain_dspsetparamdata (FMOD_DSP_STATE *dsp_state, int index, void *data, unsigned int length);
FMOD_RESULT F_CALLBACK FMOD_Gain_dspgetparamfloat(FMOD_DSP_STATE *dsp_state, int index, float *value, char *valuestr);
FMOD_RESULT F_CALLBACK FMOD_Gain_dspgetparamint  (FMOD_DSP_STATE *dsp_state, int index, int *value, char *valuestr);
FMOD_RESULT F_CALLBACK FMOD_Gain_dspgetparambool (FMOD_DSP_STATE *dsp_state, int index, FMOD_BOOL *value, char *valuestr);
FMOD_RESULT F_CALLBACK FMOD_Gain_dspgetparamdata (FMOD_DSP_STATE *dsp_state, int index, void **value, unsigned int *length, char *valuestr);
FMOD_RESULT F_CALLBACK FMOD_Gain_shouldiprocess  (FMOD_DSP_STATE *dsp_state, FMOD_BOOL inputsidle, unsigned int length, FMOD_CHANNELMASK inmask, int inchannels, FMOD_SPEAKERMODE speakermode);
FMOD_RESULT F_CALLBACK FMOD_Gain_sys_register    (FMOD_DSP_STATE *dsp_state);
FMOD_RESULT F_CALLBACK FMOD_Gain_sys_deregister  (FMOD_DSP_STATE *dsp_state);
FMOD_RESULT F_CALLBACK FMOD_Gain_sys_mix         (FMOD_DSP_STATE *dsp_state, int stage);

static bool                    FMOD_Gain_Running = false;
static FMOD_DSP_PARAMETER_DESC p_gain;
static FMOD_DSP_PARAMETER_DESC p_invert;
    
FMOD_DSP_PARAMETER_DESC *FMOD_Gain_dspparam[FMOD_GAIN_NUM_PARAMETERS] =
{
    &p_gain,
    &p_invert
};

FMOD_DSP_DESCRIPTION FMOD_Gain_Desc =
{
    FMOD_PLUGIN_SDK_VERSION,
    "Super Gain",    // name
    0x00010000,     // plug-in version
    1,              // number of input buffers to process
    1,              // number of output buffers to process
    FMOD_Gain_dspcreate,
    FMOD_Gain_dsprelease,
    FMOD_Gain_dspreset,
#ifndef FMOD_GAIN_USEPROCESSCALLBACK
    FMOD_Gain_dspread,
#else
    0,
#endif
#ifdef FMOD_GAIN_USEPROCESSCALLBACK
    FMOD_Gain_dspprocess,
#else
    0,
#endif
    0,
    FMOD_GAIN_NUM_PARAMETERS,
    FMOD_Gain_dspparam,
    FMOD_Gain_dspsetparamfloat,
    0, // FMOD_Gain_dspsetparamint,
    FMOD_Gain_dspsetparambool,
    0, // FMOD_Gain_dspsetparamdata,
    FMOD_Gain_dspgetparamfloat,
    0, // FMOD_Gain_dspgetparamint,
    FMOD_Gain_dspgetparambool,
    0, // FMOD_Gain_dspgetparamdata,
    FMOD_Gain_shouldiprocess,
    0,                                      // userdata
    FMOD_Gain_sys_register,
    FMOD_Gain_sys_deregister,
    FMOD_Gain_sys_mix
};

extern "C"
{

F_EXPORT FMOD_DSP_DESCRIPTION* F_CALL FMODGetDSPDescription()
{
    FMOD_DSP_INIT_PARAMDESC_FLOAT(p_gain, "Pitch", "", "Pitch in ratio. 0.5 to 2.0 Default = 1", FMOD_GAIN_PARAM_GAIN_MIN, FMOD_GAIN_PARAM_GAIN_MAX, FMOD_GAIN_PARAM_GAIN_DEFAULT);
    FMOD_DSP_INIT_PARAMDESC_BOOL(p_invert, "Invert", "", "Invert signal. Default = off", false, 0);
    return &FMOD_Gain_Desc;
}

}

class FMODGainState
{
public:
    FMODGainState();

    void read(float *inbuffer, float *outbuffer, unsigned int length, int channels);
    void reset();
    void setGain(float);
    void setInvert(bool);
    float gain() const { return m_target_gain; } //LINEAR_TO_DECIBELS(m_invert ? -m_target_gain : m_target_gain); }
    FMOD_BOOL invert() const { return m_invert; }

    std::unique_ptr<RubberBand::RubberBandStretcher> m_stretcher_stereo;
    std::unique_ptr<RubberBand::RubberBandStretcher> m_stretcher_mono;
    std::unique_ptr<RubberBand::RubberBandStretcher>* m_stretcher;
    float tempbuf[2][8192*2];
    float tempbuf2[2][8192 * 2];

    RubberBand::RingBuffer<float>** m_outputBuffer;
    const float** m_input;
    float** m_output;
    float* m_latency;
    float* m_cents;
    float* m_semitones;
    float* m_octaves;
    float* m_crispness;
    float* m_formant;
    float* m_fast;
    double m_ratio;
    double m_prevRatio;
    int m_currentCrispness;
    bool m_currentFormant;
    bool m_currentFast;

    size_t m_blockSize;
    size_t m_reserve;
    size_t m_minfill;

    float** m_scratch;

    int m_sampleRate;
    size_t m_channels;

private:
    double m_target_gain;
    double m_current_gain;
    int   m_ramp_samples_left;
    bool  m_invert;
};

FMODGainState::FMODGainState()
{
    m_target_gain = FMOD_GAIN_PARAM_GAIN_DEFAULT;// DECIBELS_TO_LINEAR(FMOD_GAIN_PARAM_GAIN_DEFAULT);
    m_invert = 0;
    reset();
}

void FMODGainState::read(float *inbuffer, float *outbuffer, unsigned int length, int channels)
{

    // Note: buffers are interleaved


    switch (channels) {
    case 2:
        m_stretcher = &m_stretcher_stereo;
        break;
    case 1:
        m_stretcher = &m_stretcher_mono;
        break;

    default: // Bail out if anything else than mono or stereo
        return;
    }

    if (m_current_gain != m_target_gain) {
        m_current_gain = m_target_gain;
        double stretchRatio = OCTAVE_TO_PITCH_RATIO(m_current_gain);
        m_stretcher_mono->setPitchScale(stretchRatio);
        m_stretcher_stereo->setPitchScale(stretchRatio);
    }

   
    // De-interleave
    for (int i = 0; i < length; i++) {
        for (int c = 0; c < channels; c++) {
            tempbuf[c][i] = inbuffer[c + i * channels];
        }
    }
    
    const float* ptrs[2];

    for (size_t c = 0; c < channels; ++c) {
        ptrs[c] = &(tempbuf[c][0]);
    }

    //From Ladspa
    const int samples = length;
    int processed = 0;
    size_t outTotal = 0;
    long offset = 0;

    int rs = m_outputBuffer[0]->getReadSpace();
    if (rs < int(m_minfill)) {
        //        cerr << "temporary expansion (have " << rs << ", want " << m_reserve << ")" << endl;
        (*m_stretcher)->setTimeRatio(1.1); // fill up temporarily
    }
    else if (rs > 8192) {
        //        cerr << "temporary reduction (have " << rs << ", want " << m_reserve << ")" << endl;
        (*m_stretcher)->setTimeRatio(0.9); // reduce temporarily
    }
    else {
        (*m_stretcher)->setTimeRatio(1.0);
    }

    while (processed < samples) {

        // never feed more than the minimum necessary number of
        // samples at a time; ensures nothing will overflow internally
        // and we don't need to call setMaxProcessSize

        int toCauseProcessing = (*m_stretcher)->getSamplesRequired();
        int inchunk = std::min((samples -processed), toCauseProcessing);
        for (size_t c = 0; c < channels; ++c) {
            ptrs[c] = &(tempbuf[c][offset + processed]);
        }
        (*m_stretcher)->process(ptrs, inchunk, false);
        processed += inchunk;

        int avail = (*m_stretcher)->available();
        int writable = m_outputBuffer[0]->getWriteSpace();
        int outchunk = std::min(int(avail), writable);
        size_t actual = (*m_stretcher)->retrieve(m_scratch, outchunk);
        outTotal += actual;

        //        cout << "avail: " << avail << ", outchunk = " << outchunk;
        //        if (actual != outchunk) cout << " (" << actual << ")";
        //        cout << endl;

        outchunk = actual;

        for (size_t c = 0; c < channels; ++c) {
            if (int(m_outputBuffer[c]->getWriteSpace()) < outchunk) {
                //std::cerr << "RubberBandPitchShifter::runImpl: buffer overrun: chunk = " << outchunk << ", space = " << m_outputBuffer[c]->getWriteSpace() << std::endl;
            }
            m_outputBuffer[c]->write(m_scratch[c], outchunk);
        }
    }

    for (size_t c = 0; c < channels; ++c) {
        int toRead = m_outputBuffer[c]->getReadSpace();
        if (toRead < samples && c == 0) {
            //std::cerr << "RubberBandPitchShifter::runImpl: buffer underrun: required = " << samples << ", available = " << toRead << std::endl;
        }
        int chunk = std::min(int(toRead), samples);
        m_outputBuffer[c]->read(&(tempbuf2[c][offset]), chunk);

        // Re-interleave
        for (size_t i = 0; i < chunk; i++) {
                outbuffer[c + i * (size_t)channels] = tempbuf2[c][i];
        }
    }

    if (m_minfill == 0) {
        m_minfill = m_outputBuffer[0]->getReadSpace();
        //        cerr << "minfill = " << m_minfill << endl;
    }

    //for (int i = 0; i < length; i++) {
    //    for (int c = 0; c < channels; c++) {
    //    
    //        outbuffer[c  + i * channels] = tempbuf2[c][i];


    //    }
    //}


}

void FMODGainState::reset()
{
    m_current_gain = m_target_gain;
    m_ramp_samples_left = 0;
}

void FMODGainState::setGain(float gain)
{
    m_target_gain = gain;// m_invert ? -DECIBELS_TO_LINEAR(gain) : DECIBELS_TO_LINEAR(gain);
    m_ramp_samples_left = FMOD_GAIN_RAMPCOUNT;
}

void FMODGainState::setInvert(bool invert)
{
    //if (invert != m_invert)
    //{
        //m_target_gain = -m_target_gain;
        //m_ramp_samples_left = FMOD_GAIN_RAMPCOUNT;
    //}
    m_invert = invert;
}

FMOD_RESULT F_CALLBACK FMOD_Gain_dspcreate(FMOD_DSP_STATE *dsp_state)
{
    //Some defaults... will be overriden below
    int channels = 2;
    int sampleRate = 44100;
    unsigned int m_blockSize = 1024;

    dsp_state->plugindata = (FMODGainState *)FMOD_DSP_ALLOC(dsp_state, sizeof(FMODGainState));
    if (!dsp_state->plugindata)
    {
        return FMOD_ERR_MEMORY;
    }

    dsp_state->functions->getsamplerate(dsp_state, &sampleRate);
    dsp_state->functions->getblocksize(dsp_state, &m_blockSize);
    
    FMODGainState* state = (FMODGainState*)dsp_state->plugindata;

    // I first tried to initialize a rubberbandstretcher based on how many channels
    // the dsp had. This led to crashes when changing channels at runtime.
    // Instead I know create 2 stretchers, one for mono and one for stereo
    //FMOD_SPEAKERMODE* numInputChannels = (FMOD_SPEAKERMODE*)FMOD_DSP_ALLOC(dsp_state, sizeof(FMOD_SPEAKERMODE));
    //FMOD_SPEAKERMODE* numOutputChannels = (FMOD_SPEAKERMODE*)FMOD_DSP_ALLOC(dsp_state, sizeof(FMOD_SPEAKERMODE));
    //dsp_state->functions->getspeakermode(dsp_state, numInputChannels, numOutputChannels);
    state->m_stretcher_stereo = std::make_unique<RubberBand::RubberBandStretcher>(sampleRate,
        2,
        RubberBand::RubberBandStretcher::Option::OptionProcessRealTime
        || RubberBand::RubberBandStretcher::Option::OptionFormantPreserved
        || RubberBand::RubberBandStretcher::Option::OptionStretchElastic
        || RubberBand::RubberBandStretcher::Option::OptionPitchHighQuality
        || RubberBand::RubberBandStretcher::Option::OptionTransientsSmooth
        || RubberBand::RubberBandStretcher::Option::OptionSmoothingOn
        || RubberBand::RubberBandStretcher::Option::OptionWindowLong 
        || RubberBand::RubberBandStretcher::Option::OptionPhaseLaminar
        || RubberBand::RubberBandStretcher::Option::OptionChannelsApart
        || RubberBand::RubberBandStretcher::Option::OptionDetectorSoft,
        1.0,
        1.0);

    state->m_stretcher_mono = std::make_unique<RubberBand::RubberBandStretcher>(sampleRate,
        1,
        RubberBand::RubberBandStretcher::Option::OptionProcessRealTime
        || RubberBand::RubberBandStretcher::Option::OptionFormantPreserved
        || RubberBand::RubberBandStretcher::Option::OptionStretchElastic
        || RubberBand::RubberBandStretcher::Option::OptionPitchHighQuality
        || RubberBand::RubberBandStretcher::Option::OptionTransientsSmooth
        || RubberBand::RubberBandStretcher::Option::OptionSmoothingOn
        || RubberBand::RubberBandStretcher::Option::OptionWindowLong
        || RubberBand::RubberBandStretcher::Option::OptionPhaseLaminar
        || RubberBand::RubberBandStretcher::Option::OptionChannelsApart
        || RubberBand::RubberBandStretcher::Option::OptionDetectorSoft,
        1.0,
        1.0);



    state->m_currentFast = false;
    state->m_blockSize = m_blockSize;
    state->m_reserve = 1024;
    state->m_minfill = 0;
    state->m_sampleRate = sampleRate;
    state->m_outputBuffer = new RubberBand::RingBuffer<float> * [channels];
    state->m_scratch = new float* [channels];
    // Some more values that can be set:
    //m_latency(0),
    //m_cents(0),
    //m_semitones(0),
    //m_octaves(0),
    //m_crispness(0),
    //m_formant(0),
    //m_fast(0),
    //m_ratio(1.0),
    //m_prevRatio(1.0),
    //m_currentCrispness(-1),
    //m_currentFormant = false;

    for (size_t c = 0; c < channels; ++c) {

        int bufsize = m_blockSize + state->m_reserve + 8192;

        state->m_outputBuffer[c] = new RubberBand::RingBuffer<float>(bufsize);

        state->m_scratch[c] = new float[bufsize];
        for (int i = 0; i < bufsize; ++i) state->m_scratch[c][i] = 0.f;
    }

    return FMOD_OK;
}

FMOD_RESULT F_CALLBACK FMOD_Gain_dsprelease(FMOD_DSP_STATE *dsp_state)
{
    FMODGainState *state = (FMODGainState *)dsp_state->plugindata;
    FMOD_DSP_FREE(dsp_state, state);
    return FMOD_OK;
}

#ifdef FMOD_GAIN_USEPROCESSCALLBACK

FMOD_RESULT F_CALLBACK FMOD_Gain_dspprocess(FMOD_DSP_STATE *dsp_state, unsigned int length, const FMOD_DSP_BUFFER_ARRAY *inbufferarray, FMOD_DSP_BUFFER_ARRAY *outbufferarray, FMOD_BOOL inputsidle, FMOD_DSP_PROCESS_OPERATION op)
{
    FMODGainState *state = (FMODGainState *)dsp_state->plugindata;

    if (op == FMOD_DSP_PROCESS_QUERY)
    {
        if (outbufferarray && inbufferarray)
        {
            outbufferarray[0].buffernumchannels[0] = inbufferarray[0].buffernumchannels[0];
            outbufferarray[0].speakermode       = inbufferarray[0].speakermode;
        }

        if (inputsidle)
        {
            return FMOD_ERR_DSP_DONTPROCESS;
        }
    }
    else
    {
        state->read(inbufferarray[0].buffers[0], outbufferarray[0].buffers[0], length, inbufferarray[0].buffernumchannels[0]); // input and output channels count match for this effect
    }

    return FMOD_OK;
}

#else

FMOD_RESULT F_CALLBACK FMOD_Gain_dspread(FMOD_DSP_STATE *dsp_state, float *inbuffer, float *outbuffer, unsigned int length, int inchannels, int * /*outchannels*/)
{
    FMODGainState *state = (FMODGainState *)dsp_state->plugindata;
    state->read(inbuffer, outbuffer, length, inchannels); // input and output channels count match for this effect
    return FMOD_OK;
}

#endif

FMOD_RESULT F_CALLBACK FMOD_Gain_dspreset(FMOD_DSP_STATE *dsp_state)
{
    FMODGainState *state = (FMODGainState *)dsp_state->plugindata;
    state->reset();
    return FMOD_OK;
}

FMOD_RESULT F_CALLBACK FMOD_Gain_dspsetparamfloat(FMOD_DSP_STATE *dsp_state, int index, float value)
{
    FMODGainState *state = (FMODGainState *)dsp_state->plugindata;

    switch (index)
    {
    case FMOD_GAIN_PARAM_GAIN:
        state->setGain(value);
        return FMOD_OK;
    }

    return FMOD_ERR_INVALID_PARAM;
}

FMOD_RESULT F_CALLBACK FMOD_Gain_dspgetparamfloat(FMOD_DSP_STATE *dsp_state, int index, float *value, char *valuestr)
{
    FMODGainState *state = (FMODGainState *)dsp_state->plugindata;

    switch (index)
    {
    case FMOD_GAIN_PARAM_GAIN:
        *value = state->gain();
        if (valuestr) sprintf(valuestr, "%.1f dB", state->gain());
        return FMOD_OK;
    }

    return FMOD_ERR_INVALID_PARAM;
}

FMOD_RESULT F_CALLBACK FMOD_Gain_dspsetparambool(FMOD_DSP_STATE *dsp_state, int index, FMOD_BOOL value)
{
    FMODGainState *state = (FMODGainState *)dsp_state->plugindata;

    switch (index)
    {
      case FMOD_GAIN_PARAM_INVERT:
        state->setInvert(value ? true : false);
        return FMOD_OK;
    }

    return FMOD_ERR_INVALID_PARAM;
}

FMOD_RESULT F_CALLBACK FMOD_Gain_dspgetparambool(FMOD_DSP_STATE *dsp_state, int index, FMOD_BOOL *value, char *valuestr)
{
    FMODGainState *state = (FMODGainState *)dsp_state->plugindata;

    switch (index)
    {
    case FMOD_GAIN_PARAM_INVERT:
        *value = state->invert();
        if (valuestr) sprintf(valuestr, state->invert() ? "Inverted" : "Off" );
        return FMOD_OK;
    }

    return FMOD_ERR_INVALID_PARAM;
}

FMOD_RESULT F_CALLBACK FMOD_Gain_shouldiprocess(FMOD_DSP_STATE * /*dsp_state*/, FMOD_BOOL inputsidle, unsigned int /*length*/, FMOD_CHANNELMASK /*inmask*/, int /*inchannels*/, FMOD_SPEAKERMODE /*speakermode*/)
{
    if (inputsidle)
    {
        return FMOD_ERR_DSP_DONTPROCESS;
    }

    return FMOD_OK;
}


FMOD_RESULT F_CALLBACK FMOD_Gain_sys_register(FMOD_DSP_STATE * /*dsp_state*/)
{
    FMOD_Gain_Running = true;
    // called once for this type of dsp being loaded or registered (it is not per instance)
    return FMOD_OK;
}

FMOD_RESULT F_CALLBACK FMOD_Gain_sys_deregister(FMOD_DSP_STATE * /*dsp_state*/)
{
    FMOD_Gain_Running = false;
    // called once for this type of dsp being unloaded or de-registered (it is not per instance)
    return FMOD_OK;
}

FMOD_RESULT F_CALLBACK FMOD_Gain_sys_mix(FMOD_DSP_STATE * /*dsp_state*/, int /*stage*/)
{
    // stage == 0 , before all dsps are processed/mixed, this callback is called once for this type.
    // stage == 1 , after all dsps are processed/mixed, this callback is called once for this type.
    return FMOD_OK;
}
