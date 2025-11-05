# Usage Guide

This guide describes the two main workflows supported by the ODT GPU Library and demonstrates their usage with example function sequences.

---

## 1. Raw-data Workflow

This workflow begins directly from raw holograms (detector data). It performs all preprocessing operations, such as reference correction, FFT windowing, and normalization, on the GPU.

### Example sequence:
```c
HL_addReference(...);                           // optional reference
HL00to02_FromPreprocToGenKO(...);               // hologram -> preprocessing -> K-space
HL03_setParamsAndStartDIandGP(..., nGPi);       // reconstruction (Direct Inverse or iterative)
HL04_takeReconstructionAndFreeMemory(...);      // get final reconstruction
```

**Key features:**
- Input: raw holograms (`short int`)
- Automatic preprocessing on GPU (FFT, filtering, windowing, normalization)
- Optional reference handling (`HL_addReference()`)
- Fastest and simplest path for real-time or automated processing
- No access to intermediate data (sinograms, K-space)

---

## 2. Preprocessed-data Workflow

This workflow begins from preprocessed sinograms (amplitude and phase). It allows precise control over all stages of the reconstruction and access to intermediate data.

### Example sequence:
```c
HL_addReference(...);                           // optional reference
HL01_setParams(...);
HL02_sendDataAndGenerateKO(...);                // sinograms -> K-space
HL03_setParamsAndStartDIandGP(..., nGPi);
HL04_takeReconstructionAndFreeMemory(...);
```

**Key features:**
- Input: preprocessed sinograms (`float`)
- User-controlled preprocessing (external or custom)
- Optional reference handling (`HL_addReference()`)
- Access to intermediate data (e.g., sinograms, K-space)
- Ideal for research, analysis, and debugging

---

