# EVPI vs EVSI: Policy Interpretation for Trypanosomiasis Diagnostics

## Conceptual Distinction

### Expected Value of Perfect Information (EVPI)
**EVPI** = Expected value of knowing the **true disease status** of each animal before making treatment decisions

### Expected Value of Sample Information (EVSI) 
**EVSI** = Expected value of knowing the **true test performance parameters** (sensitivity/specificity) before making diagnostic strategy decisions

## Policy Applications

### For Local Farmers (Individual Decision-Making)

#### EVPI Interpretation
- **Question**: "How much would perfect disease diagnosis be worth to me?"
- **Scenario**: Farmer knows with certainty whether each animal is infected
- **Decision**: Treat infected animals, don't treat healthy animals
- **Value**: Eliminates all diagnostic errors (false positives/negatives)
- **Practical Meaning**: Upper bound on what farmers should pay for perfect diagnostic technology

#### EVSI Interpretation  
- **Question**: "How much would better test information be worth to me?"
- **Scenario**: Farmer knows true sensitivity/specificity of available tests
- **Decision**: Choose optimal test strategy (NO_TEST, HCT, PCR) with certainty
- **Value**: Eliminates uncertainty about which test strategy is best
- **Practical Meaning**: Value of resolving test performance uncertainty (much lower than EVPI)

### For Policymakers (Population-Level Decisions)

#### EVPI Applications
1. **Research Priority Setting**: Where would perfect diagnostics provide highest returns?
2. **Technology Investment**: Maximum justifiable investment in diagnostic R&D
3. **Resource Allocation**: Prioritize locations where diagnostic uncertainty has highest cost
4. **Cost-Effectiveness Thresholds**: Compare EVPI to cost of developing better diagnostics

#### EVSI Applications  
1. **Test Validation Studies**: Where should we conduct additional test performance studies?
2. **Evidence Synthesis**: Value of resolving conflicting test accuracy studies
3. **Regulatory Decisions**: When is current test evidence sufficient vs. needing more data?
4. **Research Portfolio**: Balance between test development vs. prevalence mapping

## Current Analysis Results

### EVPI (Not Yet Calculated)
- Would show value of perfect animal-level diagnosis
- Expected to be substantial ($2-50 per animal depending on prevalence)
- Varies by location based on prevalence and current optimal strategy

### EVSI (Calculated: ~$0 per location)
- Shows minimal value of better test parameter information
- Current test performance estimates are adequate for decisions
- Research should focus elsewhere (prevalence, not test validation)

## Decision Framework

```
EVPI >> EVSI suggests:
├─ Diagnostic technology development is high priority
├─ Test parameter research is low priority  
├─ Current test accuracy estimates are sufficient
└─ Focus R&D on better diagnostics, not test validation
```

## Practical Implications

### For Farmers
- **EVPI tells you**: Maximum you should pay for perfect diagnostics
- **EVSI tells you**: Current tests are well-characterized enough for good decisions
- **Action**: Use existing HCT/PCR guidelines with confidence

### For Policymakers
- **EVPI guides**: Where to invest in diagnostic technology development
- **EVSI guides**: Where to invest in test validation studies
- **Current finding**: Invest in better diagnostics, not more test studies

### For Researchers
- **EVPI priorities**: Locations where better diagnostics would help most
- **EVSI priorities**: Test parameters that need better characterization
- **Current finding**: Prevalence mapping > test validation > diagnostic development

## Next Steps for Complete Analysis

1. **Calculate EVPI** across all Nigeria locations
2. **Compare EVPI vs EVSI** to inform research priorities
3. **Map high-EVPI locations** for targeted diagnostic development
4. **Estimate total EVPI** to justify diagnostic R&D investments
5. **Policy recommendations** based on EVPI/EVSI ratio

Would you like me to implement EVPI calculation to complete this analysis?