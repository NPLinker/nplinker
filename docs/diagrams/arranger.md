# Dataset Arranging Pipeline

The [DatasetArranger][nplinker.arranger.DatasetArranger] is implemented according to the following flowcharts.

## Strain mappings file
``` mermaid
flowchart TD
    StrainMappings[`strain_mappings.json`] --> SM{Is the mode PODP?}
    SM --> |No |SM0[Validate the file]
    SM --> |Yes|SM1[Generate the file] --> SM0
```

## Strain selection file
``` mermaid
flowchart TD
    StrainsSelected[`strains_selected.json`] --> S{Does the file exist?}
    S --> |No | S0[Nothing to do]
    S --> |Yes| S1[Validate the file]
```

## PODP project metadata json file
``` mermaid
flowchart TD
    podp[PODP project metadata json file] --> A{Is the mode PODP?}
    A --> |No | A0[Nothing to do]
    A --> |Yes| P{Does the file exist?}
    P --> |No | P0[Download the file] --> P1
    P --> |Yes| P1[Validate the file]
```

## GNPS, AntiSMASH and BigScape 
``` mermaid
flowchart TD
    ConfigError[Dynaconf config validation error]
    DataError[Data validation error]
    UseIt[Use the data]
    Download[First remove existing data if relevent, then download or generate data]

    A[GNPS, antiSMASH and BigSCape] --> B{Pass Dynaconf config validation?}
    B -->|No | ConfigError
    B -->|Yes| G{Is the mode PODP?}
    
    G -->|No, local mode| G1{Does data dir exist?}
    G1 -->|No | DataError
    G1 -->|Yes| H{Pass data validation?}
    H --> |No | DataError
    H --> |Yes| UseIt 

    G -->|Yes, podp mode| G2{Does data dir exist?}
    G2 --> |No | Download
    G2 --> |Yes | J{Pass data validation?}
    J -->|No | Download --> |try max 2 times| J
    J -->|Yes| UseIt
```

## MIBiG Data
MIBiG data is always downloaded automatically. Users cannot provide their own MIBiG data.

``` mermaid
flowchart TD
    Mibig[MIBiG] --> M0{Pass Dynaconf config validation?}
    M0 -->|No | M01[Dynaconf config validation error]
    M0 -->|Yes | MibigDownload[First remove existing data if relevant and then download data]
```
