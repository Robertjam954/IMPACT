library(yaml)

outputs <- yaml.load("
est:
  CL:    0.482334
  VC:    0.0592686
  CL_WT: 0.750000
  VC_WT: 1.00000
  nCL:   0.315414
  nVC:   0.536025
  ERRP:  0.0508497
se:
  CL:    0.0138646
  VC:    0.0055512
  nCL:   0.0188891
  nVC:   0.0900352
  ERRP:  0.0018285
fixed:
  CL:    no
  VC:    no
  CL_WT: yes
  VC_WT: yes
  nCL:   no
  nVC:   no
  ERRP:  no
shrinkage:
  nCL:  9.54556
  nVC:  47.8771
")

meta <- yaml.load("
parameters:
- name:  CL
  label: 'Clearance'
  units: 'L/h'
  type:  Structural

- name:  VC
  label: 'Volume'
  units: 'L'
  type:  Structural
  trans: 'exp'
  
- name:  CL_WT
  label: 'Weight on Clearance'
  type:  CovariateEffect

- name:  VC_WT
  label: 'Weight on Volume'
  type:  CovariateEffect
  
- name:  nCL
  label: 'On Clearance'
  type:  IIV
  trans: 'SD (CV%)'
  
- name:  nVC
  label: 'On Volume'
  type:  IIV
  trans: 'SD (CV%)'
  
- name:  ERRP
  label: 'Proportional Error'
  units: '%'
  type:  RUV
  trans: '%'
")

parframe <- pmxparframe(outputs, meta)